// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_foam.hh>
#include <dumux/adaptive/markelements.hh>
#include <dumux/porousmediumflow/velocityoutput.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>

#include "properties.hh"
#include "particles.hh"

int main(int argc, char** argv) try
{
    using namespace Dumux;

    // maybe initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc, argv);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // Define the sub problem type tags
    using BulkTypeTag = Properties::TTag::BULKTYPETAG;
    using LowDimTypeTag = Properties::TTag::LOWDIMTYPETAG;

    // create grids
    using BulkGridManager = Dumux::GridManager<GetPropType<BulkTypeTag, Properties::Grid>>;
    BulkGridManager bulkGridManager;
    bulkGridManager.init("Tissue"); // pass parameter group

    using LowDimGridManager = Dumux::GridManager<GetPropType<LowDimTypeTag, Properties::Grid>>;
    LowDimGridManager lowDimGridManager;
    lowDimGridManager.init("Vessel"); // pass parameter group
    auto lowDimGridData = lowDimGridManager.getGridData();

    // we compute on the leaf grid view
    const auto& bulkGridView = bulkGridManager.grid().leafGridView();
    const auto& lowDimGridView = lowDimGridManager.grid().leafGridView();

    // refine grid elements large than a threshold
    {
        auto& lowDimGrid = lowDimGridManager.grid();
        const auto maxElementLength = getParam<double>("Vessel.Grid.MaxElementLength", 1e-6);
        std::cout << "Refining grid with " << lowDimGridView.size(0)
                  << " elements until max length: " << maxElementLength << std::endl;

        bool markedElements = true;
        while (markedElements)
        {
            markedElements = markElements(lowDimGrid,[maxElementLength](const auto& element){
                return element.geometry().volume() > maxElementLength ? 1 : 0;
            }, false);

            if (markedElements)
            {
                lowDimGrid.preAdapt();
                lowDimGrid.adapt();
                lowDimGrid.postAdapt();
            }
        }

        auto maxLength = 0.0;
        auto minLength = 1.0;
        for (const auto& element : elements(lowDimGridView))
        {
            const auto length = element.geometry().volume();
            maxLength = std::max(maxLength, length);
            minLength = std::min(minLength, length);
        }

        std::cout << "After refinement: " << lowDimGridView.size(0) << " elements:" << std::endl
                  << "  -- max length: " << maxLength << std::endl
                  << "  -- min length: " << minLength << std::endl
                  << "  -- max level: " << lowDimGridView.grid().maxLevel() << std::endl;
    }

    // create the finite volume grid geometry
    using BulkGridGeometry = GetPropType<BulkTypeTag, Properties::GridGeometry>;
    auto bulkGridGeometry = std::make_shared<BulkGridGeometry>(bulkGridView);
    bulkGridGeometry->update();
    using LowDimGridGeometry = GetPropType<LowDimTypeTag, Properties::GridGeometry>;
    auto lowDimGridGeometry = std::make_shared<LowDimGridGeometry>(lowDimGridView);
    lowDimGridGeometry->update();

    // the mixed dimension type traits
    using Traits = MultiDomainTraits<BulkTypeTag, LowDimTypeTag>;
    constexpr auto bulkIdx = Traits::template SubDomain<0>::Index();
    constexpr auto lowDimIdx = Traits::template SubDomain<1>::Index();

    // the coupling manager
    using CouplingManager = GetPropType<BulkTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>(bulkGridGeometry, lowDimGridGeometry);

    // the problem (initial and boundary conditions)
    using BulkProblem = GetPropType<BulkTypeTag, Properties::Problem>;
    auto bulkProblem = std::make_shared<BulkProblem>(bulkGridGeometry, couplingManager);
    using LowDimProblem = GetPropType<LowDimTypeTag, Properties::Problem>;
    auto lowDimProblem = std::make_shared<LowDimProblem>(lowDimGridGeometry, couplingManager, lowDimGridData);

    // the solution vector
    Traits::SolutionVector sol;
    sol[bulkIdx].resize(bulkGridGeometry->numDofs());
    sol[lowDimIdx].resize(lowDimGridGeometry->numDofs());
    bulkProblem->applyInitialSolution(sol[bulkIdx]);
    lowDimProblem->applyInitialSolution(sol[lowDimIdx]);
    auto oldSol = sol;

    couplingManager->init(bulkProblem, lowDimProblem, sol);
    bulkProblem->computePointSourceMap();
    lowDimProblem->computePointSourceMap();

    // the grid variables
    using BulkGridVariables = GetPropType<BulkTypeTag, Properties::GridVariables>;
    auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkGridGeometry);
    bulkGridVariables->init(sol[bulkIdx]);
    using LowDimGridVariables = GetPropType<LowDimTypeTag, Properties::GridVariables>;
    auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimGridGeometry);
    lowDimGridVariables->init(sol[lowDimIdx]);

    // intialize the vtk output module
    using BulkSolutionVector = std::decay_t<decltype(sol[bulkIdx])>;
    VtkOutputModule<BulkGridVariables, BulkSolutionVector> bulkVtkWriter(*bulkGridVariables, sol[bulkIdx], bulkProblem->name());
    GetPropType<BulkTypeTag, Properties::IOFields>::initOutputModule(bulkVtkWriter);
    auto bulkVelocityOutput = std::make_shared<GetPropType<BulkTypeTag, Properties::VelocityOutput>>(*bulkGridVariables);
    bulkVtkWriter.addVelocityOutput(bulkVelocityOutput);
    bulkVtkWriter.write(0.0);

    using LowDimSolutionVector = std::decay_t<decltype(sol[lowDimIdx])>;
    VtkOutputModule<LowDimGridVariables, LowDimSolutionVector> lowDimVtkWriter(*lowDimGridVariables, sol[lowDimIdx], lowDimProblem->name());
    GetPropType<LowDimTypeTag, Properties::IOFields>::initOutputModule(lowDimVtkWriter);
    auto lowDimVelocityOutput = std::make_shared<GetPropType<LowDimTypeTag, Properties::VelocityOutput>>(*lowDimGridVariables);
    lowDimVtkWriter.addVelocityOutput(lowDimVelocityOutput);
    lowDimVtkWriter.addField(lowDimProblem->spatialParams().getRadii(), "radius");
    lowDimVtkWriter.write(0.0);

    // the assembler with time loop for instationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(bulkProblem, lowDimProblem),
                                                 std::make_tuple(bulkGridGeometry, lowDimGridGeometry),
                                                 std::make_tuple(bulkGridVariables, lowDimGridVariables),
                                                 couplingManager);

    // the linear solver
    using LinearSolver = BlockDiagILU0BiCGSTABSolver;
    auto linearSolver = std::make_shared<LinearSolver>();

    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);
    nonLinearSolver.solve(sol);
    bulkVtkWriter.write(1.0);
    lowDimVtkWriter.write(1.0);

    // compute velocity fields
    using BulkVelocityBackend = PorousMediumFlowVelocity<BulkGridVariables, GetPropType<BulkTypeTag, Properties::FluxVariables>>;
    BulkVelocityBackend bulkVelocityComputer(*bulkGridVariables);
    using LowDimVelocityBackend = PorousMediumFlowVelocity<LowDimGridVariables, GetPropType<LowDimTypeTag, Properties::FluxVariables>>;
    LowDimVelocityBackend lowDimVelocityComputer(*lowDimGridVariables);
    BulkVelocityBackend::VelocityVector bulkVelocity(bulkGridView.size(0));
    LowDimVelocityBackend::VelocityVector lowDimVelocity(lowDimGridView.size(0));

    for (const auto& element : elements(bulkGridView))
    {
        auto fvGeometry = localView(*bulkGridGeometry); fvGeometry.bind(element);
        auto elemVolVars = localView(bulkGridVariables->curGridVolVars()); elemVolVars.bind(element, fvGeometry, sol[bulkIdx]);
        auto elemFluxVarsCache = localView(bulkGridVariables->gridFluxVarsCache()); elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);
        bulkVelocityComputer.calculateVelocity(bulkVelocity, element, fvGeometry, elemVolVars, elemFluxVarsCache, 0);
    }

    for (const auto& element : elements(lowDimGridView))
    {
        auto fvGeometry = localView(*lowDimGridGeometry); fvGeometry.bind(element);
        auto elemVolVars = localView(lowDimGridVariables->curGridVolVars()); elemVolVars.bind(element, fvGeometry, sol[lowDimIdx]);
        auto elemFluxVarsCache = localView(lowDimGridVariables->gridFluxVarsCache()); elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);
        lowDimVelocityComputer.calculateVelocity(lowDimVelocity, element, fvGeometry, elemVolVars, elemFluxVarsCache, 0);
    }

    const auto extravasationProbability = lowDimProblem->extravasationProbability(lowDimVelocity, sol[lowDimIdx], *lowDimGridVariables);

    // particle simulation
    const auto tEnd = getParam<double>("Particles.EndTime");
    const auto dt = getParam<double>("Particles.TimeStepSize");
    auto timeLoop = std::make_shared<TimeLoop<double>>(0.0, dt, tEnd);

    MicrovesselParticleAlgorithm<LowDimGridGeometry, BulkGridGeometry> particleAlgorithm(
        lowDimGridGeometry, bulkGridGeometry, lowDimVelocity, bulkVelocity, extravasationProbability, lowDimProblem->spatialParams().getRadii()
    );
    ParticlePVDWriter particleWriter(particleAlgorithm.cloud(), "test_particles");

    timeLoop->start();
    while (!timeLoop->finished())
    {
        particleAlgorithm.run(timeLoop->timeStepSize());
        std::cout << "Number of particles: " << particleAlgorithm.cloud()->size() << std::endl;
        timeLoop->advanceTimeStep();
        particleWriter.write(timeLoop->time());
        timeLoop->reportTimeStep();
    }
    timeLoop->finalize();

    // print parameter report
    Parameters::print();

    return 0;
}
/////////////////////////////////////////////
///////////// Error handler /////////////////
/////////////////////////////////////////////
catch (const Dumux::ParameterException& e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (const Dune::DGFException& e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (const Dune::Exception& e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
