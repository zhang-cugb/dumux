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
/*!
 * \file
 * \ingroup RANSTests
 * \brief Pipe flow test for the staggered grid RANS model,
 *
 * This test simulates is based on pipe flow experiments by
 * John Laufers experiments in 1954 \cite Laufer1954a.
 */
#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include <dumux/assembly/staggeredfvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/staggeredvtkoutputmodule.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/discretization/method.hh>

#include "problem.hh"

/*!
 * \brief Provides an interface for customizing error messages associated with
 *        reading in parameters.
 *
 * \param progName  The name of the program, that was tried to be started.
 * \param errorMsg  The error message that was issued by the start function.
 *                  Comprises the thing that went wrong and a general help message.
 */
void usage(const char *progName, const std::string &errorMsg)
{
    if (errorMsg.size() > 0) {
        std::string errorMessageOut = "\nUsage: ";
                    errorMessageOut += progName;
                    errorMessageOut += " [options]\n";
                    errorMessageOut += errorMsg;
                    errorMessageOut += "\nPlease use the provided input files.\n";
        std::cout << errorMessageOut
                  << "\n";
    }
}
/*!
* \brief Creates analytical solution.
* Returns a tuple of the analytical solution for the pressure, the velocity and the velocity at the faces
* \param time the time at which to evaluate the analytical solution
* \param problem the problem for which to evaluate the analytical solution
*/
template<class Scalar, class Problem>
auto createAnalyticalSolution(const Scalar time, const Problem& problem)
{
    const auto& gridGeometry = problem.gridGeometry();
    using GridView = typename std::decay_t<decltype(gridGeometry)>::GridView;

    static constexpr auto dim = GridView::dimension;
    static constexpr auto dimWorld = GridView::dimensionworld;

    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    std::vector<Scalar> analyticalPressure;
    std::vector<Scalar> analyticalTurbulentKineticEnergy;
    std::vector<Scalar> analyticalDissipation;
    std::vector<VelocityVector> analyticalVelocity;
    std::vector<VelocityVector> analyticalVelocityOnFace;

    analyticalPressure.resize(gridGeometry.numCellCenterDofs());
    analyticalTurbulentKineticEnergy.resize(gridGeometry.numCellCenterDofs());
    analyticalDissipation.resize(gridGeometry.numCellCenterDofs());
    analyticalVelocity.resize(gridGeometry.numCellCenterDofs());
    analyticalVelocityOnFace.resize(gridGeometry.numFaceDofs());

    using Indices = typename Problem::Indices;
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        auto fvGeometry = localView(gridGeometry);
        fvGeometry.bindElement(element);
        for (auto&& scv : scvs(fvGeometry))
        {
            auto ccDofIdx = scv.dofIndex();
            auto ccDofPosition = scv.dofPosition();
            auto analyticalSolutionAtCc = problem.analyticalSolution(ccDofPosition, time);

            // velocities on faces
            for (auto&& scvf : scvfs(fvGeometry))
            {
                const auto faceDofIdx = scvf.dofIndex();
                const auto faceDofPosition = scvf.center();
                const auto dirIdx = scvf.directionIndex();
                const auto analyticalSolutionAtFace = problem.analyticalSolution(faceDofPosition, time);
                analyticalVelocityOnFace[faceDofIdx][dirIdx] = analyticalSolutionAtFace[Indices::velocity(dirIdx)];
            }

            analyticalPressure[ccDofIdx] = analyticalSolutionAtCc[Indices::pressureIdx];
            analyticalTurbulentKineticEnergy[ccDofIdx] = analyticalSolutionAtCc[Indices::turbulentKineticEnergyIdx];
            analyticalDissipation[ccDofIdx] = analyticalSolutionAtCc[Indices::dissipationIdx];

            for(int dirIdx = 0; dirIdx < dim; ++dirIdx)
                analyticalVelocity[ccDofIdx][dirIdx] = analyticalSolutionAtCc[Indices::velocity(dirIdx)];
        }
    }

    return std::make_tuple(analyticalPressure, analyticalTurbulentKineticEnergy, analyticalDissipation, analyticalVelocity, analyticalVelocityOnFace);
}

int main(int argc, char** argv) try
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::TYPETAG;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv, usage);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);
    problem->setTimeLoop(timeLoop);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    x[GridGeometry::cellCenterIdx()].resize(gridGeometry->numCellCenterDofs());
    x[GridGeometry::faceIdx()].resize(gridGeometry->numFaceDofs());
    problem->applyInitialSolution(x);
    problem->updateStaticWallProperties();
    problem->updateDynamicWallProperties(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // intialize the vtk output module
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    StaggeredVtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.write(0.0);

    auto analyticalSolution = createAnalyticalSolution(timeLoop->time(), *problem);
    vtkWriter.addField(std::get<0>(analyticalSolution), "pressureExact");
    vtkWriter.addField(std::get<1>(analyticalSolution), "turbulentKineticEnergyExact");
    vtkWriter.addField(std::get<2>(analyticalSolution), "dissipationExact");
    vtkWriter.addField(std::get<3>(analyticalSolution), "velocityExact");
    vtkWriter.addFaceField(std::get<4>(analyticalSolution), "faceVelocityExact");
    vtkWriter.write(0.0);

    // the assembler with time loop for instationary problem
    using Assembler = StaggeredFVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    // the linear solver
    using LinearSolver = Dumux::UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // time loop
    timeLoop->start(); do
    {
        // solve the non-linear system with time step control
        nonLinearSolver.solve(x, *timeLoop);

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        // update wall properties
        problem->updateDynamicWallProperties(x);
        assembler->updateGridVariables(x);

        if (printL2Error)
        {
            using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
            using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
            using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

            using L2Error = NavierStokesTestL2Error<Scalar, ModelTraits, PrimaryVariables>;
            const auto l2error = L2Error::calculateL2Error(*problem, x);
            const int numCellCenterDofs = gridGeometry->numCellCenterDofs();
            const int numFaceDofs = gridGeometry->numFaceDofs();
            std::cout << std::setprecision(8) << "** L2 error (abs/rel) for "
                    << std::setw(6) << numCellCenterDofs << " cc dofs and " << numFaceDofs << " face dofs (total: " << numCellCenterDofs + numFaceDofs << "): "
                    << std::scientific
                    << "L2(p) = " << l2error.first[Indices::pressureIdx] << " / " << l2error.second[Indices::pressureIdx]
                    << ", L2(k) = " << l2error.first[Indices::turbulentKineticEnergyIdx] << " / " << l2error.second[Indices::turbulentKineticEnergyIdx]
                    << ", L2(w) = " << l2error.first[Indices::dissipationIdx] << " / " << l2error.second[Indices::dissipationIdx]
                    << ", L2(vx) = " << l2error.first[Indices::velocityXIdx] << " / " << l2error.second[Indices::velocityXIdx]
                    << ", L2(vy) = " << l2error.first[Indices::velocityYIdx] << " / " << l2error.second[Indices::velocityYIdx]
                    << std::endl;
        }

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();
        problem->updateTime(timeLoop->time());
        analyticalSolution = createAnalyticalSolution(timeLoop->time(), *problem);

        // write vtk output
        vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton solver
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
        problem->updateTimeStepSize(timeLoop->timeStepSize());

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main
catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
