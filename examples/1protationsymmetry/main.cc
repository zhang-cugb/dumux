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
// ## The main program (`main.cc`)
// We look now at the main file for the tracer problem. We set up two problems in this file and solve them sequentially, first the 1p problem and afterwards the tracer problem. The result of the 1p problem is the pressure distribution in the problem domain. We use it to calculate the volume fluxes, which act as an input for the tracer problem. Based on this volume fluxes, we calculate the transport of a tracer in the following tracer problem.
// [[content]]
// ### Includes
// [[details]] includes
#include <config.h>

#include <iostream>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/integrate.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/linear/pdesolver.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include "properties.hh"
// [[/details]]

// ### Beginning of the main function
int main(int argc, char** argv) try
{
    using namespace Dumux;

    // Convenience aliases for the type tag of the problem.
    using TypeTag = Properties::TTag::OnePRotSym;

    // We initialize MPI. Finalization is done automatically on exit.
    Dune::MPIHelper::instance(argc, argv);

    // We parse the command line arguments.
    Parameters::init(argc, argv);

    // ### Create the grid and the grid geometry
    // [[codeblock]]
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // We compute on the leaf grid view.
    const auto& leafGridView = gridManager.grid().leafGridView();
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();
    // [[/codeblock]]

    // ### Initialise the problem
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // a function to update the discrete analytical solution vector
    const auto updateAnalyticalSolution = [&](auto& pExact)
    {
        pExact.resize(gridGeometry->numDofs());
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            auto fvGeometry = localView(*gridGeometry);
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
                pExact[scv.dofIndex()] = problem->exactSolution(scv.dofPosition());
        }
    };

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector p(gridGeometry->numDofs());
    SolutionVector pExact; updateAnalyticalSolution(pExact);

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(p);

    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, p, problem->name());
    GetPropType<TypeTag, Properties::IOFields>::initOutputModule(vtkWriter);
    vtkWriter.addField(pExact, "pExact");

    using Assembler = FVAssembler<TypeTag, DiffMethod::analytic>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LinearSolver = UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();
    LinearPDESolver<Assembler, LinearSolver> solver(assembler,  linearSolver);
    solver.setVerbose(false);

    // solve once with initial refinement
    solver.solve(p);

    // compute initial L2 error
    constexpr bool isBox = GridGeometry::discMethod == Dumux::DiscretizationMethod::box;
    constexpr int orderQuadratureRule = isBox ? 3 : 1;
    const int numRefinements = getParam<int>("Grid.RefinementSteps");
    std::vector<double> l2Errors(numRefinements);
    l2Errors[0] = integrateL2Error(*gridGeometry, p, pExact, orderQuadratureRule);

    // repeat for several refinements
    // [[codeblock]]
    for (int stepIdx = 1; stepIdx < numRefinements; stepIdx++)
    {
        // Globally refine the grid once
        gridManager.grid().globalRefine(1);
        gridGeometry->update();
        p.resize(gridGeometry->numDofs());
        updateAnalyticalSolution(pExact);
        gridVariables->updateAfterGridAdaption(p);
        assembler->setLinearSystem();

        // solve problem on refined grid
        solver.solve(p);

        // #### Post-processing and output
        // We calculate the L2 errors using the numerical solution
        l2Errors[stepIdx] = integrateL2Error(*gridGeometry, p, pExact, orderQuadratureRule);
        const auto numDofs = gridGeometry->numDofs();
        std::cout << std::setprecision(8) << std::scientific
                  << "-- L2 error for " << std::setw(5) << numDofs << " dofs: " << l2Errors[stepIdx]
                  << ", rate: " << std::log(l2Errors[stepIdx]/l2Errors[stepIdx-1])/std::log(0.5)
                  << std::endl;
    }
    // write vtk output on the finest grid
    vtkWriter.write(0.0);
    // [[/codeblock]]

    // program end, return with 0 exit code (success)
    return 0;
}
// ### Exception handling
// In this part of the main file we catch and print possible exceptions that could
// occur during the simulation.
// [[details]] error handler
catch (const Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (const Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (const Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
// [[/details]]
// [[/content]]
