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

// We look now at the main file for the tracer problem. We set up two problems in this file and solve them sequentially, first the 1p problem and afterwards the tracer problem. The result of the 1p problem is the pressure distribution in the problem domain. We use it to calculate the volume fluxes, which act as an input for the tracer problem. Based on this volume fluxes, we calculate the transport of a tracer in the following tracer problem.
// ### Includes
#include <config.h>

// This includes the `TypeTags` and properties to be used for the single-phase rotation symmetry example.
#include "properties.hh"

// Further, we include a standard header file for C++, to get time and date information
#include <ctime>
// and another one for in- and output.
#include <iostream>

// Dumux is based on DUNE, the Distributed and Unified Numerics Environment, which provides several grid managers and linear solvers.
// Here, we include classes related to parallel computations, time measurements and file I/O.
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>

// In Dumux, the property system is used to specify classes and compile-time options to be used by the model.
// For this, different properties are defined containing type definitions, values and methods.
// All properties are declared in the file `properties.hh`.
#include <dumux/common/properties.hh>
// The following file contains the parameter class, which manages the definition and retrieval of input
// parameters by a default value, the inputfile or the command line.
#include <dumux/common/parameters.hh>
// The file `dumuxmessage.hh` contains the class defining the start and end message of the simulation.
#include <dumux/common/dumuxmessage.hh>

// The following file contains the class, which defines the sequential linear solver backends.
#include <dumux/linear/seqsolverbackend.hh>
// The following linear pde solver allows to solve linear equations
#include <dumux/linear/pdesolver.hh>
// Further we include the assembler, which assembles the linear systems for finite volume schemes (box-scheme, tpfa-approximation, mpfa-approximation).
#include <dumux/assembly/fvassembler.hh>
// The containing class in the following file defines the different differentiation methods used to compute the derivatives of the residual.
#include <dumux/assembly/diffmethod.hh>

// We need the following class to simplify the writing of dumux simulation data to VTK format.
#include <dumux/io/vtkoutputmodule.hh>
// The gridmanager constructs a grid from the information in the input or grid file. There is a specification for the different supported grid managers.
#include <dumux/io/grid/gridmanager.hh>

// ### Beginning of the main function
int main(int argc, char** argv) try
{
    using namespace Dumux;

    // Convenience aliases for the type tag of the problem.
    using TypeTag = Properties::TTag::OnePRotSymBox;

    // We initialize MPI. Finalization is done automatically on exit.
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // We print the dumux start message.
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // We parse the command line arguments.
    Parameters::init(argc, argv);

    // ### Create the grid
    // The `GridManager` class creates the grid from information given in the input file.
    // This can either be a grid file, or in the case of structured grids, by specifying the coordinates
    // of the corners of the grid and the number of cells to be used to discretize each spatial direction.
    // Here, we solve both the single-phase and the tracer problem on the same grid.
    // Hence, the grid is only created once using the grid type defined by the type tag of the 1p problem.
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // start timer
    Dune::Timer timer;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;

    int numRefinements = getParam<int>("Grid.RefinementSteps");
    constexpr bool isBox = GridGeometry::discMethod == Dumux::DiscretizationMethod::box;
    const int orderQuadratureRule = isBox ? 3 : 1;
    std::vector<Scalar> l2Errors(numRefinements, 0.0);
    for(int i = 0; i < numRefinements; i++)
    {
        // We compute on the leaf grid view.
        const auto& leafGridView = gridManager.grid().leafGridView();

        // ### Set-up and solving of the 1p rotation symmetry problem
        // In the following section, we set up and solve the problem. As the result of this problem, we obtain the pressure distribution in the domain.
        // #### Set-up
        // We create and initialize the finite volume grid geometry, the problem, the linear system, including the jacobian matrix, the residual and the solution vector and the gridvariables.
        // We need the finite volume geometry to build up the subcontrolvolumes (scv) and subcontrolvolume faces (scvf) for each element of the grid partition.
        auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
        gridGeometry->update();

        // In the problem, we define the boundary and initial conditions.
        using Problem = GetPropType<TypeTag, Properties::Problem>;
        auto problem = std::make_shared<Problem>(gridGeometry);

        // The solution vector containing the pressure values
        using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
        SolutionVector x(gridGeometry->numDofs());

        // The grid variables store variables (primary and secondary variables) on sub-control volumes and faces (volume and flux variables).
        using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        gridVariables->init(x);

        // #### Assembling the linear system
        // We create the assembler.
        using Assembler = FVAssembler<TypeTag, DiffMethod::analytic>;
        auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

        // #### Solution
        // We set the linear solver "UMFPack" as the linear solver. Afterwards we solve the linear system.
        // The linear system of equations is solved
        using LinearSolver = UMFPackBackend;
        auto linearSolver = std::make_shared<LinearSolver>();
        LinearPDESolver<Assembler, LinearSolver> solver(assembler,  linearSolver);
        solver.setVerbose(0);
        solver.solve(x);

        // #### Post-processing and output
        // We calculate the L2 errors using the numerical solution
        l2Errors[i] = problem->calculateL2Error(x, orderQuadratureRule);
        std::cout.precision(8);
        std::cout << "L2 error for "
                  << std::setw(6) << gridGeometry->numDofs()
                  << " dofs: "
                  << std::scientific
                  << l2Errors[i]
                  << " rate: "
                  << ((i>0) ? std::to_string(std::log(l2Errors[i]/l2Errors[i-1])/(std::log(0.5))) : " - ")
                  << std::endl;


        // The vtk file is written on the finest grid.
        if(i == numRefinements-1)
        {
            // output result to vtk
            VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
            using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
            vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
            using IOFields = GetPropType<TypeTag, Properties::IOFields>;
            IOFields::initOutputModule(vtkWriter); // Add model specific output fields
            problem->addVtkFields(vtkWriter); //!< Add problem specific output fields
            vtkWriter.write(0.0);
        }

        // Globally refine the grid and repeate the solution procedure
        gridManager.grid().globalRefine(1);
    }

    // We stop the timer and display the total time of the simulation as well as the cumulative CPU time.
    timer.stop();

    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    if (mpiHelper.rank() == 0)
        std::cout << "Simulation took " << timer.elapsed() << " seconds on "
                  << comm.size() << " processes.\n"
                  << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";

    // ### Final Output, print parameters
    if (mpiHelper.rank() == 0)
        Parameters::print();

    return 0;
}
// ### Exception handling
// In this part of the main file we catch and print possible exceptions that could
// occur during the simulation.
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
