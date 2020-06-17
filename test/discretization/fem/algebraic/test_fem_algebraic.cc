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
 * \brief Test solving a time-dependent algebraic problem with the
 *        finite element method.
 */
#include <config.h>
#include <iostream>

#include <dune/istl/bvector.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>

#include <dumux/io/grid/gridmanager.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/integrate.hh>

#include <dumux/assembly/instationaryassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/timestepping/multistagemethods.hh>
#include <dumux/timestepping/multistagetimestepper.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/discretization/fem/fegridgeometry.hh>
#include <dumux/discretization/fem/fegridvariables.hh>

#include "problem.hh"
#include "operators.hh"
#include <dumux/assembly/fem/localoperator.hh> // TODO: order is like this because of the problem traits issue

template<class Scalar, class GridVariables, class TimeStepper>
std::vector<Scalar> computeErrors(GridVariables& gridVariables,
                                  TimeStepper& timeStepper)
{
    auto x = gridVariables.dofs(); x = 0.0;
    gridVariables.update(x, /*timeLevel*/0.0);
unsigned int count = 0;
    std::vector<Scalar> errors;
    Dumux::TimeLoop<Scalar> timeLoop(/*tInit*/0.0,
                                     Dumux::getParam<Scalar>("TimeLoop.Dt"),
                                     Dumux::getParam<Scalar>("TimeLoop.TEnd"));
    timeLoop.start(); do
    {
        // do time integration
        timeStepper.step(gridVariables, timeLoop.time(), timeLoop.timeStepSize());

        // advance to the time loop to the next step
        timeLoop.advanceTimeStep();

        // compute error
        const auto time = timeLoop.time();
        const auto& problem = gridVariables.problem();
        const auto& gridGeometry = gridVariables.gridGeometry();
        const auto& gridView = gridGeometry.gridView();

        using namespace Dune::Functions;
        using PrimaryVariables = typename GridVariables::PrimaryVariables;
        auto evalExact = [&] (const auto& pos) { return problem.exactSolution(time); };
        auto exactGF = makeAnalyticGridViewFunction(evalExact, gridView);
        auto numericGF = makeDiscreteGlobalBasisFunction<PrimaryVariables>(gridGeometry.feBasis(), gridVariables.dofs());

        Dune::VTKWriter<typename GridVariables::GridGeometry::GridView> vtkWriter(gridView);
        vtkWriter.addVertexData(exactGF, Dune::VTK::FieldInfo("exact", Dune::VTK::FieldInfo::Type::scalar, 1));
        vtkWriter.addVertexData(numericGF, Dune::VTK::FieldInfo("numeric", Dune::VTK::FieldInfo::Type::scalar, 1));
        vtkWriter.write("output_" + std::to_string(count++));

        const auto error = Dumux::integrateL2Error(gridView, exactGF, numericGF, 2);

        // normalize error by grid volume
        Scalar volume = 0.0;
        for (const auto& e : elements(gridView))
            volume += e.geometry().volume();

        errors.push_back( error/volume );

    } while (!timeLoop.finished());

    return errors;
}

int main (int argc, char *argv[]) try
{
    using namespace Dumux;

    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // initialize parameter tree
    Parameters::init(argc, argv);

    using Grid = Dune::YaspGrid<2>;
    using GridManager = GridManager<Grid>;
    GridManager gridManager;
    gridManager.init();

    using GridView = typename Grid::LeafGridView;
    const GridView& gridView = gridManager.grid().leafGridView();

    using Scalar = double;
    static constexpr int basisOrder = 1;

    // make the grid geometry
    using FEBasis = Dune::Functions::LagrangeBasis<GridView, basisOrder, Scalar>;
    using GridGeometry = FEGridGeometry<FEBasis>;
    auto feBasis = std::make_shared<FEBasis>(gridView);
    auto gridGeometry = std::make_shared<GridGeometry>(feBasis);
    gridGeometry->update();

    // make the problem
    using Problem = FEAlgebraicProblem<Scalar, GridGeometry>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // make the grid variables
    using PrimaryVariables = typename ProblemTraits<Problem>::PrimaryVariables;
    using SolutionVector = Dune::BlockVector<PrimaryVariables>;
    using GridVariables = FEGridVariables<Problem, SolutionVector>;
    auto gridVariables = std::make_shared<GridVariables>(problem, [&] (auto& dofs) {
        dofs.resize(gridGeometry->numDofs());
    });

    // the local operator type
    using Operators = FEAlgebraicOperators<typename GridVariables::LocalView>;
    using LocalOperator = FELocalOperator<typename GridVariables::LocalView, Operators>;

    // make the assembler
    using Assembler = InstationaryAssembler<LocalOperator, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(gridGeometry);

    // non-linear solver (jacobian is inexact due to numeric differentiation)
    using LinearSolver = ILU0BiCGSTABBackend;
    using NewtonSolver = NewtonSolver<Assembler, LinearSolver>;
    auto linearSolver = std::make_shared<LinearSolver>();
    auto newtonSolver = std::make_shared<NewtonSolver>(assembler, linearSolver);

    // lambda for printing error values
    auto printErrors = [] (const auto& schemeErrors, const std::string& name)
    {
        std::cout << "Errors for " + name << ": ";
        for (unsigned int i = 0; i < schemeErrors.size(); ++i)
            std::cout << std::left << std::setw(9) << std::setfill(' ')
                      << schemeErrors[i] << (i < schemeErrors.size()-1 ? ", " : "\n");
    };

    using TimeStepper = MultiStageTimeStepper<NewtonSolver>;
    // 1. explicit Euler
    {
        auto expEuler = std::make_shared< MultiStage::ExplicitEuler<Scalar> >();
        TimeStepper timeStepper(newtonSolver, expEuler);
        printErrors(computeErrors<Scalar>(*gridVariables, timeStepper), "explicit Euler    ");
    }

    // 2. Theta
    {
        auto theta = std::make_shared< MultiStage::Theta<Scalar> >(0.5);
        TimeStepper timeStepper(newtonSolver, theta);
        printErrors(computeErrors<Scalar>(*gridVariables, timeStepper), "theta scheme (0.5)");
    }

    // 3. implicit Euler
    {
        auto impEuler = std::make_shared< MultiStage::ImplicitEuler<Scalar> >();
        TimeStepper timeStepper(newtonSolver, impEuler);
        printErrors(computeErrors<Scalar>(*gridVariables, timeStepper), "implicit Euler    ");
    }

    // 4. Runge-Kutta 4th order (explicit)
    {
        auto rk4 = std::make_shared< MultiStage::RungeKuttaExplicitFourthOrder<Scalar> >();
        TimeStepper timeStepper(newtonSolver, rk4);
        printErrors(computeErrors<Scalar>(*gridVariables, timeStepper), "Runge-Kutta 4th   ");
    }

    return 0;
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Dune::Exception &e) {
    std::cout << e << std::endl;
    return 1;
}
