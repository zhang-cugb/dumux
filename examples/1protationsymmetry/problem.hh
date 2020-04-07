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

#ifndef DUMUX_ONEP_ROTATION_SYMMETRY_PROBLEM_HH
#define DUMUX_ONEP_ROTATION_SYMMETRY_PROBLEM_HH

// ## The problem class (`problem.hh`)
// This file contains the __problem class__ which defines the initial and boundary
// conditions for the single-phase flow simulation.
// [[content]]
// ### Includes
#include <cmath> // for `std::log`
#include <dumux/common/properties.hh> // for `GetPropType`
#include <dumux/common/parameters.hh> // for `getParam`
#include <dumux/porousmediumflow/problem.hh>  // for `PorousMediumFlowProblem`

// ### The problem class
// We enter the problem class where all necessary boundary conditions and initial conditions are set for our simulation.
// As this is a porous medium flow problem, we inherit from the base class `PorousMediumFlowProblem`.
namespace Dumux {

template<class TypeTag>
class RotSymExampleProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    // This is the constructor of our problem class:
    RotSymExampleProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        k_ = getParam<Scalar>("SpatialParams.Permeability");
        nu_ = getParam<Scalar>("Component.LiquidKinematicViscosity");
        q_ = getParam<Scalar>("Problem.Source");
        pW_ = getParam<Scalar>("Problem.WellPressure");
        rW_ = gridGeometry->bBoxMin()[0];
    }

    // First, we define the type of boundary conditions depending on the location. Two types of boundary  conditions
    // can be specified: Dirichlet or Neumann boundary condition. On a Dirichlet boundary, the values of the
    // primary variables need to be fixed. On a Neumann boundary condition, values for derivatives need to be fixed.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    // Second, we specify the values for the Dirichlet boundaries. We need to fix values of our primary variable
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    { return exactSolution(globalPos); }

    // We need to specify a constant temperature for our isothermal problem.
    // Fluid properties that depend on temperature will be calculated with this value.
    Scalar temperature() const
    { return 283.15; }

    // This function defines the exact pressure solution
    PrimaryVariables exactSolution(const GlobalPosition &globalPos) const
    {
        const auto r = globalPos[0];
        const auto p = pW_ - 1.0/(2*M_PI)*nu_/k_*q_*std::log(r/rW_);
        return p;
    }

private:
    Scalar q_, k_, nu_, rW_, pW_;
};

} // end namespace Dumux
// [[/content]]
#endif
