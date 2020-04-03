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

// ### Header guard
#ifndef DUMUX_ONEP_ROTATION_SYMMETRY_PROBLEM_HH
#define DUMUX_ONEP_ROTATION_SYMMETRY_PROBLEM_HH

// This file contains the __problem class__ which defines the initial and boundary
// conditions for the single-phase flow simulation.
//
// ### Include files
// This header contains the porous medium problem class that this class is derived from:
#include <dumux/porousmediumflow/problem.hh>
// This header contains a convience function to calculate L2 errors
#include <dumux/common/integrate.hh>
// This header contains the class that specifies all spatially variable parameters
// related to this problem.
#include "spatialparams.hh"

// ### The problem class
// We enter the problem class where all necessary boundary conditions and initial conditions are set for our simulation.
// As this is a porous medium flow problem, we inherit from the base class `PorousMediumFlowProblem`.
namespace Dumux {

template<class TypeTag>
class OnePTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    // We use convenient declarations that we derive from the property system.
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

public:
    // This is the constructor of our problem class:
    OnePTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry), q_(0.0), pExact_(gridGeometry->numDofs())
    {
        k_ = getParam<Scalar>("SpatialParams.Permeability");
        nu_ = getParam<Scalar>("Component.LiquidKinematicViscosity");
        q_ = getParam<Scalar>("Problem.Source");
        pW_ = getParam<Scalar>("Problem.WellPressure");
        rW_ = gridGeometry->bBoxMin()[0];

        for (const auto& element : elements(gridGeometry->gridView()))
        {
            auto fvGeometry = localView(*gridGeometry);
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                pExact_[scv.dofIndex()] = exactSolution(scv.dofPosition());
            }
        }
    }

    // First, we define the type of boundary conditions depending on the location. Two types of boundary  conditions
    // can be specified: Dirichlet or Neumann boundary condition. On a Dirichlet boundary, the values of the
    // primary variables need to be fixed. On a Neumann boundary condition, values for derivatives need to be fixed.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        // We specify Dirichlet boundaries everywhere
        values.setAllDirichlet();

        return values;
    }

    // Second, we specify the values for the Dirichlet boundaries. We need to fix values of our primary variable
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        // The exact solution values are set as Dirichlet values
        return exactSolution(globalPos);
    }

    // We need to specify a constant temperature for our isothermal problem.
    // Fluid properties that depend on temperature will be calculated with this value.
    Scalar temperature() const
    {
        return 283.15; // 10°C
    }

    // This method add the exact pressure values to the vtk output
    template<class VTKWriter>
    void addVtkFields(VTKWriter& vtk) const
    {
        vtk.addField(pExact_, "pExact");
    }

    // The L2 error between the exact and numerical solution is calculated using this function,
    // using a specific order for the quadrature rule.
    template<class SolutionVector>
    Scalar calculateL2Error(const SolutionVector& curSol, const int order)
    {
        return integrateL2Error(this->gridGeometry(), curSol, pExact_, order);
    }

private:
    // This function defines the exact pressure solution
    PrimaryVariables exactSolution(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        const auto r = globalPos[0];
        priVars[0] = pW_ - 1.0/(2*M_PI)*nu_/k_*q_*std::log(r/rW_);
        return priVars;
    }

    Scalar q_, k_, nu_, rW_;
    GlobalPosition pW_;
    static constexpr Scalar eps_ = 1.5e-7;
    SolutionVector pExact_;

    // This is everything the one phase rotation symmetry problem class contains.
};

// We leave the namespace Dumux.
} // end namespace Dumux

#endif
