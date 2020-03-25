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
 * \ingroup OnePTests
 * \brief The properties and problem setup for rotation symmetry test
 */
#ifndef DUMUX_ONEP_ROTATION_SYMMETRY_PROBLEM_HH
#define DUMUX_ONEP_ROTATION_SYMMETRY_PROBLEM_HH

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/common/integrate.hh>

#include "spatialparams.hh"

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief  Test problem for the incompressible one-phase model
 */
template<class TypeTag>
class OnePTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

public:
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

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     *
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return exactSolution(globalPos);
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
     *
     * This is not specific to the discretization. By default it just
     * throws an exception so it must be overloaded by the problem if
     * no energy equation is used.
     */
    Scalar temperature() const
    {
        return 283.15; // 10°C
    }

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    template<class VTKWriter>
    void addVtkFields(VTKWriter& vtk) const
    {
        vtk.addField(pExact_, "pExact");
    }

    /*!
     * \brief Writes the L2 error
     */
    template<class SolutionVector>
    Scalar calculateL2Error(const SolutionVector& curSol, const int order)
    {
        return integrateL2Error(this->gridGeometry(), curSol, pExact_, order);
    }

private:
    /*!
     * \brief The exact solution
     * The exact solution is calculated such that the mass flux over the surface of a circular disc with radius rW is q
     * and the pressure at this surface is given by pW.
     */
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
};

} // end namespace Dumux

#endif
