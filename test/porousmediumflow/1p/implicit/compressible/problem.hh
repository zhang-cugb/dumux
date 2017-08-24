// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief The properties for the incompressible test
 */
#ifndef DUMUX_INCOMPRESSIBLE_ONEP_TEST_PROBLEM_HH
#define DUMUX_INCOMPRESSIBLE_ONEP_TEST_PROBLEM_HH

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/implicit/cellcentered/assembler.hh>
#include <dumux/implicit/cellcentered/localassembler.hh>
#include <dumux/implicit/gridvariables.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/porousmediumflow/1p/implicit/propertydefaults.hh>

#include "spatialparams.hh"

namespace Dumux
{
// forward declarations
template<class TypeTag> class OnePTestProblem;

namespace Properties
{

NEW_PROP_TAG(EnableFVGridGeometryCache);
NEW_PROP_TAG(FVGridGeometry);

NEW_TYPE_TAG(IncompressibleTestProblem, INHERITS_FROM(CCTpfaModel, OneP));

// Set the grid type
SET_TYPE_PROP(IncompressibleTestProblem, Grid, Dune::YaspGrid<2>);

// Set the problem type
SET_TYPE_PROP(IncompressibleTestProblem, Problem, OnePTestProblem<TypeTag>);
SET_TYPE_PROP(IncompressibleTestProblem, SpatialParams, OnePTestSpatialParams<TypeTag>);

// the grid variables
SET_TYPE_PROP(IncompressibleTestProblem, GridVariables, GridVariables<TypeTag>);

// the grid variables
SET_TYPE_PROP(IncompressibleTestProblem, JacobianAssembler, CCImplicitAssembler<TypeTag>);

// linear solver
SET_TYPE_PROP(IncompressibleTestProblem, LinearSolver, ILU0BiCGSTABBackend<TypeTag>);

// the local assembler
SET_TYPE_PROP(IncompressibleTestProblem, LocalAssembler, CCImplicitLocalAssembler<TypeTag, DifferentiationMethods::numeric>);

// the fluid system
SET_PROP(IncompressibleTestProblem, FluidSystem)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    using type = FluidSystems::LiquidPhase<Scalar, TabulatedComponent<Scalar, H2O<Scalar>>>;
};

// Enable caching
SET_BOOL_PROP(IncompressibleTestProblem, EnableGlobalVolumeVariablesCache, false);
SET_BOOL_PROP(IncompressibleTestProblem, EnableGlobalFluxVariablesCache, false);
SET_BOOL_PROP(IncompressibleTestProblem, EnableFVGridGeometryCache, false);

} // end namespace Properties

template<class TypeTag>
class OnePTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimension>;

    static constexpr int dimWorld = GridView::dimensionworld;

public:
    OnePTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        TabulatedComponent<Scalar, H2O<Scalar>>::init(272.15, 294.15, 10,
                                                      1.0e4, 1.0e6, 200);
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes values;
        const auto globalPos = scvf.ipGlobal();

        Scalar eps = 1.0e-6;
        if (globalPos[dimWorld-1] < eps || globalPos[dimWorld-1] > this->fvGridGeometry().bBoxMax()[dimWorld-1] - eps)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    PrimaryVariables dirichlet(const Element &element,
                               const SubControlVolumeFace &scvf) const
    {
        const auto& pos = scvf.ipGlobal();
        PrimaryVariables values(0);
        values[0] = 1.0e+5*(2.0 - pos[dimWorld-1]);
        return values;
    }

    /*!
     * \brief Evaluate the initial conditions
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        return PrimaryVariables(1.0e5);
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
};

} // end namespace Dumux

#endif
