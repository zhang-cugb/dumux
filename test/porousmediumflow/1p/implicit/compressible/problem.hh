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

namespace Dumux
{
// forward declarations
template<class TypeTag> class OnePTestProblem;
template<class TypeTag> class OnePTestSpatialParams;

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
SET_PROP(IncompressibleTestProblem, Fluid)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    using type = FluidSystems::LiquidPhase<Scalar, SimpleH2O<Scalar> >;
};

// Enable caching
SET_BOOL_PROP(IncompressibleTestProblem, EnableGlobalVolumeVariablesCache, false);
SET_BOOL_PROP(IncompressibleTestProblem, EnableGlobalFluxVariablesCache, false);
SET_BOOL_PROP(IncompressibleTestProblem, EnableFVGridGeometryCache, false);

} // end namespace Properties

template<class TypeTag>
class OnePTestSpatialParams
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimension>;

public:
    using PermeabilityType = Scalar;
    PermeabilityType permeability(const Element &element,
                        const SubControlVolume &scv,
                        const ElementSolutionVector &elemSol) const
    { return 1e-12; }

    PermeabilityType permeabilityAtPos(const GlobalPosition &globalPos) const
    { return 1e-12; }

    Scalar porosity(const Element &element,
                        const SubControlVolume &scv,
                        const ElementSolutionVector &elemSol) const
    { return 0.2; }
};


template<class TypeTag>
class OnePTestProblem
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ResidualVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimension>;

public:
    OnePTestProblem(const GridView& gridView)
    : gridView_(gridView)
    {
        // TabulatedComponent<Scalar, H2O<Scalar>>::init(273.15 - 10,
        //                                               273.15 + 10,
        //                                               3,
        //                                               1e4,
        //                                               1e6,
        //                                               100);
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes values;
        const auto& pos = scvf.ipGlobal();
        values.setAllDirichlet();
        if (pos[0] > 1.0 - 1e-8 || pos[0] < 1e-8)
            values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichlet(const Element &element,
                               const SubControlVolumeFace &scvf) const
    {
        const auto& pos = scvf.ipGlobal();
        return PrimaryVariables(pos[1]*1e5 + 1e5);
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    ResidualVector neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolvars,
                           const SubControlVolumeFace& scvf) const
    {
        return ResidualVector(0.0);
    }

    /*!
     * \brief Applies the initial solution for all degrees of freedom of the grid.
     *
    */
    void applyInitialSolution(SolutionVector& sol) const
    {
        sol = 1e5;
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param values The source and sink values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^3 \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The subcontrolvolume
     *
     * For this method, the \a values parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
    ResidualVector source(const Element &element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolVars,
                          const SubControlVolume &scv) const
    {
        if (scv.dofIndex() == gridView_.size(0)/2)
            return ResidualVector(0.1/scv.volume());
        else
            return ResidualVector(0.0);
    }

    /*!
     * \brief Adds contribution of point sources for a specific sub control volume
     *        to the values.
     *        Caution: Only overload this method in the implementation if you know
     *                 what you are doing.
     */
    ResidualVector scvPointSources(const Element &element,
                                   const FVElementGeometry& fvGeometry,
                                   const ElementVolumeVariables& elemVolVars,
                                   const SubControlVolume &scv) const
    {
        return ResidualVector(0.0);
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at a given global position.
     *
     * This is not specific to the discretization. By default it just
     * calls temperature().
     *
     * \param globalPos The position in global coordinates where the temperature should be specified.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    { return 273.15 + 10; }

    /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * This is discretization independent interface. By default it
     * just calls gravity().
     */
    GlobalPosition gravityAtPos(const GlobalPosition &pos) const
    { return GlobalPosition({0.0, -9.81}); }

    /*!
     * \brief Return how much the domain is extruded at a given sub-control volume.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     */
    Scalar extrusionFactor(const Element &element,
                           const SubControlVolume &scv,
                           const ElementSolutionVector &elemSol) const
    {
        return 1.0;
    }

    const OnePTestSpatialParams<TypeTag>& spatialParams() const
    { return spatialParams_; }

private:
    GridView gridView_;
    OnePTestSpatialParams<TypeTag> spatialParams_;


};

} // end namespace Dumux

#endif
