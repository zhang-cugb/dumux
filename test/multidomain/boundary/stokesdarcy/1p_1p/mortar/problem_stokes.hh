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
 * \ingroup BoundaryTests
 * \brief A simple Stokes test problem for the staggered grid (Navier-)Stokes model.
 */

#ifndef DUMUX_STOKES_SUBPROBLEM_HH
#define DUMUX_STOKES_SUBPROBLEM_HH

#include <utility>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/simpleh2o.hh>

#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/model.hh>

namespace Dumux {
template <class TypeTag>
class StokesSubProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct StokesOneP { using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::StokesOneP>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Dumux::Components::Constant<0, Scalar> > ;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::StokesOneP> { using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::StokesOneP> { using type = Dumux::StokesSubProblem<TypeTag> ; };

template<class TypeTag>
struct EnableFVGridGeometryCache<TypeTag, TTag::StokesOneP> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::StokesOneP> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::StokesOneP> { static constexpr bool value = true; };
} // end namespace Properties

/*!
 * \ingroup BoundaryTests
 * \brief Test problem for the one-phase (Navier-) Stokes problem.
 *
 * Horizontal flow from left to right with a parabolic velocity profile.
 */
template <class TypeTag>
class StokesSubProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

    // extract sub-vector for cell-centered values (used for mortar data)
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using CellSolutionVector = std::decay_t<decltype( std::declval<SolutionVector>()[FVGridGeometry::cellCenterIdx()] )>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    StokesSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry, const std::string& paramGroup)
    : ParentType(fvGridGeometry, paramGroup)
    , eps_(1e-6)
    , isOnNegativeMortarSide_(getParamFromGroup<bool>(paramGroup, "Problem.IsOnNegativeMortarSide"))
    {
        const auto mortarVariable = getParamFromGroup<std::string>("Mortar", "VariableType");
        if (mortarVariable == "Pressure")
            useDirichletAtInterface_ = false;
        else if (mortarVariable == "Flux")
            DUNE_THROW(Dune::NotImplemented, "Flux coupling stokes problem");

        problemName_  =  getParamFromGroup<std::string>(paramGroup, "Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    { return problemName_; }

   /*!
     * \name Problem parameters
     */
    // \{

   /*!
     * \brief Returns the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10°C

    // \}

   /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scvf The sub control volume face
     */
    BoundaryTypes boundaryTypes(const Element& element,
                                const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;

        const auto& globalPos = scvf.dofPosition();

        if(onLeftBoundary_(globalPos) || onRightBoundary_(globalPos))
        {
            values.setDirichlet(Indices::pressureIdx);
            // values.setDirichlet(Indices::velocityXIdx);
            // values.setDirichlet(Indices::velocityYIdx);
        }
        // values.setDirichlet(Indices::pressureIdx);

        else if (isOnMortarInterface(globalPos))
        {
            if (useDirichletAtInterface_)
                values.setDirichlet(scvf.directionIndex());
            else
                values.setNeumann(scvf.directionIndex());
            values.setBJS(1 - scvf.directionIndex());
        }
        else
        {
            // values.setDirichlet(Indices::pressureIdx);
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }

        return values;
    }

    /*!
     * \brief Evaluates the velocity boundary conditions for a Dirichlet sub-control volume face.
     */
    PrimaryVariables dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        if (!useHomogeneousSetup_)
        {
            static const Scalar deltaP = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.PressureDifference");

            PrimaryVariables values(0.0);
            if(onLeftBoundary_(scvf.ipGlobal()))
                values[Indices::pressureIdx] = deltaP;
            return values;
        }
        else
            return initialAtPos(scvf.ipGlobal());
    }

    /*!
     * \brief Evaluates the pressure boundary conditions for a sub-control volume
     *        touchgin a dirichlet boundary segment.
     */
    PrimaryVariables dirichlet(const Element& element, const SubControlVolume& scv) const
    {
        if (!useHomogeneousSetup_)
        {
            static const Scalar deltaP = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.PressureDifference");

            PrimaryVariables values(0.0);
            if(onLeftBoundary_(scv.dofPosition()))
                values[Indices::pressureIdx] = deltaP;
            return values;
        }
        else
            return initialAtPos(scv.dofPosition());
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables, class ElementFaceVariables>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFaceVariables& elemFaceVars,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        if (isOnMortarInterface(scvf.ipGlobal()) && !useDirichletAtInterface_)
        {
            // apply mortar pressure to momentum balance
            assert(mortarProjection_.size() == this->fvGridGeometry().gridView().size(0));
            const auto eIdx = this->fvGridGeometry().elementMapper().index(element);
            values[scvf.directionIndex()] = mortarProjection_[eIdx]*scvf.directionSign();

            // mass flux
            const auto v = elemFaceVars[scvf].velocitySelf()*scvf.directionSign();
            values[Indices::conti0EqIdx] = v*elemVolVars[scvf.insideScvIdx()].density();
        }

        return values;
    }

   /*!
     * \name Volume terms
     */
    // \{

   /*!
     * \brief Evaluates the initial value for a control volume.
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    //! returns the intrinsic permeability
    Scalar permeability(const Element& element, const SubControlVolumeFace& scvf) const
    {
        static const Scalar k = getParam<Scalar>("SpatialParams.Permeability");
        return k;
    }

    //! returns the alpha value
    Scalar alphaBJ(const SubControlVolumeFace& scvf) const
    {
        static const Scalar alpha = getParam<Scalar>("SpatialParams.AlphaBeaversJoseph");
        return alpha;
    }

    //! set the pointer to the projector class
    void setMortarProjection(CellSolutionVector p)
    {
        mortarProjection_ = p;
    }

    //! Set whether or not the homogeneous system is solved
    void setUseHomogeneousSetup(bool value)
    {
        useHomogeneousSetup_ = value;
    }

    //! Returns true if a position if on the mortar interface
    bool isOnMortarInterface(const GlobalPosition& globalPos) const
    {
        return (isOnNegativeMortarSide_ && onLowerBoundary_(globalPos))
               || (!isOnNegativeMortarSide_ && onUpperBoundary_(globalPos));
    }

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->fvGridGeometry().bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - eps_; }


    Scalar eps_;
    std::string problemName_;
    CellSolutionVector mortarProjection_;

    bool isOnNegativeMortarSide_;
    bool useHomogeneousSetup_;
    bool useDirichletAtInterface_;
};
} // end namespace Dumux

#endif // DUMUX_STOKES_SUBPROBLEM_HH
