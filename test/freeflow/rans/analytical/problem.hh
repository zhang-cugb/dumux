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
 * \brief Pipe flow test for the staggered grid RANS model
 *
 * This test simulates pipe flow experiments performed by John Laufer in 1954
 * \cite Laufer1954a.
 */
#ifndef DUMUX_PIPE_LAUFER_PROBLEM_HH
#define DUMUX_PIPE_LAUFER_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>
#include <dune/common/hybridutilities.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>

#include <dumux/freeflow/turbulencemodel.hh>
#include <dumux/freeflow/turbulenceproperties.hh>
#include <dumux/freeflow/rans/problem.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/freeflow/rans/twoeq/komega/problem.hh>
#include <dumux/freeflow/rans/twoeq/komega/model.hh>

#include "komegaanalytical.hh"

namespace Dumux {

template <class TypeTag>
class AnalyticalRANSProblem;

namespace Properties {

// Create new type tags
namespace TTag {

struct RANSModel { using InheritsFrom = std::tuple<StaggeredFreeFlowModel>; };
    struct AnalyticalKOmega { using InheritsFrom = std::tuple<RANSModel, KOmega>; };
//  struct AnalyticalKEpsilon { using InheritsFrom = std::tuple<RANSModel, KEpsilon>; };

} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RANSModel>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::RANSModel>
{ using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::RANSModel>
{ using type = Dumux::AnalyticalRANSProblem<TypeTag>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::RANSModel> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::RANSModel> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::RANSModel> { static constexpr bool value = true; };
} // end namespace Properties

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the one-phase (Navier-) Stokes problem in a channel.
 *
 * This test simulates is based on pipe flow experiments by
 * John Laufers experiments in 1954 \cite Laufer1954a.
 */
template <class TypeTag>
class AnalyticalRANSProblem : public RANSProblem<TypeTag>
{
    using ParentType = RANSProblem<TypeTag>;

    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

    struct TurbulentConstants
    {
        Scalar alpha_, betaK_, sigmaK_, betaOmega_, sigmaOmega_;
        TurbulentConstants(const Scalar alpha = 0.0, const Scalar betaK = 0.0, const Scalar sigmaK = 0.0, const Scalar betaOmega = 0.0, const Scalar sigmaOmega = 0.0)
        : alpha_(alpha)
        , betaK_(betaK)
        , sigmaK_(sigmaK)
        , betaOmega_(betaOmega)
        , sigmaOmega_(sigmaOmega)
        {}
    };

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using AnalyticalModel = typename Dumux::KOmegaAnalytical<Scalar, GlobalPosition, PrimaryVariables, ModelTraits, TurbulentConstants, dimWorld>;

public:
    AnalyticalRANSProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry), eps_(1e-6)
    {
        TurbulentConstants turbulentConstants_(ParentType::alpha(), ParentType::betaK(), ParentType::sigmaK(), ParentType::betaOmega(), ParentType::sigmaOmega());

        turbulenceModelName_ = turbulenceModelToString(ModelTraits::turbulenceModel());
        std::cout << "Using the "<< turbulenceModelName_ << " Turbulence Model. \n";
        std::cout << std::endl;
    }

   /*!
     * \name Problem parameters
     */
    // \{

    NumEqVector analyticalSolutionAtPos(const GlobalPosition globalPos) const
    {
        const AnalyticalModel analyticalModel_(density_, constantKinematicViscosity_, turbulentConstants_, this->gravity());
        return analyticalModel_.solution(globalPos);
    }

    bool isOnWallAtPos(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_
               || globalPos[0] > this->gridGeometry().bBoxMax()[0] + eps_
               || globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_
               || globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_;
    }

   /*!
     * \brief Returns the temperature [K] within the domain for the isothermal model.
     */
    Scalar temperature() const
    { return temperature_; }

   /*!
     * \brief Returns the sources within the domain.
     *
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        return NumEqVector(analyticalSolutionAtPos(globalPos));
    }
    // \}

   /*!
     * \name Boundary conditions
     */
    // \{

   /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setDirichlet(Indices::pressureIdx);
        values.setNeumann(Indices::velocityXIdx);
        values.setNeumann(Indices::velocityYIdx);
        values.setDirichlet(Indices::turbulentKineticEnergyIdx);
        values.setDirichlet(Indices::dissipationIdx);

        return values;
    }

     /*!
      * \brief Evaluate the boundary conditions for a dirichlet values at the boundary.
      *
      * \param element The finite element
      * \param scvf the sub control volume face
      * \note used for cell-centered discretization schemes
      */
    NumEqVector dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        const auto globalPos = scvf.ipGlobal();
        NumEqVector values(initialAtPos(globalPos));

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for fixed values at cell centers
     *
     * \param element The finite element
     * \param scv the sub control volume
     * \note used for cell-centered discretization schemes
     */
    NumEqVector dirichlet(const Element& element, const SubControlVolume& scv) const
    {
        const auto globalPos = scv.center();
        NumEqVector values(initialAtPos(globalPos));
        return values;
    }

   /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    NumEqVector initialAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector values(analyticalSolutionAtPos(globalPos));
        return values;
    }

    // \}

private:
    Scalar eps_;
    Scalar density_ = getParam<Scalar>("Problem.Density", 1.0);
    Scalar constantKinematicViscosity_ = getParam<Scalar>("Problem.KinematicViscosity", 1.0);
    Scalar temperature_ = getParam<Scalar>("Problem.Temperature", 283.15);
    TurbulentConstants turbulentConstants_;
    std::string turbulenceModelName_;
};
} // end namespace Dumux

#endif
