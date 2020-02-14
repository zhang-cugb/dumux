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
 * \brief Test for the Reynolds Averaged Navier-Stokes model with the Komega turbulence model as closure.
 *        Comparison to an analytical solution.
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
 * \file
 * \ingroup RANSTests
 * \brief Test for the Reynolds Averaged Navier-Stokes model with the Komega turbulence model as closure.
 *        Comparison to an analytical solution.
 */
template <class TypeTag>
class AnalyticalRANSProblem : public RANSProblem<TypeTag>
{
    using ParentType = RANSProblem<TypeTag>;

    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    struct TurbulentConstants
    { Scalar alpha, betaK, sigmaK, betaOmega, sigmaOmega; };

    using AnalyticalModel = typename Dumux::KOmegaAnalytical<Scalar, GlobalPosition, PrimaryVariables, ModelTraits, TurbulentConstants, dimWorld>;

public:
    AnalyticalRANSProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry), eps_(1e-6)
    {
        createAnalyticalSolution_();

        turbulenceModelName_ = turbulenceModelToString(ModelTraits::turbulenceModel());
        std::cout << "Using the "<< turbulenceModelName_ << " Turbulence Model. \n";
        std::cout << std::endl;
    }

    /*!
     * \name Problem parameters
     */
    // \{
    bool isOnWallAtPos(const GlobalPosition &globalPos) const
    {
        bool isOnWallAtPos(0);
        isOnWallAtPos = globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_
                         || globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
        if (dimWorld > 1)
            isOnWallAtPos = globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_
                         || globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_
                         || globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_
                         || globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_;
        return isOnWallAtPos;
    }

    /*!
     * \brief Returns the temperature [K] within the domain for the isothermal model.
     */
    Scalar temperature() const
    { return temperature_; }
    // \}

   /*!
     * \name Boundary and Initial Conditions
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
        values.setDirichlet(Indices::velocityXIdx);
        if (dimWorld > 1)
            values.setDirichlet(Indices::velocityYIdx);
        values.setDirichlet(Indices::turbulentKineticEnergyIdx);
        values.setDirichlet(Indices::dissipationIdx);

        return values;
    }

    /*!
     * \brief Returns the sources within the domain.
     *
     * \param globalPos The global position
     */
    PrimaryVariables sourceAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(analyticalSolutionAtPos_(globalPos));
        return values;
    }

    /*!
     * \brief Returns whether a fixed Dirichlet value shall be used at a given cell.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scv The sub control volume
     * \param pvIdx The primary variable index in the solution vector
     */
    bool isDirichletCell(const Element& element,
                         const typename GridGeometry::LocalView& fvGeometry,
                         const typename GridGeometry::SubControlVolume& scv,
                         int pvIdx) const
    {
        bool onBoundary = false;
        for (const auto& scvf : scvfs(fvGeometry))
            onBoundary = std::max(onBoundary, scvf.boundary());
        return onBoundary;
    }

    /*!
     * \brief Return dirichlet boundary values at a given position
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        // use the values of the analytical solution from the initial condition
        return initialAtPos(globalPos);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        // use the values of the analytical solution
        PrimaryVariables values(analyticalSolutionAtPos_(globalPos));
        return values;
    }

    // \}

    /*!
     * \name Convenience Functions
     */
    // \{

    /*!
     * \brief Returns the analytical solution for the Pressure
     */
    auto& getAnalyticalPressureSolution() const
    { return analyticalPressure_; }

    /*!
     * \brief Returns the analytical solution for the Turbulent Kinetic Energy
     */
    auto& getAnalyticalTurbulentKineticEnergySolution() const
    { return analyticalTurbulentKineticEnergy_; }

    /*!
     * \brief Returns the analytical solution for the Dissipation
     */
    auto& getAnalyticalDissipationSolution() const
    { return analyticalDissipiation_; }

    /*!
     * \brief Returns the analytical solution for the Velocity at the cell center
     */
    auto& getAnalyticalVelocitySolution() const
    { return analyticalVelocity_; }

    /*!
     * \brief Returns the analytical solution for the Velocity at the faces
     */
    auto& getAnalyticalVelocitySolutionOnFace() const
    { return analyticalVelocityOnFace_; }

    // \}

private:

    NumEqVector analyticalSolutionAtPos_(const GlobalPosition globalPos) const
    {
        const AnalyticalModel analyticalModel_(density_, constantKinematicViscosity_, turbulentConstants_);
        return analyticalModel_.solution(globalPos);
    }

    void createAnalyticalSolution_()
    {
        turbulentConstants_.alpha = ParentType::alpha();
        turbulentConstants_.betaK = ParentType::betaK();
        turbulentConstants_.sigmaK = ParentType::sigmaK();
        turbulentConstants_.betaOmega = ParentType::betaOmega();
        turbulentConstants_.sigmaOmega = ParentType::sigmaOmega();

        analyticalPressure_.resize(this->gridGeometry().numCellCenterDofs());
        analyticalTurbulentKineticEnergy_.resize(this->gridGeometry().numCellCenterDofs());
        analyticalDissipiation_.resize(this->gridGeometry().numCellCenterDofs());
        analyticalVelocity_.resize(this->gridGeometry().numCellCenterDofs());
        analyticalVelocityOnFace_.resize(this->gridGeometry().numFaceDofs());

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                auto ccDofIdx = scv.dofIndex();
                auto ccDofPosition = scv.dofPosition();
                auto analyticalSolutionAtCc = analyticalSolutionAtPos_(ccDofPosition);

                // Pressure at the CC
                analyticalPressure_[ccDofIdx] = analyticalSolutionAtCc[Indices::pressureIdx];

                // Turbulence Variables at the CC
                analyticalPressure_[ccDofIdx] = analyticalSolutionAtCc[Indices::turbulentKineticEnergyIdx];
                analyticalPressure_[ccDofIdx] = analyticalSolutionAtCc[Indices::dissipationIdx];

                // velocities at the cell center
                for(int dirIdx = 0; dirIdx < ModelTraits::dim(); ++dirIdx)
                    analyticalVelocity_[ccDofIdx][dirIdx] = analyticalSolutionAtCc[Indices::velocity(dirIdx)];

                // Velocities on the faces
                for (auto&& scvf : scvfs(fvGeometry))
                {
                    const auto faceDofIdx = scvf.dofIndex();
                    const auto faceDofPosition = scvf.center();
                    const auto dirIdx = scvf.directionIndex();
                    const auto analyticalSolutionAtFace = analyticalSolutionAtPos_(faceDofPosition);
                    analyticalVelocityOnFace_[faceDofIdx][dirIdx] = analyticalSolutionAtFace[Indices::velocity(dirIdx)];
                }
            }
        }
    }

    Scalar eps_;
    Scalar density_ = getParam<Scalar>("Component.LiquidDensity");
    Scalar constantKinematicViscosity_ = getParam<Scalar>("Component.LiquidKinematicViscosity");
    Scalar temperature_ = getParam<Scalar>("Component.Temperature");
    TurbulentConstants turbulentConstants_;

    std::vector<Scalar> analyticalPressure_;
    std::vector<Scalar> analyticalTurbulentKineticEnergy_;
    std::vector<Scalar> analyticalDissipiation_;
    std::vector<VelocityVector> analyticalVelocity_;
    std::vector<VelocityVector> analyticalVelocityOnFace_;

    std::string turbulenceModelName_;
};
} // end namespace Dumux

#endif
