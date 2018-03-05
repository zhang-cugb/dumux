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
/**
 * \file
 * \ingroup OnePTests
 * \brief Definition of a problem, for the 1p2c problem:
 * Component transport of oxygen in interstitial fluid.
 */
#ifndef DUMUX_TISSUE_PROBLEM_HH
#define DUMUX_TISSUE_PROBLEM_HH

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>

#include <dumux/porousmediumflow/richardsnc/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/liquidphase2c.hh>

#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include "soilspatialparams.hh"

namespace Dumux {

template <class TypeTag>
class SoilProblem;

namespace Properties {

NEW_TYPE_TAG(SoilTypeTag, INHERITS_FROM(CCTpfaModel, RichardsNC));

// Set the grid type
SET_TYPE_PROP(SoilTypeTag, Grid, Dune::UGGrid<3>);

SET_BOOL_PROP(SoilTypeTag, EnableFVGridGeometryCache, true);
SET_BOOL_PROP(SoilTypeTag, EnableGridVolumeVariablesCache, true);
SET_BOOL_PROP(SoilTypeTag, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(SoilTypeTag, SolutionDependentAdvection, false);
SET_BOOL_PROP(SoilTypeTag, SolutionDependentMolecularDiffusion, false);
SET_BOOL_PROP(SoilTypeTag, SolutionDependentHeatConduction, false);

// Set the problem property
SET_TYPE_PROP(SoilTypeTag, Problem, SoilProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(SoilTypeTag, SpatialParams, SoilSpatialParams<TypeTag>);

// Set the fluid system
SET_PROP(SoilTypeTag, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::LiquidPhaseTwoC<Scalar, Components::SimpleH2O<Scalar>,
                                                       Components::Constant<1, Scalar>>;
};

// Set the material law
SET_PROP(SoilTypeTag, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    // define the material law parameterized by absolute saturations
    using type = EffToAbsLaw<RegularizedVanGenuchten<Scalar>>;
};

SET_BOOL_PROP(SoilTypeTag, UseMoles, true);

} // end namespace Properties


/*!
 * \ingroup OnePTests
 */
template <class TypeTag>
class SoilProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PointSource = typename GET_PROP_TYPE(TypeTag, PointSource);

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;

    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);

public:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    enum Indices {
        // world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        pressureIdx = 0,
        transportCompIdx = 1,

        conti0EqIdx = 0,
        transportEqIdx = 1,

        wPhaseIdx = FluidSystem::wPhaseIdx
    };

    SoilProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(fvGridGeometry)
    , couplingManager_(couplingManager)
    {
        //read parameters from input file
        name_ = getParam<std::string>("Problem.Name") + "_3d";
        contaminantMoleFraction_ = getParam<Scalar>("Problem.ContaminantMoleFraction");

        // for initial conditions
        const Scalar sw = getParam<Scalar>("Problem.InitTopSaturation", 0.3); // start with 30% saturation on top
        using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
        pcTop_ = MaterialLaw::pc(this->spatialParams().materialLawParamsAtPos(fvGridGeometry->bBoxMax()), sw);
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Returns the temperature within the domain [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // in [K]

    /*
      * \brief Returns the reference pressure [Pa] of the non-wetting
     *        fluid phase within a finite volume
     *
     * This problem assumes a constant reference pressure of 1 bar.
     */
    Scalar nonWettingReferencePressure() const
    { return 1.0e5; }


    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    { return initialAtPos(globalPos); }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Applies a vector of point sources. The point sources
     *        are possibly solution dependent.
     *
     * \param pointSources A vector of Dumux::PointSource s that contain
              source values for all phases and space positions.
     *
     * For this method, the \a values method of the point source
     * has to return the absolute mass rate in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void addPointSources(std::vector<PointSource>& pointSources) const
    { pointSources = this->couplingManager().bulkPointSources(); }

    /*!
     * \brief Evaluate the point sources (added by addPointSources)
     *        for all phases within a given sub-control-volume.
     *
     * This is the method for the case where the point source is
     * solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param pointSource A single point source
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub-control volume within the element
     *
     * For this method, the \a values() method of the point sources returns
     * the absolute rate mass generated or annihilate in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void pointSource(PointSource& source,
                     const Element &element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolume &scv) const
    {
        NumEqVector sourceValues;

        // compute source at every integration point
        const auto priVars3D = this->couplingManager().bulkPriVars(source.id());
        const auto priVars1D = this->couplingManager().lowDimPriVars(source.id());
        const Scalar pressure3D = priVars3D[pressureIdx];
        const Scalar pressure1D = priVars1D[pressureIdx];

        const auto& spatialParams = this->couplingManager().problem(Dune::index_constant<1>{}).spatialParams();
        const auto lowDimElementIdx = this->couplingManager().pointSourceData(source.id()).lowDimElementIdx();
        const Scalar Kr = spatialParams.Kr(lowDimElementIdx);
        const Scalar rootRadius = spatialParams.radius(lowDimElementIdx);

        // sink defined as radial flow Jr * density [m^2 s-1]* [kg m-3]
        const auto molarDensityH20 = 1000 / 0.018;
        const auto molarDensityD20 = 1000 / 0.020;

        sourceValues[conti0EqIdx] = 2 * M_PI * rootRadius * Kr * (pressure1D - pressure3D) * molarDensityH20;

        const Scalar x3D = priVars3D[transportCompIdx];
        const Scalar x1D = priVars1D[transportCompIdx];

        //! advective transport over root wall
        // compute correct upwind concentration
        if (sourceValues[conti0EqIdx] > 0)
            sourceValues[transportEqIdx] = sourceValues[conti0EqIdx]*x1D/molarDensityH20*molarDensityD20;
        else
            sourceValues[transportEqIdx] = sourceValues[conti0EqIdx]*x3D/molarDensityH20*molarDensityD20;

        //! diffusive transport over root wall
        sourceValues[transportEqIdx] += 2 * M_PI * rootRadius * 1.0e-8 * (x1D - x3D) * molarDensityD20;

        sourceValues *= source.quadratureWeight()*source.integrationElement();
        source = sourceValues;
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        const auto& gg = this->fvGridGeometry();
        const auto xTracer = [&]()
        {
            auto contaminationPos = gg.bBoxMax()-gg.bBoxMin();
            contaminationPos[0] *= 0.25;
            contaminationPos[1] *= 0.55;
            contaminationPos[2] *= 0.25;
            contaminationPos += gg.bBoxMin();

            static const Scalar extend = 0.15*(gg.bBoxMax()[0]-gg.bBoxMin()[0]);
            if ((globalPos - contaminationPos).infinity_norm() <  extend + eps_)
                return contaminantMoleFraction_;
            else
                return 0.0;
        }();

        PrimaryVariables priVars(0.0);
        //! hydrostatic pressure profile
        priVars[pressureIdx] = (nonWettingReferencePressure() - pcTop_)
                                -9.81*1000*(globalPos[dimWorld-1] - gg.bBoxMax()[dimWorld-1]);
        priVars[transportCompIdx] = xTracer;
        return priVars;
    }

    //! Called after every time step
    //! Output the total global exchange term
    void computeSourceIntegral(const SolutionVector& sol, const GridVariables& gridVars)
    {
        PrimaryVariables source(0.0);
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVars.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, sol);

            for (auto&& scv : scvs(fvGeometry))
            {
                auto pointSources = this->scvPointSources(element, fvGeometry, elemVolVars, scv);
                // conversion to kg/s
                using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
                const auto& volVars = elemVolVars[scv];
                pointSources *= scv.volume()*volVars.extrusionFactor()
                                * volVars.density(FluidSystem::wPhaseIdx) / volVars.molarDensity(FluidSystem::wPhaseIdx);

                source += pointSources;
            }
        }

        std::cout << "Global integrated source (soil): " << source[conti0EqIdx] << " (kg/s) / "
                  <<                           source[conti0EqIdx]*3600*24*1000 << " (g/day)" << '\n';
    }

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    Scalar pcTop_, contaminantMoleFraction_;

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} //end namespace Dumux

#endif
