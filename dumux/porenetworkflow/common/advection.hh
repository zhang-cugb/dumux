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
 * \brief This file contains the data which is required to calculate
 *        the fluxes of the pore network model over a face of a finite volume.
 *
 * This means pressure gradients, phase densities at the integration point, etc.
 */
#ifndef DUMUX_PNM_ADVECTIVE_FLUX_HH
#define DUMUX_PNM_ADVECTIVE_FLUX_HH

#include <dumux/discretization/method.hh>

namespace Dumux
{

// forward declaration
template <class ScalarT, class TransmissibilityLaw, DiscretizationMethod Method>
class WashburnAdvection
{};


template<class ScalarT, class TransmissibilityLawType>
class WashburnAdvection<ScalarT, TransmissibilityLawType, DiscretizationMethod::box>
{


public:
    //! Export the Scalar type
    using Scalar = ScalarT;
    using TransmissibilityLaw = TransmissibilityLawType;
    using Cache = typename TransmissibilityLaw::Cache;

    template<class Problem, class Element, class FVElementGeometry,
             class ElementVolumeVariables, class SubControlVolumeFace, class ElemFluxVarsCache>
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const int phaseIdx,
                       const ElemFluxVarsCache& elemFluxVarsCache)
    {
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        // Get the inside and outside volume variables
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[outsideScv];

        // calculate the pressure difference
        const Scalar deltaP = insideVolVars.pressure(phaseIdx) - outsideVolVars.pressure(phaseIdx);
        const Scalar transmissibility = fluxVarsCache.transmissibility(phaseIdx);

        Scalar volumeFlow = transmissibility*deltaP;

        // add gravity term
        static const bool enableGravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");
        if (enableGravity)
        {
            const Scalar rho = 0.5*insideVolVars.density(phaseIdx) + 0.5*outsideVolVars.density(phaseIdx);
            const Scalar g = problem.spatialParams().gravity(scvf.center()) * scvf.unitOuterNormal();

            // The transmissibility is with respect to the effective throat length (potentially dropping the pore body radii).
            // For gravity, we need to consider the total throat length (i.e., the cell-center to cell-center distance).
            // This might cause some inconsistencies TODO: is there a better way?
            volumeFlow += transmissibility * element.geometry().volume() * rho * g;
        }

        return volumeFlow;
    }

    /*!
     * \brief Returns the throat conductivity
     *
     * \param problem The problem
     * \param element The element
     * \param ElementVolumeVariables The element volume variables
     */
    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class FluxVariablesCache>
    static Scalar calculateTransmissibility(const Problem& problem,
                                            const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                            const ElementVolumeVariables& elemVolVars,
                                            const FluxVariablesCache& fluxVarsCache,
                                            const int phaseIdx)
    {
        if constexpr (ElementVolumeVariables::VolumeVariables::numFluidPhases() == 1)
            return TransmissibilityLaw::SinglePhase::singlePhaseTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, fluxVarsCache, phaseIdx);
        else
        {
            static_assert(ElementVolumeVariables::VolumeVariables::numFluidPhases() == 2);

            const auto& spatialParams = problem.spatialParams();
            using FluidSystem = typename ElementVolumeVariables::VolumeVariables::FluidSystem;
            const int wPhaseIdx = spatialParams.template wettingPhase<FluidSystem>(element, elemVolVars);
            const auto eIdx = problem.gridGeometry().gridView().indexSet().index(element);
            const bool invaded = fluxVarsCache.invaded();

            // TODO remove this
            static const bool blockNonWettingPhaseAtOutlet = getParam<Scalar>("SpatialParameters.BlockNonWettingPhaseAtOutlet", true);

            if (phaseIdx == wPhaseIdx)
            {
                return invaded ? TransmissibilityLaw::WettingLayer::transmissibility(element, fvGeometry, scvf, fluxVarsCache)
                               : TransmissibilityLaw::SinglePhase::singlePhaseTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, fluxVarsCache, phaseIdx);
            }
            else // non-wetting phase
            {
                if (!invaded)
                    return 0.0;
                else if (blockNonWettingPhaseAtOutlet && problem.gridGeometry().throatLabel(eIdx) == Labels::outlet)
                    return 0.0;
                else
                    return TransmissibilityLaw::NonWettingPhase::transmissibility(element, fvGeometry, scvf, fluxVarsCache);
            }
        }
    }

    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class FluxVariablesCache>
    static std::array<Scalar, 2> calculateTransmissibilities(const Problem& problem,
                                                             const Element& element,
                                                             const FVElementGeometry& fvGeometry,
                                                             const ElementVolumeVariables& elemVolVars,
                                                             const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                                             const FluxVariablesCache& fluxVarsCache)
    {
        static_assert(ElementVolumeVariables::VolumeVariables::numFluidPhases() == 1);
        const Scalar t = calculateTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, fluxVarsCache, 0);
        return std::array<Scalar, 2>{t, -t};
    }
};


} // end namespace

#endif
