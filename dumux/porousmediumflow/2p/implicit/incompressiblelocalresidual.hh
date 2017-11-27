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
 *
 * \brief Element-wise calculation of the residual and its derivatives
 *        for a two-phase, incompressible test problem.
 */
#ifndef DUMUX_2P_INCOMPRESSIBLE_TEST_LOCAL_RESIDUAL_HH
#define DUMUX_2P_INCOMPRESSIBLE_TEST_LOCAL_RESIDUAL_HH

#include <dumux/discretization/methods.hh>
#include <dumux/porousmediumflow/immiscible/localresidual.hh>

namespace Dumux
{

template<class TypeTag>
class TwoPIncompressibleLocalResidual : public ImmiscibleLocalResidual<TypeTag>
{
    using ParentType = ImmiscibleLocalResidual<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementResidualVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using EnergyLocalResidual = typename GET_PROP_TYPE(TypeTag, EnergyLocalResidual);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    // first index for the mass balance
    enum
    {
        contiWEqIdx = Indices::contiWEqIdx,
        contiNEqIdx = Indices::contiNEqIdx,

        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx
    };

    static const int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);

public:
    using ParentType::ParentType;

    template<class PartialDerivativeMatrix>
    void addStorageDerivatives(PartialDerivativeMatrix& partialDerivatives,
                               const Problem& problem,
                               const Element& element,
                               const FVElementGeometry& fvGeometry,
                               const VolumeVariables& curVolVars,
                               const SubControlVolume& scv) const
    {
        static_assert(!FluidSystem::isCompressible(FluidSystem::wPhaseIdx),
                      "2p/incompressiblelocalresidual.hh: Only incompressible fluids are allowed!");
        static_assert(!FluidSystem::isCompressible(FluidSystem::nPhaseIdx),
                      "2p/incompressiblelocalresidual.hh: Only incompressible fluids are allowed!");

        // we know that these values are constant throughout the simulation
        static const auto phi = curVolVars.porosity();
        static const auto phi_rho_w = phi*curVolVars.density(FluidSystem::wPhaseIdx);
        static const auto phi_rho_n = phi*curVolVars.density(FluidSystem::nPhaseIdx);

        const auto volume = scv.volume();

        // partial derivative of wetting phase storage term w.r.t. p_w
        partialDerivatives[contiWEqIdx][pressureIdx] += 0.0;
        // partial derivative of wetting phase storage term w.r.t. S_n
        partialDerivatives[contiWEqIdx][saturationIdx] -= volume*phi_rho_w/this->timeLoop().timeStepSize();
        // partial derivative of non-wetting phase storage term w.r.t. p_w
        partialDerivatives[contiNEqIdx][pressureIdx] += 0.0;
        // partial derivative of non-wetting phase storage term w.r.t. S_n
        partialDerivatives[contiNEqIdx][saturationIdx] += volume*phi_rho_n/this->timeLoop().timeStepSize();
    }

    template<class PartialDerivativeMatrix>
    void addSourceDerivatives(PartialDerivativeMatrix& partialDerivatives,
                              const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const VolumeVariables& curVolVars,
                              const SubControlVolume& scv) const
    { /* TODO maybe forward to problem for the user to implement the source derivatives?*/ }

    template<class PartialDerivativeMatrices, class T = TypeTag>
    std::enable_if_t<GET_PROP_VALUE(T, DiscretizationMethod) != DiscretizationMethods::Box, void>
    addFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                       const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& curElemVolVars,
                       const ElementFluxVariablesCache& elemFluxVarsCache,
                       const SubControlVolumeFace& scvf) const
    {
        static_assert(!FluidSystem::isCompressible(FluidSystem::wPhaseIdx),
                      "2p/incompressiblelocalresidual.hh: Only incompressible fluids are allowed!");
        static_assert(!FluidSystem::isCompressible(FluidSystem::nPhaseIdx),
                      "2p/incompressiblelocalresidual.hh: Only incompressible fluids are allowed!");
        static_assert(FluidSystem::viscosityIsConstant(FluidSystem::wPhaseIdx),
                      "2p/incompressiblelocalresidual.hh: Only fluids with constant viscosities are allowed!");
        static_assert(FluidSystem::viscosityIsConstant(FluidSystem::nPhaseIdx),
                      "2p/incompressiblelocalresidual.hh: Only fluids with constant viscosities are allowed!");

        using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
        using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);

        // evaluate the current wetting phase Darcy flux and resulting upwind weights
        static const Scalar upwindWeight = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, UpwindWeight);
        const auto flux_w = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars,
                                                scvf, FluidSystem::wPhaseIdx, elemFluxVarsCache);
        const auto flux_n = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars,
                                                scvf, FluidSystem::nPhaseIdx, elemFluxVarsCache);
        const auto insideWeight_w = std::signbit(flux_w) ? (1.0 - upwindWeight) : upwindWeight;
        const auto outsideWeight_w = 1.0 - insideWeight_w;
        const auto insideWeight_n = std::signbit(flux_n) ? (1.0 - upwindWeight) : upwindWeight;
        const auto outsideWeight_n = 1.0 - insideWeight_n;

        // get references to the two participating vol vars & parameters
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto outsideScvIdx = scvf.outsideScvIdx();
        const auto outsideElement = fvGeometry.fvGridGeometry().element(outsideScvIdx);
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
        const auto& insideVolVars = curElemVolVars[insideScvIdx];
        const auto& outsideVolVars = curElemVolVars[outsideScvIdx];
        const auto& insideMaterialParams = problem.spatialParams().materialLawParams(element,
                                                                                     insideScv,
                                                                                     ElementResidualVector({insideVolVars.priVars()}));
        const auto& outsideMaterialParams = problem.spatialParams().materialLawParams(outsideElement,
                                                                                      outsideScv,
                                                                                      ElementResidualVector({outsideVolVars.priVars()}));

        // get references to the two participating derivative matrices
        auto& dI_dI = derivativeMatrices[insideScvIdx];
        auto& dI_dJ = derivativeMatrices[outsideScvIdx];

        // some quantities to be reused (rho & mu are constant and thus equal for all cells)
        static const auto rho_w = insideVolVars.density(FluidSystem::wPhaseIdx);
        static const auto rho_n = insideVolVars.density(FluidSystem::nPhaseIdx);
        static const auto rhow_muw = rho_w/insideVolVars.viscosity(FluidSystem::wPhaseIdx);
        static const auto rhon_mun = rho_n/insideVolVars.viscosity(FluidSystem::nPhaseIdx);
        const auto rhowKrw_muw_inside = rho_w*insideVolVars.mobility(FluidSystem::wPhaseIdx);
        const auto rhonKrn_mun_inside = rho_n*insideVolVars.mobility(FluidSystem::nPhaseIdx);
        const auto rhowKrw_muw_outside = rho_w*outsideVolVars.mobility(FluidSystem::wPhaseIdx);
        const auto rhonKrn_mun_outside = rho_n*outsideVolVars.mobility(FluidSystem::nPhaseIdx);

        // derivative w.r.t. to Sn is the negative of the one w.r.t. Sw
        const auto insideSw = insideVolVars.saturation(FluidSystem::wPhaseIdx);
        const auto outsideSw = outsideVolVars.saturation(FluidSystem::wPhaseIdx);
        const auto dKrw_dSn_inside = MaterialLaw::dkrw_dsw(insideMaterialParams, insideSw);
        const auto dKrw_dSn_outside = MaterialLaw::dkrw_dsw(outsideMaterialParams, outsideSw);
        const auto dKrn_dSn_inside = MaterialLaw::dkrn_dsw(insideMaterialParams, insideSw);
        const auto dKrn_dSn_outside = MaterialLaw::dkrn_dsw(outsideMaterialParams, outsideSw);
        const auto dpc_dSn_inside = MaterialLaw::dpc_dsw(insideMaterialParams, insideSw);
        const auto dpc_dSn_outside = MaterialLaw::dpc_dsw(outsideMaterialParams, outsideSw);

        const auto tij = elemFluxVarsCache[scvf].advectionTij();

        // precalculate values
        const auto up_w = rhowKrw_muw_inside*insideWeight_w + rhowKrw_muw_outside*outsideWeight_w;
        const auto up_n = rhonKrn_mun_inside*insideWeight_n + rhonKrn_mun_outside*outsideWeight_n;
        const auto rho_mu_flux_w = rhow_muw*flux_w;
        const auto rho_mu_flux_n = rhon_mun*flux_n;
        const auto tij_up_w = tij*up_w;
        const auto tij_up_n = tij*up_n;

        // partial derivative of the wetting phase flux w.r.t. p_w
        dI_dI[contiWEqIdx][pressureIdx] += tij_up_w;
        dI_dJ[contiWEqIdx][pressureIdx] -= tij_up_w;

        // partial derivative of the wetting phase flux w.r.t. S_n
        dI_dI[contiWEqIdx][saturationIdx] -= rho_mu_flux_w*dKrw_dSn_inside*insideWeight_w;
        dI_dJ[contiWEqIdx][saturationIdx] -= rho_mu_flux_w*dKrw_dSn_outside*outsideWeight_w;

        // partial derivative of the non-wetting phase flux w.r.t. p_w
        dI_dI[contiNEqIdx][pressureIdx] += tij_up_n;
        dI_dJ[contiNEqIdx][pressureIdx] -= tij_up_n;

        // partial derivative of the non-wetting phase flux w.r.t. S_n (relative permeability derivative contribution)
        dI_dI[contiNEqIdx][saturationIdx] -= rho_mu_flux_n*dKrn_dSn_inside*insideWeight_n;
        dI_dJ[contiNEqIdx][saturationIdx] -= rho_mu_flux_n*dKrn_dSn_outside*outsideWeight_n;

        // partial derivative of the non-wetting phase flux w.r.t. S_n (capillary pressure derivative contribution)
        dI_dI[contiNEqIdx][saturationIdx] -= tij_up_n*dpc_dSn_inside;
        dI_dJ[contiNEqIdx][saturationIdx] += tij_up_n*dpc_dSn_outside;
    }

    template<class JacobianMatrix, class T = TypeTag>
    std::enable_if_t<GET_PROP_VALUE(T, DiscretizationMethod) == DiscretizationMethods::Box, void>
    addFluxDerivatives(JacobianMatrix& A,
                       const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& curElemVolVars,
                       const ElementFluxVariablesCache& elemFluxVarsCache,
                       const SubControlVolumeFace& scvf) const
    {
        static_assert(!FluidSystem::isCompressible(FluidSystem::wPhaseIdx),
                      "2p/incompressiblelocalresidual.hh: Only incompressible fluids are allowed!");
        static_assert(!FluidSystem::isCompressible(FluidSystem::nPhaseIdx),
                      "2p/incompressiblelocalresidual.hh: Only incompressible fluids are allowed!");
        static_assert(FluidSystem::viscosityIsConstant(FluidSystem::wPhaseIdx),
                      "2p/incompressiblelocalresidual.hh: Only fluids with constant viscosities are allowed!");
        static_assert(FluidSystem::viscosityIsConstant(FluidSystem::nPhaseIdx),
                      "2p/incompressiblelocalresidual.hh: Only fluids with constant viscosities are allowed!");

        using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
        using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);

        // evaluate the current wetting phase Darcy flux and resulting upwind weights
        static const Scalar upwindWeight = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, UpwindWeight);
        const auto flux_w = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars,
                                                scvf, FluidSystem::wPhaseIdx, elemFluxVarsCache);
        const auto flux_n = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars,
                                                scvf, FluidSystem::nPhaseIdx, elemFluxVarsCache);
        const auto insideWeight_w = std::signbit(flux_w) ? (1.0 - upwindWeight) : upwindWeight;
        const auto outsideWeight_w = 1.0 - insideWeight_w;
        const auto insideWeight_n = std::signbit(flux_n) ? (1.0 - upwindWeight) : upwindWeight;
        const auto outsideWeight_n = 1.0 - insideWeight_n;

        // get references to the two participating vol vars & parameters
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto outsideScvIdx = scvf.outsideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
        const auto& insideVolVars = curElemVolVars[insideScv];
        const auto& outsideVolVars = curElemVolVars[outsideScv];

        // we need the element solution for the material parameters
        ElementResidualVector elemSol(fvGeometry.numScv());
        for (const auto& scv : scvs(fvGeometry)) elemSol[scv.indexInElement()] = curElemVolVars[scv].priVars();
        const auto& insideMaterialParams = problem.spatialParams().materialLawParams(element, insideScv, elemSol);
        const auto& outsideMaterialParams = problem.spatialParams().materialLawParams(element, outsideScv, elemSol);

        // some quantities to be reused (rho & mu are constant and thus equal for all cells)
        static const auto rho_w = insideVolVars.density(FluidSystem::wPhaseIdx);
        static const auto rho_n = insideVolVars.density(FluidSystem::nPhaseIdx);
        static const auto rhow_muw = rho_w/insideVolVars.viscosity(FluidSystem::wPhaseIdx);
        static const auto rhon_mun = rho_n/insideVolVars.viscosity(FluidSystem::nPhaseIdx);
        const auto rhowKrw_muw_inside = rho_w*insideVolVars.mobility(FluidSystem::wPhaseIdx);
        const auto rhonKrn_mun_inside = rho_n*insideVolVars.mobility(FluidSystem::nPhaseIdx);
        const auto rhowKrw_muw_outside = rho_w*outsideVolVars.mobility(FluidSystem::wPhaseIdx);
        const auto rhonKrn_mun_outside = rho_n*outsideVolVars.mobility(FluidSystem::nPhaseIdx);

        // let the Law for the advective fluxes calculate the transmissibilities
        const auto ti = AdvectionType::calculateTransmissibilities(problem,
                                                                   element,
                                                                   fvGeometry,
                                                                   curElemVolVars,
                                                                   scvf,
                                                                   elemFluxVarsCache[scvf]);

        // get the rows of the jacobian matrix for the inside/outside scv
        auto& dI_dJ_inside = A[insideScv.dofIndex()];
        auto& dI_dJ_outside = A[outsideScv.dofIndex()];

        // precalculate values
        const auto up_w = rhowKrw_muw_inside*insideWeight_w + rhowKrw_muw_outside*outsideWeight_w;
        const auto up_n = rhonKrn_mun_inside*insideWeight_n + rhonKrn_mun_outside*outsideWeight_n;
        const auto rho_mu_flux_w = rhow_muw*flux_w;
        const auto rho_mu_flux_n = rhon_mun*flux_n;

        // add the partial derivatives w.r.t all scvs in the element
        for (const auto& scvJ : scvs(fvGeometry))
        {
            const auto globalJ = scvJ.dofIndex();
            const auto localJ = scvJ.indexInElement();

            // the transmissibily associated with the scvJ
            const auto tj = ti[localJ];

            // partial derivative of the wetting phase flux w.r.t. p_w
            const auto tj_up_w = tj*up_w;
            dI_dJ_inside[globalJ][contiWEqIdx][pressureIdx] += tj_up_w;
            dI_dJ_outside[globalJ][contiWEqIdx][pressureIdx] -= tj_up_w;

            // partial derivative of the non-wetting phase flux w.r.t. p_w
            const auto tj_up_n = tj*up_n;
            dI_dJ_inside[globalJ][contiNEqIdx][pressureIdx] += tj_up_n;
            dI_dJ_outside[globalJ][contiNEqIdx][pressureIdx] -= tj_up_n;

            // partial derivatives w.r.t. S_n (are the negative of those w.r.t sw)
            // relative permeability contributions only for inside/outside
            if (localJ == insideScvIdx)
            {
                // partial derivative of the wetting phase flux w.r.t. S_n
                const auto insideSw = insideVolVars.saturation(FluidSystem::wPhaseIdx);
                const auto dKrw_dSn_inside = MaterialLaw::dkrw_dsw(insideMaterialParams, insideSw);
                const auto dFluxW_dSnJ = rho_mu_flux_w*dKrw_dSn_inside*insideWeight_w;
                dI_dJ_inside[globalJ][contiWEqIdx][saturationIdx] -= dFluxW_dSnJ;
                dI_dJ_outside[globalJ][contiWEqIdx][saturationIdx] += dFluxW_dSnJ;

                // partial derivative of the non-wetting phase flux w.r.t. S_n (k_rn contribution)
                const auto dKrn_dSn_inside = MaterialLaw::dkrn_dsw(insideMaterialParams, insideSw);
                const auto dFluxN_dSnJ_krn = rho_mu_flux_n*dKrn_dSn_inside*insideWeight_n;
                dI_dJ_inside[globalJ][contiNEqIdx][saturationIdx] -= dFluxN_dSnJ_krn;
                dI_dJ_outside[globalJ][contiWEqIdx][saturationIdx] += dFluxN_dSnJ_krn;

                // partial derivative of the non-wetting phase flux w.r.t. S_n (p_c contribution)
                const auto dFluxN_dSnJ_pc = tj_up_n*MaterialLaw::dpc_dsw(insideMaterialParams, insideSw);
                dI_dJ_inside[globalJ][contiNEqIdx][saturationIdx] -= dFluxN_dSnJ_pc;
                dI_dJ_outside[globalJ][contiNEqIdx][saturationIdx] += dFluxN_dSnJ_pc;
            }
            else if (localJ == outsideScvIdx)
            {
                // see comments for (globalJ == insideScvIdx)
                const auto outsideSw = outsideVolVars.saturation(FluidSystem::wPhaseIdx);
                const auto dKrw_dSn_outside = MaterialLaw::dkrw_dsw(outsideMaterialParams, outsideSw);
                const auto dFluxW_dSnJ = rho_mu_flux_w*dKrw_dSn_outside*outsideWeight_w;
                dI_dJ_inside[globalJ][contiWEqIdx][saturationIdx] -= dFluxW_dSnJ;
                dI_dJ_outside[globalJ][contiWEqIdx][saturationIdx] += dFluxW_dSnJ;

                const auto dKrn_dSn_outside = MaterialLaw::dkrn_dsw(outsideMaterialParams, outsideSw);
                const auto dFluxN_dSnJ_krn = rho_mu_flux_n*dKrn_dSn_outside*outsideWeight_n;
                dI_dJ_inside[globalJ][contiNEqIdx][saturationIdx] -= dFluxN_dSnJ_krn;
                dI_dJ_outside[globalJ][contiNEqIdx][saturationIdx] += dFluxN_dSnJ_krn;

                const auto dFluxN_dSnJ_pc = tj_up_n*MaterialLaw::dpc_dsw(outsideMaterialParams, outsideSw);
                dI_dJ_inside[globalJ][contiNEqIdx][saturationIdx] -= dFluxN_dSnJ_pc;
                dI_dJ_outside[globalJ][contiNEqIdx][saturationIdx] += dFluxN_dSnJ_pc;
            }
            else
            {
                const auto& paramsJ = problem.spatialParams().materialLawParams(element, scvJ, elemSol);
                const auto swJ = curElemVolVars[scvJ].saturation(FluidSystem::wPhaseIdx);
                const auto dFluxN_dSnJ_pc = tj_up_n*MaterialLaw::dpc_dsw(paramsJ, swJ);
                dI_dJ_inside[globalJ][contiNEqIdx][saturationIdx] -= dFluxN_dSnJ_pc;
                dI_dJ_outside[globalJ][contiNEqIdx][saturationIdx] += dFluxN_dSnJ_pc;
            }
        }
    }

    template<class PartialDerivativeMatrices>
    void addCCDirichletFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                                       const Problem& problem,
                                       const Element& element,
                                       const FVElementGeometry& fvGeometry,
                                       const ElementVolumeVariables& curElemVolVars,
                                       const ElementFluxVariablesCache& elemFluxVarsCache,
                                       const SubControlVolumeFace& scvf) const
    {
        using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
        using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);

        // evaluate the current wetting phase Darcy flux and resulting upwind weights
        static const Scalar upwindWeight = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, UpwindWeight);
        const auto flux_w = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars,
                                                     scvf, FluidSystem::wPhaseIdx, elemFluxVarsCache);
        const auto flux_n = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars,
                                                        scvf, FluidSystem::nPhaseIdx, elemFluxVarsCache);
        const auto insideWeight_w = std::signbit(flux_w) ? (1.0 - upwindWeight) : upwindWeight;
        const auto outsideWeight_w = 1.0 - insideWeight_w;
        const auto insideWeight_n = std::signbit(flux_n) ? (1.0 - upwindWeight) : upwindWeight;
        const auto outsideWeight_n = 1.0 - insideWeight_n;

        // get references to the two participating vol vars & parameters
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = curElemVolVars[insideScvIdx];
        const auto& outsideVolVars = curElemVolVars[scvf.outsideScvIdx()];
        const auto& insideMaterialParams = problem.spatialParams().materialLawParams(element,
                                                                                     insideScv,
                                                                                     ElementResidualVector({insideVolVars.priVars()}));

        // some quantities to be reused (rho & mu are constant and thus equal for all cells)
        static const auto rho_w = insideVolVars.density(FluidSystem::wPhaseIdx);
        static const auto rho_n = insideVolVars.density(FluidSystem::nPhaseIdx);
        static const auto rhow_muw = rho_w/insideVolVars.viscosity(FluidSystem::wPhaseIdx);
        static const auto rhon_mun = rho_n/insideVolVars.viscosity(FluidSystem::nPhaseIdx);
        const auto rhowKrw_muw_inside = rho_w*insideVolVars.mobility(FluidSystem::wPhaseIdx);
        const auto rhonKrn_mun_inside = rho_n*insideVolVars.mobility(FluidSystem::nPhaseIdx);
        const auto rhowKrw_muw_outside = rho_w*outsideVolVars.mobility(FluidSystem::wPhaseIdx);
        const auto rhonKrn_mun_outside = rho_n*outsideVolVars.mobility(FluidSystem::nPhaseIdx);

        // get reference to the inside derivative matrix
        auto& dI_dI = derivativeMatrices[insideScvIdx];

        // derivative w.r.t. to Sn is the negative of the one w.r.t. Sw
        const auto insideSw = insideVolVars.saturation(FluidSystem::wPhaseIdx);
        const auto dKrw_dSn_inside = -1.0*MaterialLaw::dkrw_dsw(insideMaterialParams, insideSw);
        const auto dKrn_dSn_inside = -1.0*MaterialLaw::dkrn_dsw(insideMaterialParams, insideSw);
        const auto dpc_dSn_inside = -1.0*MaterialLaw::dpc_dsw(insideMaterialParams, insideSw);

        const auto tij = elemFluxVarsCache[scvf].advectionTij();
        // partial derivative of the wetting phase flux w.r.t. p_w
        const auto up_w = rhowKrw_muw_inside*insideWeight_w + rhowKrw_muw_outside*outsideWeight_w;
        dI_dI[contiWEqIdx][pressureIdx] += tij*up_w;

        // partial derivative of the wetting phase flux w.r.t. S_n
        dI_dI[contiWEqIdx][saturationIdx] += rhow_muw*flux_w*dKrw_dSn_inside*insideWeight_w;

        // partial derivative of the non-wetting phase flux w.r.t. p_w
        const auto up_n = rhonKrn_mun_inside*insideWeight_n + rhonKrn_mun_outside*outsideWeight_n;
        dI_dI[contiNEqIdx][pressureIdx] += tij*up_n;

        // partial derivative of the non-wetting phase flux w.r.t. S_n (relative permeability derivative contribution)
        dI_dI[contiNEqIdx][saturationIdx] += rhon_mun*flux_n*dKrn_dSn_inside*insideWeight_n;

        // partial derivative of the non-wetting phase flux w.r.t. S_n (capillary pressure derivative contribution)
        dI_dI[contiNEqIdx][saturationIdx] += tij*dpc_dSn_inside*up_n;
    }

    template<class PartialDerivativeMatrices>
    void addRobinFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                                 const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& curElemVolVars,
                                 const ElementFluxVariablesCache& elemFluxVarsCache,
                                 const SubControlVolumeFace& scvf) const
    { /* TODO maybe forward to problem for the user to implement the Robin derivatives?*/ }
};

} // end namespace Dumux

#endif
