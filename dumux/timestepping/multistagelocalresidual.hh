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
 * \brief A wrapper for the local residual to make it usable within a multi-stage scheme
 */
#ifndef DUMUX_TIMESTEPPING_MULTISTAGE_LOCALRESIDUAL_HH
#define DUMUX_TIMESTEPPING_MULTISTAGE_LOCALRESIDUAL_HH

namespace Dumux {

//! forward declaration
template <class Scalar>
class MultiStageParams;

/*!
 * \brief Time stepping with a multi-stage method
 * \note We limit ourselves to "diagonally" implicit multi-stage methods where solving
 *       a stage can only depend on the values of the same stage and stages before
 *       but not future stages (which would require solving larger linear systems)
 */
class<class LocalResidual, class Solutions>
class MultiStageLocalResidual
{
    MultiStageLocalResidual(const LocalResidual& localResidual,
                            const MultiStageParams<Scalar>& p,
                            const Solutions& solutions)
    : localResidual_(localResidual)
    , stageParams_(p)
    , solutions_(solutions)
    {}

    // TODO: elemVolVars, elemFluxVarsCache and bcTypes possibly depend on time
    // and have to evaluated at stageParams_.timeAtStage()
    ElementResidualVector evalResidual(const Element& element,
                                       const FVElementGeometry& fvGeometry)
    {
        ElementResidualVector residual(fvGeometry.numScv());

        for (std::size_t k = 0; k < stageParams_.size(); ++k)
        {
            auto elemVolVars = localView(gridVariables);
            elemVolVars.bind(element, fvGeometry, *solutions_[k]);

            if (!stageParams_.skipTemporal(k))
                residual.axpy(stageParams_.temporalWeight(k), evalTemporal(element, fvGeometry, elemVolVars));

            if (!stageParams_.skipSpatial(k))
            {
                auto elemFluxVarsCache = localView(gridVariables);
                elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

                auto bcTypes = get;

                residual.axpy(stageParams_.spatialWeight(k), evalSpatial(element, fvGeometry, elemVolVars, elemFluxVarsCache, bcTypes));
            }
        }

        return residual;
    }

    ElementResidualVector evalTemporal(const Element& element,
                                       const FVElementGeometry& fvGeometry,
                                       const ElementVolumeVariables& elemVolVars) const
    {
        // initialize the residual vector for all scvs in this element
        ElementResidualVector residual(fvGeometry.numScv());

        // evaluate the volume terms (storage + source terms)
        // forward to the local residual specialized for the discretization methods
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];
            auto storage = localResidual_.computeStorage(problem, scv, volVars);
            storage *= scv.volume()*volVars.extrusionFactor();
            residual[scv.localDofIndex()] += storage;
        }

        return residual;
    }

    ElementResidualVector evalSpatial(const Element& element,
                                      const FVElementGeometry& fvGeometry,
                                      const ElementVolumeVariables& elemVolVars,
                                      const ElementFluxVariablesCache& elemFluxVarsCache,
                                      const ElementBoundaryTypes &bcTypes) const
    {
        // initialize the residual vector for all scvs in this element
        ElementResidualVector residual(fvGeometry.numScv());

        // evaluate the volume terms (storage + source terms)
        // forward to the local residual specialized for the discretization methods
        for (const auto& scv : scvs(fvGeometry))
            localResidual_.evalSource(residual, this->problem(), element, fvGeometry, elemVolVars, scv);

        // forward to the local residual specialized for the discretization methods
        for (const auto& scvf : scvfs(fvGeometry))
            localResidual_.evalFlux(residual, this->problem(), element, fvGeometry, elemVolVars, bcTypes, elemFluxVarsCache, scvf);

        return residual;
    }

private:
    const LocalResidual& localResidual_;
    const TimeStepping::StageParams<Scalar>& stageParams_;
    const Solutions& solutions_;
};

} // end namespace Dumux

#endif
