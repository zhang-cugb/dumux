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

//! forward declarations
template<class Scalar> class MultiStageParams;

/*!
 * \brief Time stepping with a multi-stage method
 * \note We limit ourselves to "diagonally" implicit multi-stage methods where solving
 *       a stage can only depend on the values of the same stage and stages before
 *       but not future stages (which would require solving larger linear systems)
 */
template<class LocalResidual>
class MultiStageLocalResidual
{
public:
    //! export the types underlying the local residual we wrap around here
    using Variables = typename LocalResidual::Variables;
    using Residual = typename LocalResidual::Residual;
    using Scalar = typename LocalResidual::Scalar;

    /*!
     * \brief The constructor
     * \param localResiduals The local residuals required for this stage
     * \param params The parameters of this stage
     */
    MultiStageLocalResidual(std::vector<LocalResidual>& localResiduals,
                            const MultiStageParams<Scalar>& params)
    : localResiduals_(localResiduals)
    , stageParams_(params)
    {
        if (localResiduals.empty())
            DUNE_THROW(Dune::InvalidStateException, "At least one residual is required");
        if (localResiduals.size() != params.size())
            DUNE_THROW(Dune::InvalidStateException, "Size mismatch between residuals and stage params");
    }

    Residual eval()
    {
        // The residual returned by local residuals may be be a resizable
        // vector or simply a scalar. Thus, we require a way to construct
        // an empty residual!? This could be put elsewhere...
        Residual residual = localResiduals_[0].getEmptyResidual();

        for (std::size_t k = 0; k < stageParams_.size(); ++k)
        {
            if (!stageParams_.skipTemporal(k))
            { residual.axpy(stageParams_.temporalWeight(k), localResiduals_[k].evalStorage()); }

            if (!stageParams_.skipSpatial(k))
            { residual.axpy(stageParams_.spatialWeight(k), localResiduals_[k].evalFluxesAndSources()); }
        }

        return residual;
    }

    //! TODO: SHOULD THIS HAVE EVALTEMPORAL() AND EVALSPATIAL() AS WELL?
    //!       I.E. SHOULD IT BEHAVE LIKE A STANDARD LOCALRESIDUAL?

private:
    std::vector<LocalResidual>& localResiduals_;
    const MultiStageParams<Scalar>& stageParams_;
};

} // end namespace Dumux

#endif
