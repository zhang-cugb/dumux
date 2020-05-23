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
 * \brief A wrapper for local operators to make it usable within a multi-stage
 *        time integration schemes.
 */
#ifndef DUMUX_TIMESTEPPING_MULTISTAGE_LOCALOPERATOR_HH
#define DUMUX_TIMESTEPPING_MULTISTAGE_LOCALOPERATOR_HH

namespace Dumux {

//! forward declarations
template<class Scalar> class MultiStageParams;

/*!
 * \brief Time stepping with a multi-stage method
 * \note We limit ourselves to "diagonally" implicit multi-stage methods where solving
 *       a stage can only depend on the values of the same stage and stages before
 *       but not future stages (which would require solving larger linear systems)
 */
template<class LocalOperator>
class MultiStageLocalOperator
{
public:
    //! export the types underlying the local residual we wrap around here
    // TODO: Exporting Residual in local operators doesn't make much sense.
    //       Can we find something better!? They evaluate the terms of the
    //       equation but not the entire residual. A more generic name for
    //       the container type would be good.
    using Variables = typename LocalOperator::Variables;
    using Residual = typename LocalOperator::Residual;
    using Scalar = typename LocalOperator::Scalar;

    /*!
     * \brief The constructor
     * \param localResiduals The local residuals required for this stage
     * \param params The parameters of this stage
     */
    MultiStageLocalOperator(std::vector<LocalOperator>& localOperators,
                            const MultiStageParams<Scalar>& params)
    : localOperators_(localOperators)
    , stageParams_(params)
    {
        if (localOperators.empty())
            DUNE_THROW(Dune::InvalidStateException, "At least one local operator is required");
        if (localOperators.size() != params.size())
            DUNE_THROW(Dune::InvalidStateException, "Size mismatch between operators and stage params");
    }

    //! evaluates the entire local residual
    Residual evalLocalResidual()
    {
        // The residual returned by local residuals may be be a resizable
        // vector or simply a scalar. Thus, we require a way to construct
        // an empty residual!? This could be put elsewhere...
        Residual residual = localOperators_[0].getEmptyResidual();

        for (std::size_t k = 0; k < stageParams_.size(); ++k)
        {
            if (!stageParams_.skipTemporal(k))
            { residual.axpy(stageParams_.temporalWeight(k), localOperators_[k].evalStorage()); }

            if (!stageParams_.skipSpatial(k))
            { residual.axpy(stageParams_.spatialWeight(k), localOperators_[k].evalFluxesAndSources()); }
        }

        return residual;
    }

    Residual evalStorage()
    {
        // The residual returned by local residuals may be be a resizable
        // vector or simply a scalar. Thus, we require a way to construct
        // an empty residual!? This could be put elsewhere...
        Residual result = localOperators_[0].getEmptyResidual();

        for (std::size_t k = 0; k < stageParams_.size(); ++k)
        {
            if (!stageParams_.skipTemporal(k))
            { result.axpy(stageParams_.temporalWeight(k), localOperators_[k].evalStorage()); }
        }

        return result;
    }

    Residual evalFluxesAndSources()
    {
        // The residual returned by local residuals may be be a resizable
        // vector or simply a scalar. Thus, we require a way to construct
        // an empty residual!? This could be put elsewhere...
        Residual result = localOperators_[0].getEmptyResidual();

        for (std::size_t k = 0; k < stageParams_.size(); ++k)
        {
            if (!stageParams_.skipSpatial(k))
            { result.axpy(stageParams_.spatialWeight(k), localOperators_[k].evalFluxesAndSources()); }
        }

        return result;
    }

    //! TODO: SHOULD THIS HAVE EVALTEMPORAL() AND EVALSPATIAL() AS WELL?
    //!       I.E. SHOULD IT BEHAVE LIKE A STANDARD LOCALRESIDUAL?

private:
    std::vector<LocalOperator>& localOperators_;
    const MultiStageParams<Scalar>& stageParams_;
};

} // end namespace Dumux

#endif
