/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file TwophaseSatState.hh Specification of the state API for a twophase
 *                           saturation relation.
 */
#ifndef TWOPHASE_SAT_STATE_HH
#define TWOPHASE_SAT_STATE_HH

#include <Stupid/Auxilary/Apis.hh>

#include <Stupid/Auxilary/StateHelperMacros.hh>

namespace Stupid
{
namespace Api
{
    BEGIN_API_DEF(TwophaseSatParams)
    {
        typedef typename Implementation::Scalar Scalar;
        Scalar tmp = 0.5;
        tmp = const_impl.Swr();
        tmp = const_impl.Snr();
        tmp = const_impl.Snre();
    }
    END_API_DEF;

    BEGIN_API_DEF(TwophaseSatState)
    {
        typedef typename Implementation::Scalar Scalar;
        require<TwophaseSatParams>(impl);

        Scalar tmp = 0.5;
        impl.setSwr(tmp);
        impl.setSnr(tmp);
    }
    END_API_DEF;

}; // namespace Api

    /*!
     * \brief A reference implementation of the state API class for the
     *        twophase saturation relations.
     */
    template <class ScalarT>
    class TwophaseSatState
    {
    public:
        typedef ScalarT Scalar;

        TwophaseSatState(Scalar Swr, Scalar Snr)
            : _Swr(Swr), _Snr(Snr), _Snre(Snr/(1 - _Swr))
            {
                _setSnr(Snr);
            }

        //! Residual saturation of the wetting phase.
        PARAMETER(Scalar, Swr);
        void setSwr(Scalar val) { _Swr = val; _Snre = _Snr/(1 - _Swr); };

        //! Residual saturation of the non-wetting phase.
        PARAMETER(Scalar, Snr);
        void setSnr(Scalar val) { _Snr = val; _Snre = _Snr/(1 - _Swr); }

        //! Effective residual saturation of the non-wetting phase.
        PARAMETER(Scalar, Snre);
    };
}; // namespace Stupid

#endif
