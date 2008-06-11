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
 * \file ParkerLenhardState.hh Specification of the state API for the
 *                           Parker-Lenhard hysteresis model.
 */
#ifndef PARKER_LENHARD_STATE_HH
#define PARKER_LENHARD_STATE_HH

#include <Stupid/Auxilary/Apis.hh>

#include <Stupid/Auxilary/StateHelperMacros.hh>
#include <Stupid/Material/TwophaseSatState.hh>

#include <stdlib.h>

namespace Stupid
{
// forward declaration
template <class ScalarT>
class PLScanningCurve;


namespace Api
{
    BEGIN_API_DEF(ParkerLenhardParams)
    {
        typedef typename Implementation::Scalar Scalar;
        typedef typename Implementation::CapPressureParams PCParams;
        typedef PLScanningCurve<Scalar> ScanningCurve;

        require<TwophaseSatParams>(const_impl);

        Scalar tmp = 0.5;
        tmp = impl.Snrei();

        const PCParams *tmpPc;
        tmpPc = &const_impl.micParams();
        tmpPc = &const_impl.mdcParams();

        const ScanningCurve *tmpSC;
        tmpSC = const_impl.mdc();
        tmpSC = const_impl.pisc();
        tmpSC = const_impl.csc();
    }
    END_API_DEF;

    BEGIN_API_DEF(ParkerLenhardState)
    {
        typedef typename Implementation::Scalar Scalar;
        typedef typename Implementation::CapPressureParams PCParams;
        typedef PLScanningCurve<Scalar> ScanningCurve;

        require<ParkerLenhardParams>(impl);

        Scalar tmp = 0.5;
        impl.setSnrei(tmp);

        ScanningCurve *tmpSC = NULL;
        tmpSC = impl.mdc();   impl.setMdc(tmpSC);
        tmpSC = impl.pisc();  impl.setPisc(tmpSC);
        tmpSC = impl.csc();   impl.setCsc(tmpSC);
    }
    END_API_DEF;
}; // namespace Api


    /*!
     * \brief A reference implementation of the state API for the
     *        Parker-Lenhard hysteresis model.
     */
    template <class CapPressureParamsT>
    class ParkerLenhardState : public TwophaseSatState<typename CapPressureParamsT::Scalar>
    {
    public:
        typedef CapPressureParamsT CapPressureParams;
        typedef typename CapPressureParams::Scalar Scalar;
        typedef Stupid::TwophaseSatState<Scalar> TwophaseSatState;
        typedef Stupid::PLScanningCurve<Scalar>  ScanningCurve;

        ParkerLenhardState(Scalar Swr,
                           Scalar Snr,
                           const CapPressureParams &micParams,
                           const CapPressureParams &mdcParams)
            : TwophaseSatState(Swr, Snr),
              _micParams(micParams),
              _mdcParams(mdcParams)
            {
                Api::require<Api::ParkerLenhardState>(*this);

                _Snrei = 0;
                _mdc = new ScanningCurve();
                _pisc = _csc = NULL;
            }

        ~ParkerLenhardState()
            { delete _mdc; }

        /*!
         * \brief The parameters of the MIC for the capillary pressure
         *        model.
         */
        PROPERTY(CapPressureParams, micParams, setMicParams);

        /*!
         * \brief The parameters of the MDC for the capillary pressure
         *        model.
         */
        PROPERTY(CapPressureParams, mdcParams, setMdcParams);

        // current effective residual saturation
        PROPERTY(Scalar, Snrei, setSnrei);

        /*!
         * \brief The scanning curves
         */
        MUTABLE_PTR_PROPERTY(ScanningCurve, mdc, setMdc); // main drainage
        MUTABLE_PTR_PROPERTY(ScanningCurve, pisc, setPisc); // primary imbibition
        MUTABLE_PTR_PROPERTY(ScanningCurve, csc, setCsc); // current
    };
}; // namespace Stupid

#endif
