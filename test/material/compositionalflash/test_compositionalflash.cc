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
 * \ingroup MaterialTests
 * \brief This is a program to test the flash calculation routines
 * for compositional sequential models
 *
 * This flash calculation determines the full fluid state of all phases, including
 * pressures, saturations and composition, given the total mass of one component
 * related to the total fluid mass in the pore space.
 */

#include <config.h>

#include <dumux/material/constraintsolvers/misciblemultiphasecomposition.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>
#include <dumux/material/constraintsolvers/compositionalflash.hh>

#include <dumux/material/fluidstates/compositional.hh>

#include <dumux/material/fluidsystems/h2on2.hh>

#include <dumux/material/fluidmatrixinteractions/mp/mpadapter.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedlinearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

template <class Scalar, class FluidState>
void checkSame(const FluidState &fsRef, const FluidState &fsFlash)
{
    enum { numPhases = FluidState::numPhases };
    enum { numComponents = FluidState::numComponents };

    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        Scalar error;

        // check the pressures
        using std::abs;
        error = 1 - fsRef.pressure(phaseIdx)/fsFlash.pressure(phaseIdx);
        if (abs(error) > 1e-6) {
            std::cout << "pressure error phase " << phaseIdx << ": "
                      << fsFlash.pressure(phaseIdx)  << " flash vs "
                      << fsRef.pressure(phaseIdx) << " reference"
                      << " error=" << error << "\n";
        }

        // check the saturations
        error = fsRef.saturation(phaseIdx) - fsFlash.saturation(phaseIdx);
        if (abs(error) > 1e-6)
            std::cout << "saturation error phase " << phaseIdx << ": "
                      << fsFlash.saturation(phaseIdx) << " flash vs "
                      << fsRef.saturation(phaseIdx) << " reference"
                      << " error=" << error << "\n";

        // check the compositions for existing phases
        if(fsRef.saturation(phaseIdx) > 0.0)
            for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
                error = fsRef.moleFraction(phaseIdx, compIdx) - fsFlash.moleFraction(phaseIdx, compIdx);
                if (abs(error) > 1e-6)
                    std::cout << "composition error phase " << phaseIdx << ", component " << compIdx << ": "
                    << fsFlash.moleFraction(phaseIdx, compIdx) << " flash vs "
                    << fsRef.moleFraction(phaseIdx, compIdx) << " reference"
                    << " error=" << error << "\n";
            }
    }
}

template <class Scalar, class FluidSystem, class MaterialLaw, class FluidState>
void checkCompositionalFlash(const FluidState &fsRef)
{
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    using PhaseVector = Dune::FieldVector<Scalar, numPhases>;

    // calculate the total mass of the first component related to the total mass
    // of fluids
    Scalar Z0(0.0);
    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        Z0 +=
        fsRef.phaseMassFraction(phaseIdx)*fsRef.massFraction(phaseIdx, 0);
    }

    PhaseVector phasePressures;
    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        phasePressures[phaseIdx] += fsRef.pressure(phaseIdx);
    }

    // make a new flash
    Dumux::CompositionalFlash<Scalar, FluidSystem> flash;

    // call flash for testing
    FluidState fsFlash;
    flash.concentrationFlash2p2c(fsFlash, Z0, phasePressures, 0/*dummy*/, fsRef.temperature(/*phaseIdx=*/0));

    // compare the "flashed" fluid state with the reference one
    checkSame<Scalar>(fsRef, fsFlash);
}


template <class Scalar, class FluidSystem, class MaterialLaw, class FluidState>
void completeReferenceFluidState(FluidState &fs,
                                 typename MaterialLaw::Params &matParams,
                                 int refPhaseIdx)
{
    enum { numPhases = FluidSystem::numPhases };

    using ComputeFromReferencePhase = Dumux::ComputeFromReferencePhase<Scalar, FluidSystem>;
    using PhaseVector = Dune::FieldVector<Scalar, numPhases>;

    int otherPhaseIdx = 1 - refPhaseIdx;

    // calculate the other saturation
    fs.setSaturation(otherPhaseIdx, 1.0 - fs.saturation(refPhaseIdx));

    // calulate the capillary pressure
    PhaseVector pc;
    MaterialLaw::capillaryPressures(pc, matParams, fs, /*liquidPhaseIdx=*/0);
    fs.setPressure(otherPhaseIdx,
                   fs.pressure(refPhaseIdx)
                   + (pc[otherPhaseIdx] - pc[refPhaseIdx]));

    // make the fluid state consistent with local thermodynamic
    // equilibrium
    typename FluidSystem::ParameterCache paramCache;
    ComputeFromReferencePhase::solve(fs,
                                     paramCache,
                                     refPhaseIdx);
}


int main()
{
    using Scalar = double;
    using FluidSystem = Dumux::FluidSystems::H2ON2<Scalar, Dumux::FluidSystems::H2ON2DefaultPolicy</*fastButSimplifiedRelations=*/true>>;
    using CompositionalFluidState = Dumux::CompositionalFluidState<Scalar, FluidSystem>;

    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    enum { liquidPhaseIdx = FluidSystem::liquidPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { wPhaseIdx = liquidPhaseIdx };

    enum { H2OIdx = FluidSystem::H2OIdx };
    enum { N2Idx = FluidSystem::N2Idx };

    using EffMaterialLaw = Dumux::RegularizedBrooksCorey<Scalar>;
    using MaterialLaw = Dumux::EffToAbsLaw<EffMaterialLaw>;
    using MaterialLawParams = MaterialLaw::Params;
    using MPAdapter = Dumux::MPAdapter<MaterialLaw, numPhases>;

    Scalar T = 273.15 + 25;

    // initialize the tables of the fluid system
    Scalar Tmin = T - 1.0;
    Scalar Tmax = T + 1.0;
    int nT = 3;

    Scalar pmin = 0.0;
    Scalar pmax = 1.25 * 2e6;
    int np = 100;

    FluidSystem::init(Tmin, Tmax, nT, pmin, pmax, np);

    // set the parameters for the capillary pressure law
    MaterialLawParams matParams;
    matParams.setSwr(0.0);
    matParams.setSnr(0.0);
    matParams.setPe(0);
    matParams.setLambda(2.0);

    CompositionalFluidState fsRef;

    // create a fluid state which is consistent

    // set the fluid temperatures
    fsRef.setTemperature(T);

    ////////////////
    // only liquid
    ////////////////
    std::cout << "testing single-phase liquid\n";

    // set liquid saturation
    fsRef.setSaturation(liquidPhaseIdx, 1.0);

    // set pressure of the liquid phase
    fsRef.setPressure(liquidPhaseIdx, 2e5);

    // set the liquid composition to pure water
    fsRef.setMoleFraction(liquidPhaseIdx, N2Idx, 0.0);
    fsRef.setMoleFraction(liquidPhaseIdx, H2OIdx, 1.0);

    fsRef.setWettingPhase(wPhaseIdx);

    // "complete" the fluid state
    completeReferenceFluidState<Scalar, FluidSystem, MPAdapter>(fsRef, matParams, liquidPhaseIdx);

    // check the flash calculation
    checkCompositionalFlash<Scalar, FluidSystem, MPAdapter>(fsRef);

    ////////////////
    // only gas
    ////////////////
    std::cout << "testing single-phase gas\n";
    // set gas saturation
    fsRef.setSaturation(gasPhaseIdx, 1.0);

    // set pressure of the gas phase
    fsRef.setPressure(gasPhaseIdx, 1e6);

    // set the gas composition to 99.9% nitrogen and 0.1% water
    fsRef.setMoleFraction(gasPhaseIdx, N2Idx, 0.999);
    fsRef.setMoleFraction(gasPhaseIdx, H2OIdx, 0.001);

    // "complete" the fluid state
    completeReferenceFluidState<Scalar, FluidSystem, MPAdapter>(fsRef, matParams, gasPhaseIdx);

    // check the flash calculation
    checkCompositionalFlash<Scalar, FluidSystem, MPAdapter>(fsRef);

    ////////////////
    // both phases
    ////////////////
    std::cout << "testing two-phase\n";

    // set saturations
    fsRef.setSaturation(liquidPhaseIdx, 0.5);
    fsRef.setSaturation(gasPhaseIdx, 0.5);

    // set pressures
    fsRef.setPressure(liquidPhaseIdx, 1e6);
    fsRef.setPressure(gasPhaseIdx, 1e6);

    FluidSystem::ParameterCache paramCache;
    using MiscibleMultiPhaseComposition = Dumux::MiscibleMultiPhaseComposition<Scalar, FluidSystem>;
    MiscibleMultiPhaseComposition::solve(fsRef, paramCache);

    // check the flash calculation
    checkCompositionalFlash<Scalar, FluidSystem, MPAdapter>(fsRef);

    ////////////////
    // with capillary pressure
    ////////////////

    MaterialLawParams matParams2;
    matParams2.setSwr(0.0);
    matParams2.setSnr(0.0);
    matParams2.setPe(1e3);
    matParams2.setLambda(2.0);

    // set gas saturation
    fsRef.setSaturation(gasPhaseIdx, 0.5);
    fsRef.setSaturation(liquidPhaseIdx, 0.5);

    // set pressure of the liquid phase
    fsRef.setPressure(liquidPhaseIdx, 1e6);

    // calulate the capillary pressure
    using PhaseVector = Dune::FieldVector<Scalar, numPhases>;
    PhaseVector pc;
    MPAdapter::capillaryPressures(pc, matParams2, fsRef, /*wPhaseIdx=*/0);
    fsRef.setPressure(gasPhaseIdx,
                      fsRef.pressure(liquidPhaseIdx)
                      + (pc[gasPhaseIdx] - pc[liquidPhaseIdx]));

    using MiscibleMultiPhaseComposition = Dumux::MiscibleMultiPhaseComposition<Scalar, FluidSystem>;
    MiscibleMultiPhaseComposition::solve(fsRef, paramCache);


    // check the flash calculation
    checkCompositionalFlash<Scalar, FluidSystem, MPAdapter>(fsRef);

    return 0;
}
