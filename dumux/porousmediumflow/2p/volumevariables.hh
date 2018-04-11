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
 * \ingroup TwoPModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase model.
 */
#ifndef DUMUX_2P_VOLUME_VARIABLES_HH
#define DUMUX_2P_VOLUME_VARIABLES_HH

#include <dumux/porousmediumflow/volumevariables.hh>
#include "indices.hh"

namespace Dumux {

/*!
 * \ingroup TwoPModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase model.
 */
template <class Traits>
class TwoPVolumeVariables
: public PorousMediumFlowVolumeVariables<Traits, TwoPVolumeVariables<Traits> >
{
    using ParentType = PorousMediumFlowVolumeVariables<Traits, TwoPVolumeVariables<Traits> >;

    using PermeabilityType = typename Traits::PermeabilityType;
    using ModelTraits = typename Traits::ModelTraits;
    using Indices = typename ModelTraits::Indices;
    static constexpr auto formulation = ModelTraits::priVarFormulation();
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using FS = typename Traits::FluidSystem;

    enum
    {
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,

        phase0Idx = FS::phase0Idx,
        phase1Idx = FS::phase1Idx
    };


public:
    //! export type of fluid system
    using FluidSystem = typename Traits::FluidSystem;
    //! export type of fluid state
    using FluidState = typename Traits::FluidState;

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub control volume
    */
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol &elemSol,
                const Problem &problem,
                const Element &element,
                const Scv& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        completeFluidState(elemSol, problem, element, scv, fluidState_);

        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);

        const int wPhaseIdx = problem.spatialParams().template wettingPhase<FluidSystem>(element, scv, elemSol);
        const int nPhaseIdx = 1 - wPhaseIdx;

        mobility_[wPhaseIdx] =
            MaterialLaw::krw(materialParams, fluidState_.saturation(wPhaseIdx))
            / fluidState_.viscosity(wPhaseIdx);

        mobility_[nPhaseIdx] =
            MaterialLaw::krn(materialParams, fluidState_.saturation(wPhaseIdx))
            / fluidState_.viscosity(nPhaseIdx);

        // porosity
        porosity_ = problem.spatialParams().porosity(element, scv, elemSol);
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);

    }

    /*!
     * \brief Complete the fluid state
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The problem
     * \param element The element
     * \param scv The sub control volume
     * \param fluidState The fluid state
     *
     * Set temperature, saturations, capillary pressures, viscosities, densities and enthalpies.
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void completeFluidState(const ElemSol& elemSol,
                            const Problem& problem,
                            const Element& element,
                            const Scv& scv,
                            FluidState& fluidState)
    {
        Scalar t =  ParentType::temperature(elemSol, problem, element, scv);
        fluidState.setTemperature(t);

        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        const auto& priVars = ParentType::extractDofPriVars(elemSol, scv);

        const int wPhaseIdx = problem.spatialParams().template wettingPhase<FluidSystem>(element, scv, elemSol);
        if (formulation == TwoPFormulation::p0s1)
        {
            fluidState.setSaturation(phase1Idx, priVars[saturationIdx]);
            fluidState.setSaturation(phase0Idx, 1 - priVars[saturationIdx]);
        }
        else if (formulation == TwoPFormulation::p1s0)
        {
            fluidState.setSaturation(phase0Idx, priVars[saturationIdx]);
            fluidState.setSaturation(phase1Idx, 1 - priVars[saturationIdx]);
        }

        pc_ = MaterialLaw::pc(materialParams, fluidState.saturation(wPhaseIdx));
        if (formulation == TwoPFormulation::p0s1)
        {
            fluidState.setPressure(phase0Idx, priVars[pressureIdx]);
            fluidState.setPressure(phase1Idx, (wPhaseIdx == phase0Idx) ? priVars[pressureIdx] + pc_
                                                                       : priVars[pressureIdx] - pc_);
        }
        else if (formulation == TwoPFormulation::p1s0)
        {
            fluidState.setPressure(phase1Idx, priVars[pressureIdx]);
            fluidState.setPressure(phase0Idx, (wPhaseIdx == phase0Idx) ? priVars[pressureIdx] - pc_
                                                                       : priVars[pressureIdx] + pc_);
        }

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        for (int phaseIdx = 0; phaseIdx < ModelTraits::numPhases(); ++phaseIdx) {
            // compute and set the viscosity
            Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            fluidState.setViscosity(phaseIdx, mu);

            // compute and set the density
            Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, rho);

            // compute and set the enthalpy
            Scalar h = ParentType::enthalpy(fluidState, paramCache, phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }
    }

    /*!
     * \brief Returns the phase state for the control volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the saturation of a given phase within
     *        the control volume in \f$[-]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar saturation(int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume in \f$[kg/m^3]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume in \f$[kg/(m*s^2)=N/m^2=Pa]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns the capillary pressure within the control volume
     * in \f$[kg/(m*s^2)=N/m^2=Pa]\f$.
     */
    Scalar capillaryPressure() const
    { return pc_; }

    /*!
     * \brief Returns temperature inside the sub-control volume
     * in \f$[K]\f$.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the dynamic viscosity of the fluid within the
     *        control volume in \f$\mathrm{[Pa s]}\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar viscosity(int phaseIdx) const
    { return fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume in \f$[s*m/kg]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(int phaseIdx) const
    { return mobility_[phaseIdx]; }

    /*!
     * \brief Returns the average porosity within the control volume in \f$[-]\f$.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the permeability within the control volume in \f$[m^2]\f$.
     */
    const PermeabilityType& permeability() const
    { return permeability_; }

protected:
    FluidState fluidState_;

private:
    Scalar pc_;
    Scalar porosity_;
    PermeabilityType permeability_;
    Scalar mobility_[ModelTraits::numPhases()];
};

} // end namespace Dumux

#endif
