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
* \ingroup OnePNCMinModel
* \brief A single-phase, multi-component model considering mineralization processes.
*
* This model implements one-phase n-component flow of a compressible fluid composed of
* the n components \f$\kappa \f$ in combination with mineral precipitation and dissolution
* of the solid phases. The standard multiphase Darcy
* approach is used as the equation for the conservation of momentum:
* \f[
v = - \frac{k_{r}}{\mu} \mbox{\bf K}
\left(\text{grad}\, p - \varrho_{f} \mbox{\bf g} \right)
* \f]
*
* By inserting this into the equations for the conservation of the
* components, one gets one transport equation for each component
* \f[
 \frac{\partial ( \varrho_f X^\kappa \phi  )}
{\partial t} -  \text{div} \left\{ \varrho_f X^\kappa
\frac{k_{r}}{\mu} \mbox{\bf K}
(\text{grad}\, p - \varrho_{f}  \mbox{\bf g}) \right\}
- \text{div} \left\{{\bf D_{pm}^\kappa} \varrho_{f} \text{grad}\, X^\kappa \right\}
-  q_\kappa = 0 \qquad \kappa \in \{w, a,\cdots \}
* \f]
*
* The solid or mineral phases are assumed to consist of a single component.
* Their mass balance consist only of a storage and a source term:
* \f[
 \frac{\partial \varrho_\lambda \phi_\lambda )} {\partial t} = q_\lambda
* \f]
*
* All equations are discretized using a vertex-centered finite volume (box)
* or cell-centered finite volume scheme as spatial and the implicit Euler method as time
* discretization.
*
* The primary variables are the pressure \f$p\f$ and the mole fractions of the
* dissolved components \f$x^k\f$. The primary variable of the solid phases is the volume
* fraction
\f$\phi_\lambda = \frac{V_\lambda}{V_{total}}\f$.
*
* The source an sink terms link the mass balances of the n-transported component to the
* solid phases. The porosity \f$\phi\f$ is updated according to the reduction of the initial
* (or solid-phase-free porous medium) porosity \f$\phi_0\f$ by the accumulated volume
* fractions of the solid phases:
* \f$ \phi = \phi_0 - \sum (\phi_\lambda)\f$
* Additionally, the permeability is updated depending on the current porosity.
*/

#ifndef DUMUX_1PNCMIN_MODEL_HH
#define DUMUX_1PNCMIN_MODEL_HH

#include <dumux/porousmediumflow/1pnc/model.hh>
#include <dumux/porousmediumflow/mineralization/model.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>

namespace Dumux
{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
NEW_TYPE_TAG(OnePNCMin, INHERITS_FROM(OnePNC, Mineralization));
NEW_TYPE_TAG(OnePNCMinNI, INHERITS_FROM(OnePNCMin, NonIsothermal));

//////////////////////////////////////////////////////////////////
// Property tags for the isothermal 2pncmin model
//////////////////////////////////////////////////////////////////
SET_TYPE_PROP(OnePNCMin, NonMineralizationVolumeVariables, OnePNCVolumeVariables<TypeTag>);     //!< the VolumeVariables property

//! Set the vtk output fields specific to this model
SET_PROP(OnePNCMin, NonMineralizationVtkOutputFields)
{
private:
   using FluidSystem =  typename GET_PROP_TYPE(TypeTag, FluidSystem);
   static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
public:
    using type = OnePNCVtkOutputFields<FluidSystem, phaseIdx>;
};

//////////////////////////////////////////////////////////////////
// Properties for the non-isothermal 2pncmin model
//////////////////////////////////////////////////////////////////
SET_TYPE_PROP(OnePNCMinNI, IsothermalVolumeVariables, MineralizationVolumeVariables<TypeTag>);  //!< set isothermal VolumeVariables

//! isothermal vtkoutput
SET_PROP(OnePNCMinNI, IsothermalVtkOutputFields)
{
private:
   using NonMineralizationVtkOutputFields =  typename GET_PROP_TYPE(TypeTag, NonMineralizationVtkOutputFields);
   using FluidSystem =  typename GET_PROP_TYPE(TypeTag, FluidSystem);

public:
    using type = MineralizationVtkOutputFields<NonMineralizationVtkOutputFields, FluidSystem>;
};

SET_TYPE_PROP(OnePNCMinNI, IsothermalLocalResidual, MineralizationLocalResidual<TypeTag>);      //!< set isothermal output fields

SET_TYPE_PROP(OnePNCMinNI,
              ThermalConductivityModel,
              ThermalConductivityAverage<typename GET_PROP_TYPE(TypeTag, Scalar)>); //!< Use the average for effective conductivities

SET_PROP(OnePNCMinNI, IsothermalNumEq) {
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem));
public:
    static const int value = FluidSystem::numComponents + FluidSystem::numSPhases;
};

//! use 1pnc indices for the isothermal indices
SET_PROP(OnePNCMinNI, IsothermalIndices)
{
private:
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
public:
    using type = OnePNCIndices<phaseIdx>;
};

} // end namespace Properties
} // end namespace Dumux

#endif
