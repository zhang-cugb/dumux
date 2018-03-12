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
 * \ingroup OnePModel
 * \brief A single-phase, isothermal flow model using the fully implicit scheme.
 *
 * Single-phase, isothermal flow model, which uses a standard Darcy approach as the
 * equation for the conservation of momentum:
 * \f[
 v = - \frac{\textbf K}{\mu}
 \left(\textbf{grad}\, p - \varrho {\textbf g} \right)
 * \f]
 *
 * and solves the mass continuity equation:
 * \f[
 \phi \frac{\partial \varrho}{\partial t} + \text{div} \left\lbrace
 - \varrho \frac{\textbf K}{\mu} \left( \textbf{grad}\, p -\varrho {\textbf g} \right) \right\rbrace = q,
 * \f]
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 * The model supports compressible as well as incompressible fluids.
 */

#ifndef DUMUX_1P_MODEL_HH
#define DUMUX_1P_MODEL_HH

#include <dumux/common/properties.hh>

#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>
#include <dumux/material/fluidstates/immiscible.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "vtkoutputfields.hh"

namespace Dumux {
/*!
 * \ingroup OnePModel
 * \brief Specifies a number properties of single-phase models.
 */
struct OnePModelTraits
{
    using Indices = OnePIndices;

    static constexpr int numEq() { return 1; }
    static constexpr int numPhases() { return 1; }
    static constexpr int numComponents() { return 1; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return false; }
    static constexpr bool enableEnergyBalance() { return false; }
};

namespace Properties {
//! The type tags for the isothermal single phase model
NEW_TYPE_TAG(OneP, INHERITS_FROM(PorousMediumFlow));
//! The type tags for the non-isothermal single phase model
NEW_TYPE_TAG(OnePNI, INHERITS_FROM(OneP, NonIsothermal));

///////////////////////////////////////////////////////////////////////////
// properties for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
SET_TYPE_PROP(OneP, VtkOutputFields, OnePVtkOutputFields);            //!< default vtk output fields specific to this model
SET_TYPE_PROP(OneP, LocalResidual, ImmiscibleLocalResidual<TypeTag>); //!< the local residual function
SET_TYPE_PROP(OneP, VolumeVariables, OnePVolumeVariables<TypeTag>);   //!< the VolumeVariables property
SET_TYPE_PROP(OneP, ModelTraits, OnePModelTraits);                    //!< states some specifics of the one-phase model
SET_TYPE_PROP(OneP, Indices, OnePIndices);                            //!< Indices used by the one-phase model

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(OneP, FluidState)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = ImmiscibleFluidState<Scalar, FluidSystem>;
};

///////////////////////////////////////////////////////////////////////////
// properties for the non-isothermal single phase model
///////////////////////////////////////////////////////////////////////////
SET_TYPE_PROP(OnePNI, IsothermalVtkOutputFields, OnePVtkOutputFields);              //!< the isothermal vtk output fields
SET_TYPE_PROP(OnePNI, IsothermalVolumeVariables, OnePVolumeVariables<TypeTag>);     //!< Vol vars of the isothermal model
SET_TYPE_PROP(OnePNI, IsothermalLocalResidual, ImmiscibleLocalResidual<TypeTag>);   //!< Local residual of the isothermal model
SET_TYPE_PROP(OnePNI, IsothermalIndices, OnePIndices);                              //!< Indices of the isothermal model
SET_TYPE_PROP(OnePNI, IsothermalModelTraits, OnePModelTraits);                      //!< The model traits of the isothermal model
SET_TYPE_PROP(OnePNI,
              ThermalConductivityModel,
              ThermalConductivityAverage<typename GET_PROP_TYPE(TypeTag, Scalar)>); //!< Use the average for effective conductivities

} // end namespace Properties
} // end namespace Dumux

#endif
