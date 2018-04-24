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
 * \ingroup FreeflowNCModel
 *
 * \copydoc Dumux::NavierStokesModel
 *
 * The system is closed by a <B> component mass/mole balance equation </B> for each component \f$\kappa\f$:
 * \f[
 *    \frac{\partial \left(\varrho X^\kappa\right)}{\partial t}
 *    + \nabla \cdot \left( \varrho {\boldsymbol{v}} X^\kappa
 *    - D^\kappa \varrho \frac{M^\kappa}{M} \textbf{grad}\, x^\kappa \right)
 *    - q^\kappa = 0
 * \f]
 *
 * Alternatively, one component balance equation can be replace by a <B> total mass/mole balance equation </B>:
 *
 * \f[
 *    \frac{\partial \varrho_g}{\partial t}
 *    + \nabla \cdot \left(
 *        \varrho {\boldsymbol{v}}
 *        - \sum_\kappa D^\kappa \varrho \frac{M^\kappa}{M} \textbf{grad}\, x^\kappa
 *      \right)
 *    - q = 0
 * \f]
 *
 * So far, only the staggered grid spatial discretization (for structured grids) is available.
 */

#ifndef DUMUX_FREEFLOW_NC_MODEL_HH
#define DUMUX_FREEFLOW_NC_MODEL_HH

#include <dumux/common/properties.hh>

#include <dumux/freeflow/navierstokes/model.hh>
#include <dumux/freeflow/nonisothermal/model.hh>
#include <dumux/freeflow/nonisothermal/indices.hh>
#include <dumux/freeflow/nonisothermal/navierstokesnivtkoutputfields.hh>
#include <dumux/discretization/fickslaw.hh>
#include <dumux/discretization/fourierslaw.hh>

#include "volumevariables.hh"
#include "indices.hh"
#include "localresidual.hh"
#include "fluxvariables.hh"
#include "vtkoutputfields.hh"

#include <dumux/assembly/staggeredlocalresidual.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/material/fluidstates/compositional.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowNCModel
 * \brief Traits for the multi-component free-flow model
 */
template<int dimension, int nComp, int phaseIdx, int replaceCompEqIdx, bool useM>
struct NavierStokesNCModelTraits
{
    //! The dimension of the model
    static constexpr int dim() { return dimension; }

    //! There are as many momentum balance equations as dimensions
    //! and as many balance equations as components.
    static constexpr int numEq() { return dimension+nComp; }

    //! The number of phases is always 1
    static constexpr int numPhases() { return 1; }

    //! The number of components
    static constexpr int numComponents() { return nComp; }

    //! Use moles or not
    static constexpr bool useMoles() { return useM; }

    //! Enable advection
    static constexpr bool enableAdvection() { return true; }

    //! The one-phase model has no molecular diffusion
    static constexpr bool enableMolecularDiffusion() { return true; }

    //! The model is isothermal
    static constexpr bool enableEnergyBalance() { return false; }

    //! the indices
    using Indices = FreeflowNCIndices<dim(), numEq(), phaseIdx, replaceCompEqIdx>;
};

///////////////////////////////////////////////////////////////////////////
// properties for the single-phase, multi-component free-flow model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, multi-component isothermal free-flow model
NEW_TYPE_TAG(NavierStokesNC, INHERITS_FROM(FreeFlow));

//! The type tag for the single-phase, multi-component non-isothermal free-flow model
NEW_TYPE_TAG(NavierStokesNCNI, INHERITS_FROM(NavierStokesNC));

///////////////////////////////////////////////////////////////////////////
// default property values
///////////////////////////////////////////////////////////////////////////

//!< states some specifics of the free-flow model
SET_PROP(NavierStokesNC, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr int dim = GridView::dimension;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
    static constexpr int replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
public:
    using type = NavierStokesNCModelTraits<dim, numComponents, phaseIdx, replaceCompEqIdx, useMoles>;
};

SET_INT_PROP(NavierStokesNC, PhaseIdx, 0); //!< Defines the phaseIdx
SET_BOOL_PROP(NavierStokesNC, UseMoles, false); //!< Defines whether molar (true) or mass (false) density is used
SET_INT_PROP(NavierStokesNC, ReplaceCompEqIdx, 0); //<! Set the ReplaceCompEqIdx to 0 by default
SET_BOOL_PROP(NavierStokesNC, EnableInertiaTerms, true); //!< Consider inertia terms by default
SET_BOOL_PROP(NavierStokesNC, NormalizePressure, true); //!< Normalize the pressure term in the momentum balance by default


//! The local residual
SET_TYPE_PROP(NavierStokesNC, LocalResidual, FreeflowNCResidual<TypeTag>);

//! Set the volume variables property
SET_PROP(NavierStokesNC, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
public:
    using type = FreeflowNCVolumeVariables<Traits>;
};


//! The flux variables
SET_TYPE_PROP(NavierStokesNC, FluxVariables, FreeflowNCFluxVariables<TypeTag>);

//! The flux variables cache class, by default the one for free flow
SET_TYPE_PROP(NavierStokesNC, FluxVariablesCache, FreeFlowFluxVariablesCache<TypeTag>);

//! The specific vtk output fields
SET_PROP(NavierStokesNC, VtkOutputFields)
{
private:
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
public:
     using type = FreeflowNCVtkOutputFields<FVGridGeometry, FluidSystem, phaseIdx>;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(NavierStokesNC, FluidState)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = CompositionalFluidState<Scalar, FluidSystem>;
};

//! Use Fick's law for molecular diffusion per default
SET_TYPE_PROP(NavierStokesNC, MolecularDiffusionType, FicksLaw<TypeTag>);

//////////////////////////////////////////////////////////////////////////
// Property values for non-isothermal multi-component free-flow model
//////////////////////////////////////////////////////////////////////////

//! The model traits of the non-isothermal model
SET_PROP(NavierStokesNCNI, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr int dim = GridView::dimension;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
    static constexpr int replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    using IsothermalModelTraits = NavierStokesNCModelTraits<dim, numComponents, phaseIdx, replaceCompEqIdx, useMoles>;
public:
    using type = FreeflowNIModelTraits<IsothermalModelTraits>;
};

//! The non-isothermal vtk output fields
SET_PROP(NavierStokesNCNI, VtkOutputFields)
{
private:
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
    using IsothermalFields = FreeflowNCVtkOutputFields<FVGridGeometry, FluidSystem, phaseIdx>;
public:
     using type = NavierStokesNonIsothermalVtkOutputFields<IsothermalFields>;
};

//! Use Fourier's Law as default heat conduction type
SET_TYPE_PROP(NavierStokesNCNI, HeatConductionType, FouriersLaw<TypeTag>);

// \}
} // end namespace Properties
} // end namespace Dumux


#endif
