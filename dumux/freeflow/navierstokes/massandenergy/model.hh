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
 * \ingroup NavierStokesModel
 *
 * \brief A single-phase, isothermal Navier-Stokes model
 *
 * This model implements a single-phase, isothermal Navier-Stokes model, solving the <B> momentum balance equation </B>
 * \f[
 \frac{\partial (\varrho \textbf{v})}{\partial t} + \nabla \cdot (\varrho \textbf{v} \textbf{v}^{\text{T}}) = \nabla \cdot (\mu (\nabla \textbf{v} + \nabla \textbf{v}^{\text{T}}))
   - \nabla p + \varrho \textbf{g} - \textbf{f}
 * \f]
 * By setting the runtime parameter <code>Problem.EnableInertiaTerms</code> to <code>false</code> the Stokes
 * equation can be solved. In this case the term
 * \f[
 *    \nabla \cdot (\varrho \textbf{v} \textbf{v}^{\text{T}})
 * \f]
 * is neglected.
 *
 * The <B> mass balance equation </B>
 * \f[
       \frac{\partial \varrho}{\partial t} + \nabla \cdot (\varrho \textbf{v}) - q = 0
 * \f]
 *
 * closes the system.
 *
 *
 * So far, only the staggered grid spatial discretization (for structured grids) is available.
 */

#ifndef DUMUX_NAVIERSTOKES_MASS_AND_ENERGY_MODEL_HH
#define DUMUX_NAVIERSTOKES_MASS_AND_ENERGY_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>
// #include <dumux/freeflow/properties.hh>
// #include <dumux/freeflow/nonisothermal/model.hh>
// #include <dumux/freeflow/nonisothermal/indices.hh>
// #include <dumux/freeflow/nonisothermal/iofields.hh>

#include "localresidual.hh"
#include "volumevariables.hh"
#include "fluxvariables.hh"
#include "indices.hh"
#include "./../iofields.hh"

#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/discretization/method.hh>
#include <dumux/flux/fourierslaw.hh>
#include <dumux/flux/fluxvariablescaching.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Traits for the Navier-Stokes model
 *
 * \tparam dimension The dimension of the problem
 */
template<int dimension>
struct NavierStokesMassAndEnergyModelTraits
{
    //! The dimension of the model
    static constexpr int dim() { return dimension; }

    //! There are as many momentum balance equations as dimensions
    //! and one mass balance equation.
    static constexpr int numEq() { return 1; }

    //! The number of phases is 1
    static constexpr int numFluidPhases() { return 1; }

    //! The number of components is 1
    static constexpr int numFluidComponents() { return 1; }

    //! Enable advection
    static constexpr bool enableAdvection() { return true; }

    //! The one-phase model has no molecular diffusion
    static constexpr bool enableMolecularDiffusion() { return false; }

    //! The model is isothermal
    static constexpr bool enableEnergyBalance() { return false; }

    //! The model does not include a turbulence model
    static constexpr bool usesTurbulenceModel() { return false; }

    //! return the type of turbulence model used
    static constexpr auto turbulenceModel()
    { return TurbulenceModel::none; }

    //! the indices
    using Indices = NavierStokesMassAndEnergyIndices;
};

/*!
 * \ingroup NavierStokesModel
 * \brief Traits class for the volume variables of the Navier-Stokes model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam FST The fluid state type
 * \tparam MT The model traits
 */
template<class PV,
         class FSY,
         class FST,
         class MT>
struct NavierStokesMassAndEnergyVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using FluidSystem = FSY;
    using FluidState = FST;
    using ModelTraits = MT;
};

// \{
///////////////////////////////////////////////////////////////////////////
// properties for the single-phase Navier-Stokes model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tag for the single-phase, isothermal Navier-Stokes model
struct NavierStokesMassAndEnergy { using InheritsFrom = std::tuple<ModelProperties>; };

}


//!< states some specifics of the Navier-Stokes model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::NavierStokesMassAndEnergy>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    static constexpr auto dim = GridView::dimension;
public:
    using type = NavierStokesMassAndEnergyModelTraits<dim>;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
template<class TypeTag>
struct FluidState<TypeTag, TTag::NavierStokesMassAndEnergy>{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = Dumux::ImmiscibleFluidState<Scalar, FluidSystem>;
};

//! The local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::NavierStokesMassAndEnergy> { using type = NavierStokesMassAndEnergyResidual<TypeTag>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::NavierStokesMassAndEnergy>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");
    static_assert(!FSY::isMiscible(), "The Navier-Stokes model only works with immiscible fluid systems.");

    using Traits = NavierStokesMassAndEnergyVolumeVariablesTraits<PV, FSY, FST, MT>;
public:
    using type = NavierStokesMassAndEnergyVolumeVariables<Traits>;
};

//! The flux variables
template<class TypeTag>
struct FluxVariables<TypeTag, TTag::NavierStokesMassAndEnergy> { using type = NavierStokesMassAndEnergyFluxVariables<TypeTag>; };

template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::NavierStokesMassAndEnergy> { using type = FluxVariablesCaching::EmptyCache<GetPropType<TypeTag, Properties::Scalar>>; };

template<class TypeTag>
struct FluxVariablesCacheFiller<TypeTag, TTag::NavierStokesMassAndEnergy> { using type = FluxVariablesCaching::EmptyCacheFiller; };

// ! The specific I/O fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::NavierStokesMassAndEnergy> { using type = NavierStokesIOFields; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::NavierStokesMassAndEnergy>
{
private:
    struct EmptyCouplingManager {};
public:
    using type = EmptyCouplingManager;
};

}
// }

} // end namespace

#endif // DUMUX_NAVIERSTOKES_MODEL_HH
