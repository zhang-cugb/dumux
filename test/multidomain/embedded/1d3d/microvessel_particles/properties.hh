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
 * \ingroup EmbeddedTests
 * \brief A test problem for the mixed-dimension embedded one-phase model
 */

#ifndef DUMUX_TEST_1D3D_MICROVESSEL_PROPERTIES_HH
#define DUMUX_TEST_1D3D_MICROVESSEL_PROPERTIES_HH

#include <type_traits>

#include <dune/grid/yaspgrid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/1p/model.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>

#include "spatialparams_tissue.hh"
#include "spatialparams_blood.hh"

#include "problem_tissue.hh"
#include "problem_blood.hh"

#ifndef BULKTYPETAG
#define BULKTYPETAG BulkCC
#endif
#ifndef LOWDIMTYPETAG
#define LOWDIMTYPETAG LowDimCC
#endif
#ifndef COUPLINGMODE
#define COUPLINGMODE EmbeddedCouplingMode::cylindersources
#endif
#ifndef BULKGRIDTYPE
#define BULKGRIDTYPE Dune::YaspGrid<3,Dune::EquidistantOffsetCoordinates<double,3>>
#endif

namespace Dumux::Properties {

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////// BULK ////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
// Create new type tags
namespace TTag {
struct Bulk { using InheritsFrom = std::tuple<OneP>; };
struct BulkCC { using InheritsFrom = std::tuple<Bulk, CCTpfaModel>; };
struct BulkBox { using InheritsFrom = std::tuple<Bulk, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Bulk> { using type = BULKGRIDTYPE; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Bulk> { using type = TissueProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Bulk>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<2, Scalar> >;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Bulk>
{
    using type = TissueSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                     GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag> struct EnableGridGeometryCache<TypeTag, TTag::Bulk> { static constexpr bool value = true; };
template<class TypeTag> struct EnableGridVolumeVariablesCache<TypeTag, TTag::Bulk> { static constexpr bool value = false; };
template<class TypeTag> struct EnableGridFluxVariablesCache<TypeTag, TTag::Bulk> { static constexpr bool value = false; };
template<class TypeTag> struct SolutionDependentAdvection<TypeTag, TTag::Bulk> { static constexpr bool value = false; };
template<class TypeTag> struct SolutionDependentMolecularDiffusion<TypeTag, TTag::Bulk> { static constexpr bool value = false; };
template<class TypeTag> struct SolutionDependentHeatConduction<TypeTag, TTag::Bulk> { static constexpr bool value = false; };

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////// EMBEDDED ////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
// Create new type tags
namespace TTag {
struct LowDim { using InheritsFrom = std::tuple<OneP>; };
struct LowDimCC { using InheritsFrom = std::tuple<LowDim, CCTpfaModel>; };
struct LowDimBox { using InheritsFrom = std::tuple<LowDim, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::LowDim> { using type = Dune::FoamGrid<1, 3>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::LowDim> { using type = BloodFlowProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::LowDim>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::LowDim>
{
    using type = BloodFlowSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                        GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag> struct EnableGridGeometryCache<TypeTag, TTag::LowDim> { static constexpr bool value = true; };
template<class TypeTag> struct EnableGridVolumeVariablesCache<TypeTag, TTag::LowDim> { static constexpr bool value = false; };
template<class TypeTag> struct EnableGridFluxVariablesCache<TypeTag, TTag::LowDim> { static constexpr bool value = false; };
template<class TypeTag> struct SolutionDependentAdvection<TypeTag, TTag::LowDim> { static constexpr bool value = false; };
template<class TypeTag> struct SolutionDependentMolecularDiffusion<TypeTag, TTag::LowDim> { static constexpr bool value = false; };
template<class TypeTag> struct SolutionDependentHeatConduction<TypeTag, TTag::LowDim> { static constexpr bool value = false; };

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////// COUPLING ////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

template<class Traits>
using TheCouplingManager = EmbeddedCouplingManager1d3d<Traits, COUPLINGMODE>;

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::BULKTYPETAG>
{ using type = TheCouplingManager<MultiDomainTraits<TypeTag, Properties::TTag::LOWDIMTYPETAG>>; };

template<class TypeTag>
struct PointSource<TypeTag, TTag::BULKTYPETAG>
{ using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSource<0>; };

template<class TypeTag>
struct PointSourceHelper<TypeTag, TTag::BULKTYPETAG>
{ using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSourceHelper<0>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::LOWDIMTYPETAG>
{ using type = TheCouplingManager<MultiDomainTraits<Properties::TTag::BULKTYPETAG, TypeTag>>; };

template<class TypeTag>
struct PointSource<TypeTag, TTag::LOWDIMTYPETAG>
{ using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSource<1>; };

template<class TypeTag>
struct PointSourceHelper<TypeTag, TTag::LOWDIMTYPETAG>
{ using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSourceHelper<1>; };

} // end namespace Dumux::Properties

#endif
