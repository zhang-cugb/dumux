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

#ifndef DUMUX_ONEP_ROTATION_SYMMETRY_PROPERTIES_HH
#define DUMUX_ONEP_ROTATION_SYMMETRY_PROPERTIES_HH

// ## The properties (`properties.hh`)
// This file defines the `TypeTag` used for the single-phase rotation symmetry simulation, for
// which we then define the necessary properties.
// [[content]]
// ### Includes
// [[details]] includes
#include <dune/grid/yaspgrid.hh> // for `Dune::YaspGrid`
#include <dumux/discretization/box.hh> // for `TTag::BoxModel`

// The local residual for incompressible flow is included.
// The one-phase flow model (included above) uses a default implementation of the
// local residual for single-phase flow. However, in this example we are using an
// incompressible fluid phase. Therefore, we are including the specialized local
// residual which contains functionality to analytically compute the entries of
// the Jacobian matrix. We will use this in the main file.
#include <dumux/porousmediumflow/1p/model.hh> // for `TTag::OneP`
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

// For rotational symmetric problems we use special geometry traits
#include <dumux/discretization/rotationsymmetricgridgeometrytraits.hh>

#include "problem.hh"
#include "spatialparams.hh"
// [[/details]]

namespace Dumux::Properties {

// A `TypeTag` for our simulation is created which inherits from the one-phase flow model
// and the cell centered finite volume scheme with two-point-flux discretization scheme:
namespace TTag {
struct OnePRotSym { using InheritsFrom = std::tuple<OneP, BoxModel>; };
}

// We use a structured 1D grid with an offset:
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePRotSym>
{ using type =  Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>; };

// Special grid geometry traits are needed
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::OnePRotSym>
{
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using GGTraits = RotationSymmetricGridGeometryTraits<BoxDefaultGridGeometryTraits<GridView>, RotationPolicy::disc>;

    using type = BoxFVGridGeometry<Scalar, GridView, enableCache, GGTraits>;
};

// We use the local residual that contains analytic derivative methods for incompressible flow:
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::OnePRotSym>
{ using type = OnePIncompressibleLocalResidual<TypeTag>; };

// The problem class specifies initial and boundary conditions:
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePRotSym>
{ using type = RotSymExampleProblem<TypeTag>; };

// We define the spatial parameters for our simulation:
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePRotSym>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = RotSymExampleSpatialParams<GridGeometry, Scalar>;
};

// In the following we define the fluid system to be used:
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePRotSym>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

} // end namespace Dumux::Properties
// [[/content]]
#endif
