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

#ifndef DUMUX_ONEP_ROTATION_SYMMETRY_SPATIAL_PARAMS_HH
#define DUMUX_ONEP_ROTATION_SYMMETRY_SPATIAL_PARAMS_HH
// [[content]]
#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {

template<class GridGeometry, class Scalar>
class RotSymExampleSpatialParams
: public FVSpatialParamsOneP<GridGeometry, Scalar, RotSymExampleSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = RotSymExampleSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
public:
    using PermeabilityType = Scalar;
    RotSymExampleSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    { permeability_ = getParam<Scalar>("SpatialParams.Permeability"); }

    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }
private:
    Scalar permeability_;
};

} // end namespace Dumux
// [[/content]]
#endif
