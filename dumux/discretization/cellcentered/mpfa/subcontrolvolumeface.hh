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
 * \brief Class for the sub-control volume face in mpfa schemes
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_SUBCONTROLVOLUMEFACE_HH

#include "methods.hh"

namespace Dumux
{
//! Forward declaration of the method specific implementations
//! Available implementations have to be included at the end of this file.
template<MpfaMethods M, class G, class GT, typename I>
class CCMpfaSubControlVolumeFaceImplementation;

/*!
 * \ingroup Mpfa
 * \brief Class for a sub control volume face in mpfa methods, i.e a part of the boundary
 *        of a control volume we compute fluxes on. This class inherits from the actual implementations
 *        and defines the constructor interface.
 *
 * \param M the mpfa method used
 * \param G the geometry type for the scvf geometries
 * \param GT the traits class for the geometry type
 * \param I the type used for indices
 */
template<MpfaMethods M, class G, class GT, typename I>
class CCMpfaSubControlVolumeFace : public CCMpfaSubControlVolumeFaceImplementation<M, G, GT, I>
{
    using ParentType = CCMpfaSubControlVolumeFaceImplementation<M, G, GT, I>;
    using Geometry = G;
    using IndexType = I;

    using Scalar = typename Geometry::ctype;
    static const int dim = Geometry::mydimension;
    static const int dimWorld = Geometry::coorddimension;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    //! state the traits class puclicly
    using Traits = GT;

    //! For convenience, state the corner storage vector from the traits
    using CornerVector = typename Traits::template CornerStorage<dim, dimWorld>::Type;

    /*!
     * \brief Constructor
     *
     * \param geomHelper The mpfa geometry helper
     * \param corners The corners of the scv face
     * \param unitOuterNormal The unit outer normal vector of the scvf
     * \param vIdxGlobal The global vertex index the scvf is connected to
     * \param localIndex Some element local index (e.g. the local vertex index in mpfa-o-fps scheme)
     * \param scvfIndex The global index of this scv face
     * \param insideScvIdx The inside scv index connected to this face
     * \param outsideScvIndices The outside scv indices connected to this face
     * \param q The parameterization of the quadrature point on the scvf for flux calculation
     * \param boundary Boolean to specify whether or not the scvf is on a boundary
     */
    template<class MpfaHelper>
    CCMpfaSubControlVolumeFace(const MpfaHelper& helper,
                               CornerVector&& corners,
                               GlobalPosition&& unitOuterNormal,
                               IndexType vIdxGlobal,
                               unsigned int localIndex,
                               IndexType scvfIndex,
                               IndexType insideScvIdx,
                               const std::vector<IndexType>& outsideScvIndices,
                               Scalar q,
                               bool boundary)
    : ParentType(helper,
                 std::forward<CornerVector>(corners),
                 std::forward<GlobalPosition>(unitOuterNormal),
                 vIdxGlobal,
                 localIndex,
                 scvfIndex,
                 insideScvIdx,
                 outsideScvIndices,
                 q,
                 boundary)
    {}
};
} // end namespace

//! The available implementations should be included here
// #include <dumux/discretization/cellcentered/mpfa/lmethod/subcontrolvolumeface.hh>
#include <dumux/discretization/cellcentered/mpfa/omethod/subcontrolvolumeface.hh>
// #include <dumux/discretization/cellcentered/mpfa/omethodfps/subcontrolvolumeface.hh>

#endif
