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
 * \ingroup Concepts
 * \brief Contains concepts that can be used together with the
 *        Dune infracstructure for concepts
 * \note See https://gitlab.dune-project.org/core/dune-common/-/blob/master/dune/common/concept.hh
 */
#ifndef DUMUX_CONCEPTS_HH
#define DUMUX_CONCEPTS_HH

namespace Dumux::Concept {

struct Resizable
{
    template<class V>
    auto require(const V& v) -> decltype(std::declval<V>().resize(0))
    {}
};

struct MatrixResizable
{
    template<class M>
    auto require(const M& m) -> decltype(std::declval<M>().resize(0, 0))
    {}
};

} // end namespace Dumux::Concept

#endif
