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
 * \ingroup MortarCoupling.
 * \brief Defines the interface for projectors used in the
 *        context of coupling different domains via mortar spaces.
 */
#ifndef DUMUX_MORTAR_PROJECTOR_INTERFACE_HH
#define DUMUX_MORTAR_PROJECTOR_INTERFACE_HH

namespace Dumux {

/*!
 * \ingroup MortarCoupling.
 * \brief Defines the interface for projectors used in the
 *        context of coupling different domains via mortar spaces.
 */
template< class SubDomainSolutionVector, class MortarSolutionVector >
struct MortarProjectorInterface
{
    //! Projects solution from the mortar space to sub-domain space
    virtual SubDomainSolutionVector projectMortarToSubDomain(const MortarSolutionVector& x) const = 0;

    //! Projects solution from the sub-domain space to mortar space
    virtual MortarSolutionVector projectSolutionToMortar(const SubDomainSolutionVector& x) const = 0;

    //! Projects sub-domain fluxes from the sub-domain space to mortar space
    virtual MortarSolutionVector projectFluxToMortar(const SubDomainSolutionVector& x) const = 0;
};

} // end namespace Dumux

#endif
