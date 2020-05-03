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
 * \ingroup Particles
 * \brief A simple particle
 */
#ifndef DUMUX_PARTICLES_PARTICLE_HH
#define DUMUX_PARTICLES_PARTICLE_HH

#include <dune/common/fvector.hh>

namespace Dumux {

/*!
 * \ingroup Particles
 * \brief a basic particle
 */
template<int dimWorld, class ctype = double>
class Particle
{
public:
    static constexpr int coordDimension = dimWorld;
    using CoordScalar = ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

    Particle(std::size_t id, const GlobalPosition& p, bool active = true)
    : id_(id), position_(p), active_(active) {}

    Particle(std::size_t id, bool active = true)
    : Particle(id, GlobalPosition(0.0), active) {}

    const GlobalPosition& position() const
    { return position_; }

    //! set a new global position
    void setPosition(const GlobalPosition& pos)
    { position_ = pos; }

    //! move by length in direction
    void move(const CoordScalar length, const GlobalPosition& direction)
    { position_.axpy(length, direction); }

    //! the particle identifier
    std::size_t id() const
    { return id_; }

    //! if the particle is active
    bool isActive() const
    { return active_; }

    //! deactivate the particle
    void deactivate()
    { active_ = false; }

    //! activate the particle
    void activate()
    { active_ = true; }

private:
    std::size_t id_;
    GlobalPosition position_;
    bool active_;
};

} // end namespace Dumux

#endif
