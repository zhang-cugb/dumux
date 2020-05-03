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
 * \brief Particle iterators
 */
#ifndef DUMUX_PARTICLES_PARTICLE_ITERATOR_HH
#define DUMUX_PARTICLES_PARTICLE_ITERATOR_HH

#include <type_traits>
#include <dune/common/iteratorfacades.hh>

namespace Dumux {

/*!
 * \ingroup Particles
 * \brief Iterator over active particles
 */
template<class Container, class Particle>
class ActiveParticleIterator
: public Dune::ForwardIteratorFacade<ActiveParticleIterator<Container, Particle>, Particle>
{
    using Iterator = std::conditional_t<std::is_const_v<Particle>,
                                        typename Container::const_iterator,
                                        typename Container::iterator>;
public:
    //! default construtor
    ActiveParticleIterator() = default;

    //! create from iterator
    ActiveParticleIterator(const Iterator& it, const Iterator& endIt)
    : it_{it}, endIt_{endIt}
    {
        while (it_ != endIt_ && !it_->isActive())
            ++it_;
    }

    //! dereferencing yields a particle
    Particle& dereference() const
    {
        return *it_;
    }

    //! test for equality
    bool equals(const ActiveParticleIterator& other) const
    {
        return it_ == other.it_;
    }

    //! increment until the next active particle
    void increment()
    {
        ++it_;
        while (it_ != endIt_ && !it_->isActive())
            ++it_;
    }

private:
    Iterator it_;
    const Iterator endIt_;
};

} // end namespace Dumux

#endif
