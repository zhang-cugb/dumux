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
 * \brief A simple implementation of a cloud of particles
 */
#ifndef DUMUX_PARTICLES_SIMPLE_PARTICLE_CLOUD_HH
#define DUMUX_PARTICLES_SIMPLE_PARTICLE_CLOUD_HH

#include <vector>

#include <dune/common/iteratorrange.hh>

#include <dumux/particles/particle.hh>
#include <dumux/particles/particleiterator.hh>

namespace Dumux {

template<class P>
class SimpleParticleCloud
{
public:
    using Particle = P;

    SimpleParticleCloud() = default;

    /*!
     * \brief Resize the cloud to contain numParticles particles
     */
    void resize(std::size_t numParticles)
    {
        // pointers are invalidated when resizing
        inactiveParticles_.clear();
        if (numParticles > particles_.size())
        {
            particles_.reserve(numParticles);
            const auto newParticles = numParticles - particles_.size();
            for (std::size_t i = 0; i < newParticles; ++i)
                particles_.emplace_back(particles_.size(), false);
            if (particles_.size() != numParticles)
                DUNE_THROW(Dune::Exception, "");
        }
        else
        {
            particles_.erase(particles_.begin()+numParticles, particles_.end());
            particles_.shrink_to_fit();
            if (particles_.size() != numParticles)
                DUNE_THROW(Dune::Exception, "");
        }

        // reinitialize deactived particle list
        const auto pEndIt = particles_.end();
        for (auto pIt = particles_.begin(); pIt != pEndIt; ++pIt)
            if (!pIt->isActive())
                inactiveParticles_.push_back(pIt);
    }

    /*!
     * \brief the number of particles in the cloud
     */
    std::size_t size (bool onlyActive = true) const
    {
        if (onlyActive)
            return particles_.size()-inactiveParticles_.size();
        else
            return particles_.size();
    }

    /*!
     * \brief Bring one dorming particle to life
     * Activates one particle from the pool of inactive particles
     * \note the particle is automatically activated
     * \note Take care to place it somewhere afterwards! Otherwise it's hanging somewhere in space.
     * \note If no particles are left new particles are added (this invalidates pointers and may be costly)
     */
    typename std::vector<Particle>::iterator activateParticle()
    {
        if (inactiveParticles_.empty())
        {
            std::cout << "-- running out of particles -- adding 10 new particles." << std::endl;
            this->resize(particles_.size() + 10);
        }

        typename std::vector<Particle>::iterator pIt(inactiveParticles_.back());
        pIt->activate();
        inactiveParticles_.pop_back();
        return pIt;
    }

    /*!
     * \brief Deactivate a given particle
     */
    void deactivateParticle(const Particle& p)
    {
        particles_[p.id()].deactivate();
        inactiveParticles_.push_back(std::next(particles_.begin(), p.id()));
    }

    /*!
     * \brief Const range generator for iterating over all active particles
     * \note usage: for(const auto& particle : particles(cloud))
     * This is a free function found by means of ADL
     */
    friend inline Dune::IteratorRange<ActiveParticleIterator<std::vector<Particle>, const Particle>>
    particles(const SimpleParticleCloud& cloud)
    { return { {std::cbegin(cloud.particles_), std::cend(cloud.particles_)}, {std::cend(cloud.particles_), std::cend(cloud.particles_)} }; }

    /*!
     * \brief Non-const range generator for iterating over all active particles
     * \note usage: for(const auto& particle : particles(cloud))
     * This is a free function found by means of ADL
     */
    friend inline Dune::IteratorRange<ActiveParticleIterator<std::vector<Particle>, Particle>>
    particles(SimpleParticleCloud& cloud)
    { return { {std::begin(cloud.particles_), std::end(cloud.particles_)}, {std::end(cloud.particles_), std::end(cloud.particles_)} }; }

private:
    std::vector<Particle> particles_;
    std::vector<typename std::vector<Particle>::iterator> inactiveParticles_;
};

} // end namespace Dumux

#endif
