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
 *
 * \brief Particle tracker for the microvessel test
 */
#ifndef DUMUX_TEST_1D3D_MICROVESSEL_PARTICLES_HH
#define DUMUX_TEST_1D3D_MICROVESSEL_PARTICLES_HH

#include <algorithm>
#include <random>

#include <dumux/common/geometry/intersectspointgeometry.hh>

#include <dumux/particles/particle.hh>
#include <dumux/particles/simpleparticlecloud.hh>
#include <dumux/particles/particlevtkwriter.hh>

namespace Dumux {

template<class GG1D, class GG3D>
class MicrovesselParticleAlgorithm
{
    using Particle = Dumux::Particle<3>;
    using Cloud = SimpleParticleCloud<Particle>;
    using GlobalPosition = Particle::GlobalPosition;

    struct ParticleData
    {
        std::size_t domainIdx; // in which domain this particle is (1D=0, 3D=1)
        std::size_t elementIdx; // in which element this particle is
        double time; // at which time the particle lives
    };

    struct Inlet
    {
        std::size_t elementIdx;
        GlobalPosition position;
        double massFlowRate; // in particles per second
    };

    enum class FaceType1D { outflow, inflow, inner, converging, diverging, none };

    struct Line
    {
        static constexpr int mydimension = 1;
        static constexpr int coorddimension = GlobalPosition::dimension;
        const GlobalPosition& corner(int i) const { return c[i]; }
        std::array<std::reference_wrapper<const GlobalPosition>, 2> c;
    };

public:
    MicrovesselParticleAlgorithm(std::shared_ptr<const GG1D> gg1D,
                                 std::shared_ptr<const GG3D> gg3D,
                                 const std::vector<GlobalPosition>& velocities1D,
                                 const std::vector<GlobalPosition>& velocities3D,
                                 const std::vector<double>& extravasationProbability,
                                 const std::vector<double>& radius)
    : gridGeometry1D_(gg1D)
    , gridGeometry3D_(gg3D)
    , velocities1D_(velocities1D)
    , velocities3D_(velocities3D)
    , extravasationProbability_(extravasationProbability)
    , radius_(radius)
    , cloud_(std::make_shared<Cloud>())
    , rndGen_(std::random_device{}())
    {
        initializeData_();

        particleDensity_ = getParam<double>("Particles.Density"); // in particles / kg;
        velocityPertubationPercent_ = getParam<double>("Particles.VelocityPertubationPercent", 0.1);
        velPerturbation_ = std::uniform_real_distribution<double>{1.0-velocityPertubationPercent_, 1.0+velocityPertubationPercent_};
        posPerturbation_ = std::uniform_real_distribution<double>{-1.0, 1.0};
    }

    void run(const double timeStepSize)
    {
        // activate new particles at the inlet
        activateParticlesAtInlet_(timeStepSize);

        // move the particles with the given velocity field through the grid
        for (auto& p : particles(*cloud_))
            moveParticle_(p, timeStepSize);
    }

    std::shared_ptr<const Cloud> cloud() const
    { return cloud_; }

private:
    void activateParticlesAtInlet_(const double timeStepSize)
    {
        std::uniform_real_distribution<double> timeDist(0.0, timeStepSize);
        std::uniform_real_distribution<double> prob(0.0, 1.0);
        for (const auto& inlet : inlets_)
        {
            const auto averageNumParticles = inlet.massFlowRate*timeStepSize*particleDensity_;
            int numParticles = static_cast<int>(std::floor(averageNumParticles));
            const auto p = averageNumParticles - std::floor(averageNumParticles);
            const auto roll = prob(rndGen_);
            if (roll < p)
                ++numParticles;

            for (int i = 0; i < numParticles; ++i)
            {
                auto particle = cloud_->activateParticle();
                particle->setPosition(inlet.position);
                data_.resize(cloud_->size(false));
                auto& particleData = data_[particle->id()];
                particleData.domainIdx = 0;
                particleData.elementIdx = inlet.elementIdx;
                particleData.time = timeDist(rndGen_);
            }
        }
    }

    void moveParticle_(Particle& particle, const double timeStepSize)
    {
        auto& particleData = data_[particle.id()];
        particleData.time = timeStepSize;
        if (particleData.domainIdx == 0)
            moveParticleInVessel_(particle);
        else
            moveParticleInTissue_(particle);
    }

    void moveParticleInVessel_(Particle& particle)
    {
        auto& particleData = data_[particle.id()];
        const auto eIdx = particleData.elementIdx;
        const auto& element = gridGeometry1D_->element(eIdx);
        const auto geometry = element.geometry();

        auto velocity = velocities1D_[eIdx];
        velocity *= velPerturbation_(rndGen_);

        const auto startPos = projectPointOntoSegment_(particle.position(), geometry.corner(0), geometry.corner(1));
        auto endPos = startPos;
        endPos.axpy(particleData.time, velocity);

        // we stay inside this element
        if (intersectsPointGeometry(endPos, geometry))
        {
            // add some random pertubation
            auto rndUnitVector = GlobalPosition({posPerturbation_(rndGen_), posPerturbation_(rndGen_), posPerturbation_(rndGen_)});
            rndUnitVector /= rndUnitVector.two_norm();
            endPos.axpy(std::min((endPos-startPos).two_norm()*0.5, 0.7*radius_[eIdx]), rndUnitVector);
            particle.setPosition(endPos);
        }

        // we go to a neighboring element
        else
        {
            auto fvGeometry = localView(*gridGeometry1D_);
            fvGeometry.bindElement(element);

            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (scvf.unitOuterNormal()*(endPos-startPos) >= 0.0)
                {
                    if (scvf.boundary())
                        cloud_->deactivateParticle(particle);
                    else
                    {
                        const auto& newPos = scvf.ipGlobal();
                        particle.setPosition(newPos);
                        particleData.elementIdx = scvf.outsideScvIdx();
                        particleData.time -= (newPos-startPos).two_norm()/velocity.two_norm();
                        moveParticleInVessel_(particle);
                        break;
                    }
                }
            }
        }
    }

    void moveParticleInTissue_(Particle& particle)
    {
        // do nothing for now
    }

    void initializeData_()
    {
        // initially add all particles as inactive particles
        const auto numParticles = getParam<std::size_t>("Particles.TotalCloudSize");
        cloud_->resize(numParticles);
        data_.resize(numParticles);

        const auto& gridView = gridGeometry1D_->gridView();
        orientation_.resize(gridView.size(0));
        for (const auto& element : elements(gridView))
        {
            const auto eIdx = gridGeometry1D_->elementMapper().index(element);
            const auto geometry = element.geometry();
            orientation_[eIdx] = geometry.corner(1)-geometry.corner(0);
            orientation_[eIdx] /= orientation_[eIdx].two_norm();
        }

        // we assume a stationary velocity field, so we can precompute some stuff
        // find inlets and face types
        faceType_.resize(gridView.size(1), FaceType1D::none);
        for (const auto& element : elements(gridView))
        {
            const auto eIdx = gridGeometry1D_->elementMapper().index(element);
            const auto& velocity = velocities1D_[eIdx];

            auto fvGeometry = localView(*gridGeometry1D_);
            fvGeometry.bindElement(element);

            for (const auto& scvf : scvfs(fvGeometry))
            {
                const auto flux = velocity*scvf.unitOuterNormal();
                const auto indexInInside = (scvf.unitOuterNormal()*orientation_[eIdx] > 0.0) ? 1 : 0;
                const auto vIdx = gridView.indexSet().subIndex(element, indexInInside, 1);

                if (scvf.boundary())
                {
                    if (flux >= 0.0)
                        faceType_[vIdx] = FaceType1D::outflow;
                    else
                    {
                        faceType_[vIdx] = FaceType1D::inflow;
                        const auto radius = radius_[eIdx];
                        inlets_.emplace_back(Inlet{eIdx, scvf.ipGlobal(), velocity.two_norm()*scvf.area()*M_PI*radius*radius*1030.0});
                    }
                }

                // inner nodes only have one outside neighbor
                else if (scvf.numOutsideScvs() == 1)
                    faceType_[vIdx] = FaceType1D::inner;

                // bifurcation nodes have more
                else if (scvf.numOutsideScvs() > 1)
                {
                    int countInflow = 0; int countOutflow = 0;
                    if (flux >= 0.0)
                        countInflow += 1;
                    else
                        countOutflow += 1;

                    for (int i = 0; i < scvf.numOutsideScvs(); ++i)
                    {
                        const auto& flipScvf = fvGeometry.flipScvf(scvf.index(), i);
                        const auto nIdx = flipScvf.insideScvIdx();
                        const auto nFlux = velocities1D_[nIdx]*flipScvf.unitOuterNormal();
                        if (nFlux >= 0.0)
                            countInflow += 1;
                        else
                            countOutflow += 1;
                    }

                    if (countInflow > 1 && countOutflow == 1)
                        faceType_[vIdx] = FaceType1D::converging;
                    else if (countInflow > 0 && countOutflow == 2)
                        faceType_[vIdx] = FaceType1D::diverging;
                    else if (countOutflow > 2)
                        DUNE_THROW(Dune::NotImplemented, "Node has more than two outflow faces!");
                    else
                        DUNE_THROW(Dune::InvalidStateException, "Wrong faces configuration! "
                                    << "outflow: "  << countOutflow << ", inflow: " << countInflow << ".");
                }
            }
        }

    }

    GlobalPosition projectPointOntoSegment_(const GlobalPosition& p, const GlobalPosition& a, const GlobalPosition& b) const
    {
        const auto v = b - a;
        const auto w = p - a;

        const auto proj1 = v*w;
        if (proj1 <= 0.0)
            return a;

        const auto proj2 = v.two_norm2();
        if (proj2 <= proj1)
            return b;

        const auto t = proj1 / proj2;
        auto x = a;
        x.axpy(t, v);
        return x;
    }


    std::shared_ptr<const GG1D> gridGeometry1D_;
    std::shared_ptr<const GG3D> gridGeometry3D_;

    std::vector<GlobalPosition> velocities1D_, velocities3D_;
    std::vector<double> extravasationProbability_, radius_;
    std::vector<FaceType1D> faceType_;
    std::vector<GlobalPosition> orientation_;

    std::vector<std::size_t> outflowElements_;
    std::vector<std::array<std::size_t, 2>> offsetSize_;

    std::vector<Inlet> inlets_;
    std::vector<ParticleData> data_;
    std::shared_ptr<Cloud> cloud_;

    double particleDensity_;
    double velocityPertubationPercent_;
    double posPerturbationMagnitude_;

    std::mt19937 rndGen_;
    std::uniform_real_distribution<double> velPerturbation_;
    std::uniform_real_distribution<double> posPerturbation_;
};

} //end namespace Dumux

#endif
