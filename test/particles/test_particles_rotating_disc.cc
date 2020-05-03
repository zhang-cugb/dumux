#include <config.h>

#include <memory>
#include <cmath>
#include <random>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/math.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/particles/particle.hh>
#include <dumux/particles/simpleparticlecloud.hh>
#include <dumux/particles/particlevtkwriter.hh>

//! initialize particles on a disc
template<class Cloud>
void initialize(Cloud& cloud)
{
    const double radius = 4.0;
    const std::size_t numParticles = 100;
    cloud.resize(numParticles);

    std::mt19937 gen(0.0);
    std::uniform_real_distribution<double> distRadiusSq(0.0, radius*radius);
    std::uniform_real_distribution<double> distAngle(0.0, 2.0*M_PI);
    for (std::size_t pIdx = 0; pIdx < numParticles; ++pIdx)
    {
        auto particleIt = cloud.activateParticle();
        const auto radius = std::sqrt(distRadiusSq(gen));
        const auto angle = distAngle(gen);
        particleIt->setPosition({radius*std::cos(angle), radius*std::sin(angle)});
    }
}

//! evolve particles on rotating disc
template<class Cloud>
void evolve(Cloud& cloud, const double dt)
{
    const double angularVelocity = 2.0*M_PI/5.0;
    std::mt19937 gen(0.0);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (auto& particle : particles(cloud))
    {
        const auto& pos = particle.position();
        const auto radius = std::hypot(pos[0], pos[1]) + 0.1*dist(gen);
        const auto angle = std::atan2(pos[1], pos[0]) + dt*angularVelocity;
        particle.setPosition({radius*std::cos(angle), radius*std::sin(angle)});

        // deactivate the particle with a probability of 1%
        if (dist(gen) < 0.05)
            cloud.deactivateParticle(particle);
    }
}

int main(int argc, char** argv) try
{
    using namespace Dumux;

    Dune::MPIHelper::instance(argc, argv);

    static constexpr int dimWorld = 2;
    auto cloud = std::make_shared<SimpleParticleCloud<Particle<dimWorld>>>();
    initialize(*cloud);

    // vtk output
    ParticlePVDWriter writer(cloud, "test_particles_rotating_disc");
    std::vector<double> radius = linspace(0.5, 0.7, cloud->size());
    writer.addParticleData(radius, "radius");
    writer.write(0.0);

    TimeLoop<double> timeLoop(0.0, 0.1, 5.0);
    timeLoop.start();
    while (!timeLoop.finished())
    {
        evolve(*cloud, timeLoop.timeStepSize());
        timeLoop.advanceTimeStep();
        writer.write(timeLoop.time());
        timeLoop.reportTimeStep();

        std::cout << "-- number of active particles: " << cloud->size() << "\n"
                  << "-- number of inactive particles: " << cloud->size(false)-cloud->size(true) << std::endl;
    }

    timeLoop.finalize();

    return 0;
} // end main
catch (const Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 1;
}
