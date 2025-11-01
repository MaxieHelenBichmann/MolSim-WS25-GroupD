#ifndef CUBOID_GENERATOR_H
#define CUBOID_GENERATOR_H

#include <cstddef>

#include "particles/Generator.h"
#include "utils/Vector.h"
namespace mol_sim {
using R3 = Vector<double, 3>;
using N3 = Vector<size_t, 3>;
class CuboidGenerator : public Generator {
   private:
    R3 position;
    R3 velocity;
    N3 num_particles;
    double mass;
    double distance;
    double avg_velo;

   public:
    void generateParticles(ContainerRef particles) override;
    CuboidGenerator(R3 position, R3 velocity, N3 num_particles, double mass, double distance, double avg_velo)
        : position(position),
          velocity(velocity),
          num_particles(num_particles),
          mass(mass),
          distance(distance),
          avg_velo(avg_velo) {};
};
}  // namespace mol_sim

#endif
