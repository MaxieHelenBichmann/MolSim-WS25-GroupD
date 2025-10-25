#ifndef SIMULATION_H
#define SIMULATION_H

#include <memory>

#include "io/outputWriter/VTKWriter.h"
#include "io/outputWriter/XYZWriter.h"
#include "particles/ParticleContainer.h"
#include "physics/ForceSource.h"
#include "physics/GravitationalForce.h"

namespace mol_sim {

template <ParticleContainer containerType>
class Simulation {
 private:
  containerType& particles;
  std::unique_ptr<ForceSource> force_source;
  double delta_t;

 public:
  Simulation(containerType& particles, Force forceType, double DELTA_T) : particles(particles), delta_t(DELTA_T) {
    switch (forceType) {
      case GRAVITATIONAL:
        this->force_source = std::make_unique<GravitationalForce>();
        break;
      default:
        break;
    }
  }

  void calculateF() {
    for (auto& p1 : particles) {
      p1.getOldF() = p1.getF();
      p1.getF() = Vector<double, 3>();
      for (auto& p2 : particles) {
        if (p1 == p2) {
          continue;
        }
        p1.getF() = p1.getF() + force_source->calculateForce(p1, p2);
      }
    }
  }

  void calculateX() {
    for (auto& p : particles) {
      p.getX() = p.getX() + delta_t * p.getV() + (0.5 * delta_t * delta_t / p.getM()) * p.getF();
    }
  }

  void calculateV() {
    for (auto& p : particles) {
      p.getV() = p.getV() + (0.5 * delta_t / p.getM()) * (p.getOldF() + p.getF());
    }
  }

  void plotParticles(int iteration) {
    std::string out_name("MD_vtk");
#ifdef ENABLE_VTK_OUTPUT
    VTKWriter writer;
#else
    XYZWriter writer;
#endif
    writer.plotParticles(particles, out_name, iteration);
  }
};

}  // namespace mol_sim

#endif