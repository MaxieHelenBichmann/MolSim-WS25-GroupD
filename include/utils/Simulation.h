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
  std::unique_ptr<ForceSource> forceSource;
  double DELTA_T;

 public:
  Simulation(containerType& particles, Force forceType, double DELTA_T) : particles(particles), DELTA_T(DELTA_T) {
    switch (forceType) {
      case GRAVITATIONAL:
        this->forceSource = std::make_unique<GravitationalForce>();
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
        p1.getF() = p1.getF() + forceSource->calculateForce(p1, p2);
      }
    }
  }

  void calculateX() {
    for (auto& p : particles) {
      p.getX() = p.getX() + DELTA_T * p.getV() + (0.5 * DELTA_T * DELTA_T / p.getM()) * p.getF();
    }
  }

  void calculateV() {
    for (auto& p : particles) {
      p.getV() = p.getV() + (0.5 * DELTA_T / p.getM()) * (p.getOldF() + p.getF());
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