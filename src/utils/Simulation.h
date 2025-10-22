#pragma once

#include <memory>

#include "io/outputWriter/XYZWriter.h"
#include "io/outputWriter/VTKWriter.h"
#include "particles/ParticleContainer.h"
#include "physics/ForceSource.h"
#include "physics/GravitationalForce.h"

class Simulation {
 private:
  particle_containers::ParticleContainer& particles;
  std::unique_ptr<ForceSource> forceSource;
  double delta_t;

 public:
  Simulation(particle_containers::ParticleContainer& particles, Force forceType, double delta_t)
      : particles(particles), delta_t(delta_t) {
    switch(forceType) {
      case GRAVITATIONAL:
        this->forceSource = std::make_unique<GravitationalForce>();
        break;
      default: break;
    }
  }

  void calculateF() {
    for (auto& p1 : particles) {
      p1.getOldF() = p1.getF();
      p1.getF() = Vector<double, 3>(0.);
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
    outputWriter::VTKWriter writer;
    #else
    outputWriter::XYZWriter writer;
    #endif
    writer.plotParticles(particles, out_name, iteration);
  }
};
