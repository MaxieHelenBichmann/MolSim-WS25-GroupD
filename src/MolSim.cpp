
#include <iostream>
#include <list>
#include <memory>

#include "io/FileReader.h"
#include "io/outputWriter/VTKWriter.h"
#include "physics/ForceSource.h"
#include "physics/GravitationalForce.h"
#include "utils/ArrayUtils.h"

/**** forward declaration of the calculation functions ****/

/**
 * calculate the force for all particles
 */
void calculateF();

/**
 * calculate the position for all particles
 */
void calculateX();

/**
 * calculate the position for all particles
 */
void calculateV();

/**
 * plot the particles to a xyz-file
 */
void plotParticles(int iteration);

constexpr double start_time = 0;
constexpr double end_time = 1000;
constexpr double delta_t = 0.014;

// TODO: what data structure to pick?
std::list<Particle> particles;
std::unique_ptr<ForceSource> forceSource = std::make_unique<GravitationalForce>();

int main(int argc, char* argsv[]) {
  std::cout << "Hello from MolSim for PSE!" << std::endl;
  if (argc != 2) {
    std::cout << "Erroneous programme call! " << std::endl;
    std::cout << "./molsym filename" << std::endl;
  }

  FileReader fileReader;
  fileReader.readFile(particles, argsv[1]);

  double current_time = start_time;

  int iteration = 0;

  // for this loop, we assume: current x, current f and current v are known
  while (current_time < end_time) {
    // calculate new x
    calculateX();
    // calculate new f
    calculateF();
    // calculate new v
    calculateV();

    iteration++;
    if (iteration % 10 == 0) {
      plotParticles(iteration);
    }
    std::cout << "Iteration " << iteration << " finished." << std::endl;

    current_time += delta_t;
  }

  std::cout << "output written. Terminating..." << std::endl;
  return 0;
}

void calculateF() {
  std::list<Particle>::iterator iterator;
  iterator = particles.begin();

  for (auto& p1 : particles) {
    p1.modOldF() = {p1.getF()[0], p1.getF()[1], p1.getF()[2]};
    p1.modF() = {0., 0., 0.};
    for (auto& p2 : particles) {
      if (p1 == p2) {
        continue;
      }
      p1.modF() = forceSource->calculateForce(p1, p2);
    }
  }
}

void calculateX() {
  for (auto& p : particles) {
    p.modX() = {p.getX()[0] + delta_t * p.getV()[0] + 0.5 * delta_t * delta_t * p.getF()[0] / p.getM(),
                p.getX()[1] + delta_t * p.getV()[1] + 0.5 * delta_t * delta_t * p.getF()[1] / p.getM(),
                p.getX()[2] + delta_t * p.getV()[2] + 0.5 * delta_t * delta_t * p.getF()[2] / p.getM()};
  }
}

void calculateV() {
  for (auto& p : particles) {
    p.modV() = {
        p.getV()[0] + 0.5 * delta_t * (p.getOldF()[0] + p.getF()[0]) / p.getM(),
        p.getV()[1] + 0.5 * delta_t * (p.getOldF()[1] + p.getF()[1]) / p.getM(),
        p.getV()[2] + 0.5 * delta_t * (p.getOldF()[2] + p.getF()[2]) / p.getM(),
    };
  }
}

void plotParticles(int iteration) {
  std::string out_name("MD_vtk");

  outputWriter::VTKWriter writer;
  writer.plotParticles(particles, out_name, iteration);
}
