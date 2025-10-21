
#include <iostream>
#include <list>
#include <memory>

#include "io/FileReader.h"
#include "io/outputWriter/VTKWriter.h"
#include "particles/container/SimpleContainer.h"
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
particle_containers::SimpleContainer particles;
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
  auto iterator = particles.begin();

  for (auto& p1 : particles) {
    p1.getOldF() = p1.getF();
    p1.getF() = Vector<double, 3>(0.);
    for (auto& p2 : particles) {
      if (p1 == p2) {
        continue;
      }
      p1.getF() = forceSource->calculateForce(p1, p2);
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

  outputWriter::VTKWriter writer;
  writer.plotParticles(particles, out_name, iteration);
}
