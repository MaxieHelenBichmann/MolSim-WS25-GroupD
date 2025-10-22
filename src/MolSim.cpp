
#include <iostream>
#include <memory>

#include "io/FileReader.h"
#include "io/outputWriter/VTKWriter.h"
#include "particles/container/SimpleContainer.h"
#include "physics/ForceSource.h"
#include "physics/GravitationalForce.h"

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

constexpr double START_TIME = 0;
constexpr double END_TIME = 1000;
constexpr double DELTA_T = 0.014;

// TODO: what data structure to pick?
particle_containers::SimpleContainer particles;
std::unique_ptr<ForceSource> forceSource = std::make_unique<GravitationalForce>();

int main(int argc, char* argsv[]) {
  std::cout << "Hello from MolSim for PSE!" << '\n';
  if (argc != 2) {
    std::cout << "Erroneous programme call! " << '\n';
    std::cout << "./molsym filename" << '\n';
  }

  FileReader file_reader;
  FileReader::readFile(particles, argsv[1]);

  double current_time = START_TIME;

  int iteration = 0;

  // for this loop, we assume: current x, current f and current v are known
  while (current_time < END_TIME) {
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
    std::cout << "Iteration " << iteration << " finished." << '\n';

    current_time += DELTA_T;
  }

  std::cout << "output written. Terminating..." << '\n';
  return 0;
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
  outputWriter::VTKWriter::plotParticles(particles, out_name, iteration);
}
