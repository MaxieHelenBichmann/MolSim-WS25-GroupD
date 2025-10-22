#include <iostream>

#include "io/fileReader/XVMReader.h"
#include "particles/container/SimpleContainer.h"
#include "physics/ForceSource.h"
#include "utils/Simulation.h"
#include "io/CLIParse.h"

constexpr double start_time = 0;
particle_containers::SimpleContainer particles;

int main(int argc, char* argsv[]) {
  double delta_t;
  double end_time;
  XVMReader filereader;

  cliParse(argc, argsv, filereader, delta_t, end_time, particles);
  Simulation simulation(particles, GRAVITATIONAL, delta_t);

  double current_time = start_time;
  int iteration = 0;

  // for this loop, we assume: current x, current f and current v are known
  while (current_time < end_time) {
    simulation.calculateX();
    simulation.calculateF();
    simulation.calculateV();

    iteration++;
    if (iteration % 10 == 0) {
      simulation.plotParticles(iteration);
    }
    std::cout << "Iteration " << iteration << " finished." << '\n';

    current_time += delta_t;
  }

  std::cout << "output written. Terminating..." << '\n';
  return 0;
}