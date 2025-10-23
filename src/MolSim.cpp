#include <iostream>

#include "io/CLIParse.h"
#include "io/fileReader/XVMReader.h"
#include "particles/container/SimpleContainer.h"
#include "physics/ForceSource.h"
#include "utils/Simulation.h"

constexpr double START_TIME = 0;
particle_containers::SimpleContainer particles;

int main(int argc, char* argsv[]) {
  double DELTA_T;
  double END_TIME;
  XVMReader file_reader;

  cliParse(argc, argsv, file_reader, DELTA_T, END_TIME, particles);
  Simulation simulation(particles, GRAVITATIONAL, DELTA_T);

  double current_time = START_TIME;
  int iteration = 0;

  // for this loop, we assume: current x, current f and current v are known
  while (current_time < END_TIME) {
    // calculate new x
    simulation.calculateX();
    // calculate new f
    simulation.calculateF();
    // calculate new v
    simulation.calculateV();

    iteration++;
    if (iteration % 10 == 0) {
      simulation.plotParticles(iteration);
    }
    std::cout << "Iteration " << iteration << " finished." << '\n';

    current_time += DELTA_T;
  }

  std::cout << "output written. Terminating..." << '\n';
  return 0;
}


