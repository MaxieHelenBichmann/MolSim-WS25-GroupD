#include <iostream>

#include "io/CLIParse.h"
#include "io/fileReader/XVMReader.h"
#include "particles/container/SimpleContainer.h"
#include "physics/ForceSource.h"
#include "utils/Simulation.h"

using namespace mol_sim;

constexpr double START_TIME = 0;
SimpleContainer particles;

int main(int argc, char* argsv[]) {
  double delta_t = 0.014;
  double end_time = 1000;
  XVMReader file_reader;

  cliParse(argc, argsv, file_reader, delta_t, end_time, particles);
  Simulation<SimpleContainer> simulation(particles, GRAVITATIONAL, delta_t);

  double current_time = START_TIME;
  int iteration = 0;

  // for this loop, we assume: current x, current f and current v are known
  while (current_time < end_time) {
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

    current_time += delta_t;
  }

  std::cout << "output written. Terminating..." << '\n';
  return 0;
}
