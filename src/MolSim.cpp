#include "io/CLIParse.h"
#include "io/fileReader/XVMReader.h"
#include "io/fileReader/YAMLReader.h"
#include "particles/container/SimpleContainer.h"
#include "physics/ForceSource.h"
#include "utils/Simulation.h"

using namespace mol_sim;

constexpr double START_TIME = 0;

int main(int argc, char* argsv[]) {
    double delta_t = 0.014;
    double end_time = 1000;
    YAMLReader file_reader;
    SimpleContainer particles;

    cliParse(argc, argsv, file_reader, delta_t, end_time, particles);
    Simulation<SimpleContainer> simulation(particles, GRAVITATIONAL, delta_t, START_TIME, end_time);
    simulation.run();
}
