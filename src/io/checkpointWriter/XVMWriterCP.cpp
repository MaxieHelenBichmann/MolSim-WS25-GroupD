#include "io/checkpointWriter/XVMWriterCP.h"

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <sstream>

using namespace mol_sim;

XVMWriterCP::XVMWriterCP() = default;

XVMWriterCP::~XVMWriterCP() = default;

void XVMWriterCP::createCheckpoint([[maybe_unused]] const Domain& domain, ContainerRef particles, int iteration,
                                   [[maybe_unused]] Force force, [[maybe_unused]] double delta_t,
                                   [[maybe_unused]] double start_time, [[maybe_unused]] double end_time,
                                   [[maybe_unused]] size_t frequency_output,
                                   [[maybe_unused]] size_t frequency_checkpoint,
                                   [[maybe_unused]] const std::string& base_name, [[maybe_unused]] double cutoff_radius,
                                   size_t N) const {
    int decimal_places = 0;
    while (N >= 10) {
        N /= 10;
        decimal_places++;
    }
    std::ofstream file;
    std::stringstream strstr;
    strstr << "cp_" << std::setfill('0') << std::setw(decimal_places) << iteration << ".txt";

    file.open(strstr.str().c_str());

    // Particles
    file << "# Number Particles" << '\n';
    file << particles.size() << '\n';
    file << "# Position | Velocity | Mass" << '\n';
    for (auto& p : particles) {
        file << p.getX()[0] << " " << p.getX()[1] << " " << p.getX()[2] << " ";
        file << p.getV()[0] << " " << p.getV()[1] << " " << p.getV()[2] << " ";
        file << p.getM() << '\n';
    }
    file.close();
}