#include "io/checkpointWriter/XVMWriterCP.h"

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string_view>

#include "particles/container/SimpleContainer.h"

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

    // Settings
    file.open(strstr.str().c_str());
    file << "# Domain" << '\n';
    file << domain.getDimension()[0] << " " << domain.getDimension()[1] << " " << domain.getDimension()[2] << '\n';

    const auto boundary_type_string = [](BoundaryType type) -> std::string_view {
        switch (type) {
            case BoundaryType::OUTFLOW:
                return "OUTFLOW";
            case BoundaryType::REFLECTING:
                return "REFLECTING";
            case BoundaryType::VELOCITYREFLECT:
                return "VELOCITYREFLECT";
        }
        return "OUTFLOW";
    };

    file << "# Boundaries" << '\n';
    file << "upper " << boundary_type_string(domain.getBoundary(BoundaryLocation::UPPER).getType()) << '\n';
    file << "lower " << boundary_type_string(domain.getBoundary(BoundaryLocation::LOWER).getType()) << '\n';
    file << "left " << boundary_type_string(domain.getBoundary(BoundaryLocation::LEFT).getType()) << '\n';
    file << "right " << boundary_type_string(domain.getBoundary(BoundaryLocation::RIGHT).getType()) << '\n';
    file << "front " << boundary_type_string(domain.getBoundary(BoundaryLocation::FRONT).getType()) << '\n';
    file << "back " << boundary_type_string(domain.getBoundary(BoundaryLocation::BACK).getType()) << '\n';

    file << "# Force" << '\n';

    file << "# delta_t" << '\n';
    file << delta_t << '\n';

    file << "# start_time" << '\n';
    file << start_time << '\n';

    file << "# end_time" << '\n';
    file << end_time << '\n';

    file << "# frequency_output" << '\n';
    file << frequency_output << '\n';

    file << "# frequency_checkpoint" << '\n';
    file << frequency_checkpoint << '\n';

    file << "# base_name" << '\n';
    file << base_name << '\n';

    file << "# cutoff_radius" << '\n';
    file << cutoff_radius << '\n';

    // Particles
    file << "# Container Type" << '\n';
    file << (std::holds_alternative<SimpleContainer*>(particles.getInstance()) ? "SIMPLE" : "LINKED") << '\n';
    file << "# Particles" << '\n';
    file << particles.size() << '\n';
    file << "# Position | Velocity | Force | Old Force | Mass | Sigma | Epsilon | Type " << '\n';
    for (auto& p : particles) {
        file << p.getX()[0] << " " << p.getX()[1] << " " << p.getX()[2] << " ";
        file << p.getV()[0] << " " << p.getV()[1] << " " << p.getV()[2] << " ";
        file << p.getF()[0] << " " << p.getF()[1] << " " << p.getF()[2] << " ";
        file << p.getOldF()[0] << " " << p.getOldF()[1] << " " << p.getOldF()[2] << " ";
        file << p.getM() << " ";
        file << p.getSigma() << " ";
        file << p.getEpsilon() << " ";
        file << p.getType();
        file << '\n';
    }

    file.close();
}