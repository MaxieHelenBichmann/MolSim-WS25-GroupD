#ifndef OUTFLOW_H
#define OUTFLOW_H

#include "particles/boundaries/Boundary.h"

namespace mol_sim {

class Outflow : public Boundary {
   public:
    explicit Outflow(BoundaryLocation location) : Boundary(location, BoundaryType::OUTFLOW) {}
    ~Outflow() override = default;

    [[nodiscard]] std::optional<Particle> applyBoundary(Particle& p) override {
        (void)p;
        return std::nullopt;
    }
};

}  // namespace mol_sim

#endif
