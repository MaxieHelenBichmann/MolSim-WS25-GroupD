#ifndef GRAV_FORCE
#define GRAV_FORCE

#include "physics/ForceSource.h"
#include "utils/ArrayUtils.h"

class GravitationalForce : public ForceSource {
    public:

    [[nodiscard]] std::array<double, 3> calculateForce(const Particle& p1, const Particle& p2) const override;

};

#endif
