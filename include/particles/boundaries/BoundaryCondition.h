#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

#include <optional>
#include <vector>

#include "particles/Particle.h"
#include "particles/ParticleContainer.h"
namespace mol_sim {
/**
 * @brief Enum for boundary types in the Linked-Cell Container.
 *
 * UPPER: +z direction (x-y plane at max z)
 * LOWER: -z direction (x-y plane at min z)
 * FRONT: -y direction (x-z plane at min y)
 * BACK: +y direction (x-z plane at max y)
 * LEFT: -x direction (y-z plane at min x)
 * RIGHT: +x direction (y-z plane at max x)
 */
enum class BoundaryType : std::uint8_t { UPPER, LOWER, FRONT, BACK, LEFT, RIGHT };
enum class BoundaryConditionType : std::uint8_t { OUTFLOW, REFLECTING, PERIODIC };

class BoundaryConditionDeclaration {
   public:
    BoundaryConditionDeclaration& setType(BoundaryConditionType type);
    BoundaryConditionDeclaration& setLocation(BoundaryType location);
    BoundaryConditionDeclaration& setCounterSigma(double counter_sigma);
    BoundaryConditionDeclaration& setCounterEpsilon(double counter_epsilon);

    BoundaryConditionType getType();
    BoundaryType getLocation();
    std::optional<double> getCounterSigma();
    std::optional<double> getCounterEpsilon();

    BoundaryConditionDeclaration() = default;

    BoundaryConditionDeclaration(BoundaryType location, BoundaryConditionType type);

    ~BoundaryConditionDeclaration() = default;

   protected:
    BoundaryConditionType type;
    BoundaryType location;
    std::optional<double> counter_sigma = std::nullopt;
    std::optional<double> counter_epsilon = std::nullopt;
};

template <ParticleContainer containerType>
class BoundaryCondition : public BoundaryConditionDeclaration {
    virtual bool boundaryConditionApplies(Particle& p) = 0;
    virtual void boundaryStrategy(Particle& p) = 0;

   protected:
    containerType* particles;
    BoundaryType location;

   public:
    BoundaryCondition(const BoundaryCondition& other) {
        particles = other.particles;
        location = other.location;
        counter_sigma = other.counter_sigma;
        counter_epsilon = other.counter_epsilon;
    }
    BoundaryCondition(containerType& particles, BoundaryType location);
    virtual ~BoundaryCondition();
    /***
     * @brief Iterates over the boundary cells and applies the boundary condtion if necessary.
     */
    void applyBoundary();
};

}  // namespace mol_sim

#endif