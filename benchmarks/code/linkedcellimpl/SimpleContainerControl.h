#ifndef SIMPLE_CONTAINER_CONTROL_H
#define SIMPLE_CONTAINER_CONTROL_H

#include <particles/Particle.h>

#include <limits>
#include <set>
#include <vector>

#include "particles/ParticleContainer.h"
#include "particles/boundaries/Boundary.h"

namespace mol_sim {

/**
 * @brief Simple Container for Particles
 *
 * This container implements the concept ParticleContainer.
 *
 * Implemented as a SimpleContainer without the Newton's Third Law optimization for benchmarking purposes, to compare to
 * Linked-Cell-Implementations. Implements all methods of the SimpleContainer (so documentation is analogous), but only
 * used for benchmarking.
 *
 *
 */
class SimpleContainerControl : public std::vector<Particle> {
    R3 domain_size = R3{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity()};
    double cutoff_radius = std::numeric_limits<double>::infinity();

   public:
    // constructors
    using std::vector<Particle>::vector;

    SimpleContainerControl(R3 domain_size_arg, double cutoff_radius_arg);

    [[nodiscard]] bool fitsDomain(R3 v) const noexcept;

    void addParticle(Particle&& value);
    void addParticle(const Particle& value);
    void addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg);
    void addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg, int type);
    void addParticle(R3 x_arg, R3 old_x_arg, R3 v_arg, R3 f_arg, R3 old_f_arg, double m_arg, double epsilon_arg,
                     double sigma_arg, int type);

    std::vector<Particle>::iterator updateParticlePosition(std::vector<Particle>::iterator p, R3 new_x);
    std::vector<Particle>::iterator eraseParticle(std::vector<Particle>::iterator p);

    /**
     * @brief Iterator that iterates over all particles that apply a force on a given particle.
     */
    template <typename P>
        requires(std::is_same_v<P, Particle> || std::is_same_v<P, const Particle>)
    class proximity_iterator {
        P* cur;
        P* end;
        double radius;
        R3 center_or_domain;
        bool prox;
        std::set<BoundaryLocation> locations;

        bool fitBoundary() {  // NOLINT
            double effective_radius = -radius;
            if (cur->getX()[0] > center_or_domain[0] || cur->getX()[1] > center_or_domain[1] ||
                cur->getX()[2] > center_or_domain[2] || cur->getX()[0] < 0.0 || cur->getX()[1] < 0.0 ||
                cur->getX()[2] < 0.0) {
                return false;
            }
            for (auto location : locations) {
                switch (location) {
                    case BoundaryLocation::UPPER: {
                        if (cur->getX()[2] >= center_or_domain[2] - effective_radius &&
                            cur->getX()[2] <= center_or_domain[2]) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::LOWER: {
                        if (cur->getX()[2] >= 0.0 && cur->getX()[2] <= effective_radius) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::FRONT: {
                        if (cur->getX()[1] >= 0.0 && cur->getX()[1] <= effective_radius) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::BACK: {
                        if (cur->getX()[1] >= center_or_domain[1] - effective_radius &&
                            cur->getX()[1] <= center_or_domain[1]) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::LEFT: {
                        if (cur->getX()[0] >= 0.0 && cur->getX()[0] <= effective_radius) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::RIGHT: {
                        if (cur->getX()[0] >= center_or_domain[0] - effective_radius &&
                            cur->getX()[0] <= center_or_domain[0]) {
                            return true;
                        }
                        break;
                    }
                    default:
                        break;
                }
            }
            return false;
        }
        bool fitHalo() {  // NOLINT
            if (cur->getX()[0] > center_or_domain[0] + radius || cur->getX()[1] > center_or_domain[1] + radius ||
                cur->getX()[2] > center_or_domain[2] + radius || cur->getX()[0] < -radius || cur->getX()[1] < -radius ||
                cur->getX()[2] < -radius) {
                return false;
            }
            for (auto location : locations) {
                switch (location) {
                    case BoundaryLocation::UPPER: {
                        if (cur->getX()[2] > center_or_domain[2] && cur->getX()[2] <= center_or_domain[2] + radius) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::LOWER: {
                        if (cur->getX()[2] >= -radius && cur->getX()[2] < 0.0) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::FRONT: {
                        if (cur->getX()[1] >= -radius && cur->getX()[1] < 0.0) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::BACK: {
                        if (cur->getX()[1] > center_or_domain[1] && cur->getX()[1] <= center_or_domain[1] + radius) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::LEFT: {
                        if (cur->getX()[0] >= -radius && cur->getX()[0] < 0.0) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::RIGHT: {
                        if (cur->getX()[0] > center_or_domain[0] && cur->getX()[0] <= center_or_domain[0] + radius) {
                            return true;
                        }
                        break;
                    }
                    default:
                        break;
                }
            }
            return false;
        }
        void satisfy() {
            if (prox) {  // proximity check
                if (std::isinf(radius)) {
                    return;
                }
                while (cur != end && !((center_or_domain - cur->getX()).euclidNorm() <= radius)) {
                    ++cur;
                }
            } else {
                if (radius < 0.0) {  // boundary check
                    while (cur != end && !(fitBoundary())) {
                        ++cur;
                    }
                } else {  // halo check
                    while (cur != end && !(fitHalo())) {
                        ++cur;
                    }
                }
            }
        }

       public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = P;
        using difference_type = std::ptrdiff_t;
        using pointer = P*;
        using reference = P&;

        proximity_iterator() noexcept : cur(nullptr), end(nullptr), radius(0.0), prox(true) {}
        proximity_iterator(R3 center_or_domain, double radius, P* cur, P* end, size_t offset, bool prox = true,
                           std::set<BoundaryLocation> locations = {})
            : cur(cur + offset),
              end(end),
              radius(radius),
              center_or_domain(center_or_domain),
              prox(prox),
              locations(std::move(locations)) {
            satisfy();
        }

        reference operator*() const noexcept { return *cur; }
        pointer operator->() const noexcept { return cur; }

        proximity_iterator<P>& operator++() {
            ++cur;  // NOLINT
            satisfy();
            return *this;
        }

        proximity_iterator<P> operator++(int) {
            proximity_iterator<P> tmp = *this;
            ++(*this);
            return tmp;
        }

        friend bool operator==(const proximity_iterator<P>& a, const proximity_iterator<P>& b) noexcept {
            return a.cur == b.cur;
        }
        friend bool operator!=(const proximity_iterator<P>& a, const proximity_iterator<P>& b) noexcept {
            return !(a == b);
        }

        [[nodiscard]] R3 getCenter() const noexcept { return center_or_domain; }
        [[nodiscard]] double getRadius() const noexcept { return radius; }
        [[nodiscard]] bool isProximity() const noexcept { return prox; }
        [[nodiscard]] std::set<BoundaryLocation> getLocations() const noexcept { return locations; }
    };
    static_assert(std::forward_iterator<proximity_iterator<Particle>>);
    static_assert(std::forward_iterator<proximity_iterator<const Particle>>);

    proximity_iterator<Particle> eraseParticle(const proximity_iterator<Particle>& p);

    [[nodiscard]] proximity_iterator<Particle> proximityBegin(R3 center, size_t offset = 0);
    [[nodiscard]] proximity_iterator<const Particle> proximityBegin(R3 center, size_t offset = 0) const;
    [[nodiscard]] proximity_iterator<Particle> proximityEnd(R3 center);
    [[nodiscard]] proximity_iterator<const Particle> proximityEnd(R3 center) const;

    // boundary and halo iterators

    [[nodiscard]] proximity_iterator<Particle> haloBegin(const std::set<BoundaryLocation>& locations = {
                                                             BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                             BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                             BoundaryLocation::LEFT, BoundaryLocation::RIGHT});
    [[nodiscard]] proximity_iterator<const Particle> haloBegin(const std::set<BoundaryLocation>& locations = {
                                                                   BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                                   BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                                   BoundaryLocation::LEFT,
                                                                   BoundaryLocation::RIGHT}) const;
    [[nodiscard]] proximity_iterator<Particle> haloEnd(const std::set<BoundaryLocation>& locations = {
                                                           BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                           BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                           BoundaryLocation::LEFT, BoundaryLocation::RIGHT});
    [[nodiscard]] proximity_iterator<const Particle> haloEnd(const std::set<BoundaryLocation>& locations = {
                                                                 BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                                 BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                                 BoundaryLocation::LEFT,
                                                                 BoundaryLocation::RIGHT}) const;

    [[nodiscard]] proximity_iterator<Particle> boundaryBegin(const std::set<BoundaryLocation>& locations = {
                                                                 BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                                 BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                                 BoundaryLocation::LEFT, BoundaryLocation::RIGHT});
    [[nodiscard]] proximity_iterator<const Particle> boundaryBegin(const std::set<BoundaryLocation>& locations = {
                                                                       BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                                       BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                                       BoundaryLocation::LEFT,
                                                                       BoundaryLocation::RIGHT}) const;
    [[nodiscard]] proximity_iterator<Particle> boundaryEnd(const std::set<BoundaryLocation>& locations = {
                                                               BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                               BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                               BoundaryLocation::LEFT, BoundaryLocation::RIGHT});
    [[nodiscard]] proximity_iterator<const Particle> boundaryEnd(const std::set<BoundaryLocation>& locations = {
                                                                     BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                                     BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                                     BoundaryLocation::LEFT,
                                                                     BoundaryLocation::RIGHT}) const;

    void prepareForParallelIteration() const {}  // No-op for benchmark container
};
static_assert(ParticleContainer<SimpleContainerControl>);

}  // namespace mol_sim

#endif