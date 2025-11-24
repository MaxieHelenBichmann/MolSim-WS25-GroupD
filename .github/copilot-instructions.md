# MolSim AI Coding Guidelines

## Project Overview
MolSim is a C++20 molecular dynamics simulator using template metaprogramming with C++ concepts for polymorphic particle containers and force calculations. The project simulates particle interactions (gravitational, Lennard-Jones) with configurable parameters via YAML input files.

## File Structure
- `include/`: Public headers for core components (particles, physics, I/O, utils)
- `input/`: Sample input YAML files and format documentation
- `src/`: Implementation files for main simulation logic and CLI
- `tests/`: GoogleTest unit and integration tests
- `benchmarks/`: Google Benchmark performance tests
- `src/`, `tests/`, and `benchmarks/` all have their own CMakeLists.txt for modular builds

## Architecture

### Core Design Pattern: Compile-time Polymorphism
The codebase uses **C++ concepts instead of virtual functions** for performance-critical components:

- **ParticleContainer concept** (`include/particles/ParticleContainer.h`): Defines interface for containers like `SimpleContainer` and `LinkedCellContainer`
- **ForceSource concept** (`include/physics/ForceSource.h`): Defines interface for forces like `GravitationalForce` and `LennardJonesForce`
- **Simulation class** is templated: `template <ParticleContainer containerType, ForceSource forceType>`

**Why?** Concepts enable compile-time polymorphism, avoiding runtime virtual function overhead in tight simulation loops.

### Key Components
- **Particles**: `Particle` class stores position (x), velocity (v), force (f), old_force (old_f), mass, epsilon, sigma, type
- **Containers**: `SimpleContainer` (std::vector wrapper), `LinkedCellContainer` (spatial optimization with cells)
- **Forces**: `GravitationalForce`, `LennardJonesForce` - implement `applyForce(p1, p2) -> R3`
- **Simulation**: Main loop in `Simulation::run()` using Störmer-Verlet integration
- **I/O**: YAML input via `YAMLReader`, output via `VTKWriter`/`XYZWriter`
- **Generators**: `ParticleGenerator` for creating particle arrangements (cubes, discs)

### Namespace Convention
All project code lives in `namespace mol_sim`. Use `using namespace mol_sim;` in `.cpp` files only, never in headers.

## Build & Development Workflow

### Initial Setup
```bash
mkdir build && cd build
ccmake ..  # Configure options interactively
make -j $(nproc)
```

### Key CMake Options
- `ENABLE_TESTING=ON` - GoogleTest suite
- `ENABLE_BENCHMARK=ON` - Google Benchmark (requires `BENCHMARK_DOWNLOAD_DEPENDENCIES=ON`)
- `ENABLE_DOXYGEN=ON` - API documentation
- `COVERAGE=ON` - LLVM coverage reports
- `ENABLE_VTK_OUTPUT=ON` - VTK visualization output (requires VTK 8.9+)
- `CMAKE_BUILD_TYPE` - Debug (O0), Release (O3), RelWithDebInfo, MinSizeRel

### Running Simulations
```bash
# Lennard-Jones particles collision
./build/project/MolSim ./input/particles.yaml -d 0.0002 -t 5

# Gravitational simulation
./build/project/MolSim ./input/planets.yaml -d 0.014 -t 1000 -f GRAV
```

### Testing
```bash
./build/tests/tests
# or
ctest -V --test-dir ./build/tests
```


### Code Quality Tools
- **clang-format**: Uses Google style with 120 char line limit, 4-space indent (`.clang-format`)
- **clang-tidy**: Enabled by default, warnings as errors (`.clang-tidy`)
  - Run auto-fix: `make fix` (in build directory)
  - Enforces naming: `camelBack` functions, `lower_case` variables/members, `UPPER_CASE` constants, `CamelCase` enums

## Project-Specific Patterns

### Adding New Force Types
1. Create class in `include/physics/` implementing `ForceSource` concept
2. Implement `Vector<double, 3> applyForce(const Particle& p1, const Particle& p2)`
3. Add enum entry to `Force` in `physics/ForceSource.h`
4. Update switch statement in `src/MolSim.cpp` main function
5. Update CLI parser in `include/io/CLIParse.h` if needed

### Adding New Container Types
1. Create class in `include/particles/container/` satisfying `ParticleContainer` concept
2. Implement required methods: `operator[]`, `size()`, `empty()`, `clear()`, `reserve()`, `addParticle()`, iterators (`begin()`, `end()`), proximity iterators
3. See `SimpleContainer` for vector-based reference, `LinkedCellContainer` for spatial optimization

### YAML Input Format
All formats are documented in `input/formats`
Input files can have `settings:` block with `format: Settings`, `delta_t`, `end_time`, `start_time`. Particles defined via:
- `cube:` with `cuboids:` array (position, velocity, particle grid dimensions, mass, distance, mean_velo, epsilon, sigma)
- `disc:` with `discs:` array (similar structure for disc-shaped particle arrangements)

### Logging with spdlog
- Use `SPDLOG_INFO()`, `SPDLOG_WARN()`, `SPDLOG_ERROR()`, `SPDLOG_TRACE()`, `SPDLOG_DEBUG()`
- Log levels controlled by CMake build type (Debug/RelWithDebInfo: TRACE, Release/MinSizeRel: WARN)
- Formatted logging: `SPDLOG_INFO("Simulation: {} particles", count)`

### Testing Patterns
- Use `EXPECT_R3_EQ(expected, actual)` or `EXPECT_R3_NEAR(expected, actual, tolerance)` macros from `tests/testingUtils.h` for vector comparisons
- Default tolerance is `1e-9` for floating-point equality and 5e-8 for computational results

## CI/CD Pipeline
GitHub Actions workflow (`.github/workflows/ci.yml`):
- **build**: Clang, Release mode, clang-tidy, code coverage generation
- **sanitizer**: AddressSanitizer with Clang in Debug mode
- **documentation**: Doxygen generation and GitHub Pages deployment

Runs on custom `mini-arc-runner-set` with `mbichmann/molsim-vtk` Docker image.

## Performance Considerations
- Force calculations in `Simulation::calculateF()` use proximity iterators for spatial locality
- `LinkedCellContainer` provides O(n) force calculation vs O(n²) in `SimpleContainer` for large particle counts
- Cutoff radius limits force calculation range (default: infinity)

## Common Pitfalls
- Don't use virtual functions for containers/forces - breaks template design
- Keep proximity iteration consistent: always call `proximityBegin()` with position, cutoff radius, particle index
- YAML parsing errors caught as `YAMLReaderException` - handle in CLI parsing
- Vector type `R3` is alias for `Vector<double, 3>` (defined in `utils/Vector.h`)
