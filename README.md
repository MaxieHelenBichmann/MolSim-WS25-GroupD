# MolSim Group D

[![CI MolSim](https://github.com/MaxieHelenBichmann/MolSim-WS25-GroupD/actions/workflows/ci.yml/badge.svg)](https://github.com/MaxieHelenBichmann/MolSim-WS25-GroupD/actions/workflows/ci.yml)

**Contributors:**
- Maxie Helen Bichmann
- Georg Sebastian Eisner
- Henry Jacob Meyran

**Documentation:** [https://maxiehelenbichmann.github.io/MolSim-WS25-GroupD/](https://maxiehelenbichmann.github.io/MolSim-WS25-GroupD/)

---

## Table of Contents

- [Overview](#overview)
- [Project Structure](#project-structure)
- [Dependencies](#dependencies)
- [Building the Project](#building-the-project)
- [Running the Project](#running-the-project)
- [Documentation](#documentation)
- [Testing the Project](#testing-the-project)
- [Benchmarking the Project](#benchmarking-the-project)
- [Optional Tools](#optional-tools)
- [Profiling](#profiling)

---

## Overview

MolSim is a molecular dynamics simulation framework developed as part of the PSE Molekulardynamik course. The application simulates the physical behavior of particle systems using numerical integration methods and various force models.

**Key Features:**
- Multiple force implementations (Gravitational, Lennard-Jones)
- Efficient particle containers (Direct Sum, Linked Cell)
- Boundary condition support (Reflecting, Outflow, VelocityReflect)
- Particle generators utilities (Cuboid, Disc) with Brownian Motion
- Output formats (VTK, XYZ)
- Configurable via YAML input files

**Supported Simulations:**
- Gravitational N-body problems (planetary systems, stellar dynamics)
- Molecular dynamics with Lennard-Jones potentials (fluids, collisions)
- Large-scale particle systems using spatial optimization (Linked Cells)

---

## Project Structure

```
MolSim-WS25-GroupD/
├── src/                        # Source code
│   └── MolSim.cpp              # Main application entry point
├── include/                    # Public headers
├── tests/                      # Unit and integration tests
├── benchmarks/                 # Performance benchmarks
├── input/                      # Example simulation configurations 
│   └── formats/                # Explanation of all available input formats
├── build/                      # Build artifacts (generated)
└── doxys_documentation/        # Generated documentation (generated)
```

---

## Dependencies

### Essential
- CMake Version 3.10+
- Make 4.3+
- C++ compiler with C++20 support (tested with clang 18.1.3)

### VTK Output
- VTK Version 8.9+

### Optional Tools
- clang-format
- clang-tidy

### Doxygen Support
- Doxygen
- Graphviz

### Profiling
- Perf
- Valgrind

---

## Building the Project

**1. Navigate to project root**

**2. Configure with CMake**
```bash
mkdir build && cd build
ccmake ..
```

**3. Configure build options interactively**

<details>
<summary><b>CMake Configuration Options</b></summary>

| Option | Description |
|--------|-------------|
| `ENABLE_TESTING` | Enable all tests |
| `ENABLE_BENCHMARK` | Enable benchmarking (requires `BENCHMARK_DOWNLOAD_DEPENDENCIES`) |
| `ENABLE_DOXYGEN` | Enable Doxygen documentation generation |
| `ENABLE_PROFILING` | Enable profiling tools (perf and valgrind) |
| `COVERAGE` | Enable code coverage reports |
| `ENABLE_VTK_OUTPUT` | Enable output in VTK format |
| `CMAKE_BUILD_TYPE` | Build type: `Debug`, `Release`, `RelWithDebInfo`, `MinSizeRel` |

</details>

**4. Build the executable**
```bash
make -j $(nproc)
```

---

## Running the Project

From the project root:

**Planets simulation:**
```bash
./build/project/MolSim ./input/planets.yaml 
```

**Particles simulation:**
```bash
./build/project/MolSim ./input/particles.yaml
```
**Basin drop simulation**
```bash
./build/project/MolSim ./input/basin_drop.yaml
```
---

## Documentation

**Online API Documentation:**

Complete API and input format documentation is available at:
- **[https://maxiehelenbichmann.github.io/MolSim-WS25-GroupD/](https://maxiehelenbichmann.github.io/MolSim-WS25-GroupD/)**

**Local Documentation:**

You can also generate the documentation locally using Doxygen (see [Optional Tools](#optional-tools) section).

**Input File Format Reference:**

Additional format documentation files are available in `input/formats/`:
- `SettingsFormat.md` - Simulation settings configuration
- `CuboidFormat.md` - Cuboid particle generation
- `DiscFormat.md` - Disc particle generation
- `XVMFormat.md` - Individual particle specification

---

## Testing the Project

> **Note:** Ensure `ENABLE_TESTING` is enabled in CMake configuration.

**Run tests directly:**
```bash
./build/tests/tests
```

**Run with CTest (verbose output):**
```bash
ctest -V --test-dir ./build/tests
```

## Benchmarking the Project

> **Note:** Ensure `ENABLE_BENCHMARK` and `BENCHMARK_DOWNLOAD_DEPENDENCIES` are enabled in CMake configuration.

### Running Benchmarks

**Run all benchmarks:**
```bash
./build/benchmarks/MolSimBench
```

**Run specific benchmark suite:**
```bash
./build/benchmarks/MolSimBench --benchmark_filter=<pattern>
```

**Save results to file:**
```bash
./build/benchmarks/MolSimBench --benchmark_format=console > benchmark_results.txt
```

### Available Benchmark Suites

| Filter Pattern | Description |
|----------------|-------------|
| `Simulation/Complexity/` | Compares O(n) LinkedCell vs O(n²) DirectSum scaling |
| `Cell/` | Compares cell data structures (Vector, Set, UnorderedSet) |
| `Boundary/` | Measures boundary condition overhead (Reflecting, VelocityReflect) |
| `Simulation/` | Full end-to-end simulation benchmarks |(Do Not Work Currently)
| `LinkedCell/` | LinkedCellContainer-specific operation benchmarks | (Do Not Work Currently)

### Example Benchmark Commands

```bash
# Run only complexity benchmarks
./build/benchmarks/MolSimBench --benchmark_filter=Simulation/Complexity/

# Run cell implementation comparison
./build/benchmarks/MolSimBench --benchmark_filter=Cell/

# Run boundary benchmarks
./build/benchmarks/MolSimBench --benchmark_filter=Boundary/
```

---

## Optional Tools

### Doxygen Documentation

Generate documentation (requires `ENABLE_DOXYGEN`):
```bash
make doc_doxygen
```
Output location: `doxys_documentation/`

---

### Clang-Tidy

When enabled, clang-tidy automatically runs checks from `.clang-tidy` during build.

**Apply automatic fixes:**
```bash
make fix
```

---

### Test Coverage

Generate coverage reports (requires `COVERAGE` enabled):
```bash
make coverage_report
```
Output location: `build/coverage/`

---

## Profiling

### Prerequisites
| Tool | Purpose |
|------|---------|
| `perf` | CPU profiling (requires Linux kernel support) |
| `valgrind` | Memory debugging and profiling suite |
| `ms_print` | View massif output (included with valgrind) |


### Build Configuration for Profiling

| Profiling Type | Recommended Build Type | 
|----------------|------------------------|
| **Performance (perf)** | `RelWithDebInfo` |
| **Memory (valgrind)** | `Debug` or `RelWithDebInfo` |

### Customizing Profiling Inputs

<details>
<summary><b>Customize input files and arguments for profiling runs</b></summary>


| Option | Description |
|--------|-------------|
| `PERF_INPUT_FILE` | Input file used by all perf calls |
| `PERF_ARGS` | Extra arguments passed to perf calls |
| `VALGRIND_INPUT_FILE` | Input file used by all valgrind calls |
| `VALGRIND_ARGS` | Extra arguments passed to valgrind calls |
</details>  

> **Note:** Valgrind runs much slower, so use shorter simulations for memory checks.

### Perf Profiling

**Recommended Workflow:**
```bash
# 1. Build with profiling enabled
cd build && ccmake ..  # Set ENABLE_PROFILING=ON, CMAKE_BUILD_TYPE=RelWithDebInfo
make -j $(nproc)

# 2. Record performance data
make perf-record

# 3. Analyze results interactively
make perf-report

# 4. Check hardware counters
make perf-stat

# 5. Investigate cache issues if needed
make perf-cache
```

---

### Valgrind Memory Analysis

<details>
<summary><b>Available Targets</b></summary>

#### Memory Check
```bash
make valgrind-memcheck
```


Results: `build/valgrind/memcheck.log`

#### Heap Profiling
```bash
make valgrind-massif
```
Tracks heap memory usage over time.

**Analyze results:**
```bash
ms_print build/valgrind/massif.out
```

</details>

**Recommended Workflow:**
```bash
# 1. Build with profiling and debug symbols
cd build && ccmake ..  # Set ENABLE_PROFILING=ON, CMAKE_BUILD_TYPE=Debug
make -j $(nproc)

# 2. Check for memory leaks
make valgrind-memcheck
cat build/valgrind/memcheck.log  # Review findings

# 3. Profile heap usage
make valgrind-massif
ms_print build/valgrind/massif.out | less
```

---


