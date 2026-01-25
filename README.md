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
- Multiple force implementations (Gravitational, Lennard-Jones, Harmonic)
- Single particle forces (Target forces for membrane simulations)
- Efficient particle containers (Direct Sum, Linked Cell)
- Boundary condition support (Reflecting, Outflow, Periodic)
- Particle generators (Cuboid, Disc, Membrane) with Brownian Motion
- Output formats (VTK, XYZ)
- Checkpointing support for long simulations
- Statistical Analysis of the Simulation over its runtime
- 2D Particle Membranes 
- Configurable via YAML input files

**Supported Simulations:**
- Gravitational N-body problems (planetary systems, stellar dynamics)
- Molecular dynamics with Lennard-Jones potentials (fluids, collisions)
- Membrane simulations with harmonic spring forces
- Large-scale particle systems using spatial optimization (Linked Cells)
- Argon Crystallization simulations

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
- Intel oneAPI toolkit

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
| `ENABLE_INTEL_VTUNE` | Enable Intel VTune profiling targets |
| `ENABLE_INTEL_ADVISOR` | Enable Intel Advisor profiling targets |
| `ENABLE_IO` | Enable IO Output |
| `ENABLE_CHECKPOINTS` | Enables Checkpoint writing at specified intervals |
| `ENABLE_STATS` | Enables writing of simulation statistics at specified intervals |
| `ENABLE_UBSAN` | Enables undefined behavior sanitizer (Debug only) |
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
| `Simulation/Full` | Full end-to-end simulation benchmarks and the Contest 1 benchmark  |
| `LinkedCell/` | LinkedCellContainer-specific operation benchmarks |

### Example Benchmark Commands

```bash
# Run only complexity benchmarks
./build/benchmarks/MolSimBench --benchmark_filter=Simulation/Complexity/

# Run cell implementation comparison
./build/benchmarks/MolSimBench --benchmark_filter=Cell/

# Run boundary benchmarks
./build/benchmarks/MolSimBench --benchmark_filter=Boundary/
```

### Python Analysis
You can use the analyse_data.py file in benchmarks/data to quickly get an overview of the data produced by the benchmarks. MatPlotLib and Pandas are required for this -> see benchmarks/data/requirements.txt
```bash
pip install requirements.txt

./build/benchmarks/MolSimBench --benchmark_filter=<filter> --benchmark_out_format=json --benchmark_out=<path/to/outputfile>

python benchmarks/data/analyse_data.py <path/to/outputfile>
```

## Analyzing the Project
This project being a Molecular Dynamics Simulator, with `ENABLE_STATS` set and configured via the input file, you can additionally collect thermodynamical statistics to further analyze the simulation.

### Collecting Data
There are two supported statistics - the Diffusion and Radial Distribution Function. The required configurations can be seen in input/formats/SettingsFormat.md. While running, the Simulator will collect the necessary statistics in `diffusion.csv` or `rdf.csv`.

### Python Analysis
You can use the plot_diffusion.py or plot_rdf.py file in scripts/statistics to visualize the collected data. MatPlotLib is required for this -> see scripts/statistics/requirements.txt
```bash
cd ./scripts/statistics

pip install -r requirements.txt

# Plot data collected for Diffusion
python plot_diffusion.py <path/to/diffusion.csv> --delta-t <DELTA_T> --start-time <START_TIME> --out <FILE>.png

# Plot data of Radial Distribution Function
python plot_rdf.py <path/to/rdf.csv> --delta-t <DELTA_T> --start-time <START_TIME> --out <FILE>.png --y0
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

You can enable treating clang-tidy warnings as errors by configuring CMake with `-DCLANG_TIDY_WARNINGS_AS_ERRORS=ON` (default: `OFF`). The CI `lint` job enables this option.

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
| `vtune` | Intel VTune Profiler (optional, requires Intel oneAPI) |
| `advisor` | Intel Advisor (optional, requires Intel oneAPI) |


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

### Intel VTune and Advisor

> **Note:** Requires `ENABLE_INTEL_VTUNE` and/or `ENABLE_INTEL_ADVISOR` in CMake configuration and Intel oneAPI toolkit installed.

**VTune Profiler targets:**
```bash
make vtune-hotspots    # CPU hotspot analysis
make vtune-memory      # Memory access patterns
make vtune-uarch       # Microarchitecture exploration
make vtune-all         # Run all VTune analyses
```

**Advisor targets:**
```bash
make advisor-survey       # Vectorization opportunities
make advisor-tripcounts   # Loop iteration counts
make advisor-roofline     # Roofline performance model
make advisor-all          # Run all Advisor analyses
```

**Combined Intel profiling:**
```bash
make intel-profile-all  # Run all analyses and create archive
```

Results saved to `build/vtune/`, `build/advisor/`, or `build/profiling-results/`.

---


