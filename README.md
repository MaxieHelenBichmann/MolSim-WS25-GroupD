MolSim Group D
===

Contributors:

- Maxie Helen Bichmann
- Georg Sebastian Eisner
- Henry Jacob Meyran

[![CI MolSim](https://github.com/MaxieHelenBichmann/MolSim-WS25-GroupD/actions/workflows/ci.yml/badge.svg)](https://github.com/MaxieHelenBichmann/MolSim-WS25-GroupD/actions/workflows/ci.yml)

## Dependencies

- **Essential:**
    - CMake Version 3.10+
    - Make 4.3+
    - C++ compiler with C++20 support (tested with g++ 14.2.0)
- **VTK Output:**
    - VTK Version 8.9+
- **Optional Tools:**
    - clang-format
    - clang-tidy
- **Doxygen Support:**
    - Doxygen
    - Graphviz 
- **Profiling**
    - Perf
    - Valgrind

## Building the project

**How do you build the project?**

0. Make sure to be in the project root.

1. Configure the project with CMake by running these commands:

    ```
    mkdir build && cd build
    ccmake ..
    ```

2. Now you can interactively change the options.
    <details><summary>Configure Options</summary>
    - <b>ENABLE_TESTING</b> enables all tests.<br />
    - <b>ENABLE_BENCHMARK</b> enables benchmarking. (<b>BENCHMARK_DOWNLOAD_DEPENDENCIES</b> should also be enabled)<br />
    - <b>ENABLE_DOXYGEN</b> enables doxygen.<br />
    - <b>ENABLE_PROFILING</b> enables profiling tools (perf and valgrind).<br />
    - <b>COVERAGE</b> enables code coverage reports.<br />
    - <b>ENABLE_VTK_OUTPUT</b> enables output in the vtk format.<br />
    - <b>CMAKE_BUILD_TYPE</b> specifies the build type (Debug, Release,...).<br />
   </details>  

3. Build the executable with the Makefile by running the command:

    ```
    make -j $(nproc)
    ```

## Running the project

**How do I run the project?**

0. Make sure to be in the project root after you built the project.

1. Run the executable with the given input file and optional desired delta_T (-d), end_T (-t) and force (-f):

    ```
    ./build/project/MolSim ./input/planets.yaml -d 0.014 -t 1000 -f GRAV
    ```
    ```
    ./build/project/MolSim ./input/particles.yaml -d 0.0002 -t 5
    ```

## Testing the project

**How do I test the project?**

0. Make sure to be in the project root and to have built the project with the according CMake configuration.

1. Run the test executable or use ctest:
    ```
    ./build/tests/tests
    ```
    ```
    ctest -V --test-dir ./build/tests 
    ```

## Benchmarking the project

**How do I benchmark the project?**

0. Make sure to be in the project root and to have built the project with the according CMake configuration.

1. Run the benchmark executable:
    ```
    ./build/benchmarks/MolSimBench
    ```

## Optional Tools

**Doxygen**  
If `doxygen` support was enabled:

```
make doc_doxygen
```

Creates the `doxygen` documentation in the doxys_documentation directory.

**Clang-tidy**  
If `clang-tidy` support was enabled, `clang-tidy` will automatically run the checks specified in the .clang-tidy file.  
There is also:

```
make fix
```

Which automatically applies fixes where `clang-tidy` finds them.

**Test Coverage**  
If coverage was enabled:

```
make coverage_report
```

Creates a coverage report on the coverage directory inside the build directory.

### Profiling

If profiling is enabled both `perf` (for CPU profiling) and `valgrind` (for memory analysis) are available.

It is recommended to run `perf` on a `Release` or `RelWithDebInfo` build for accurate performance measurements.  
`Valgrind` relies on debug symblos (`-g`) so only use it with `Debug` or `RelWithDebInfo`.
 
The input and arguments for `valgrind` and `perf` can be customized using CMake Arguments `PERF_INPUT_FILE`. `PERF_ARGS`, `VALGRIND_INPUT_FILE` and `VALGRIND_ARGS`.

**Note:** Valgrind runs much slower, so use shorter simulations for memory checks.

#### Perf Profiling

**Available Targets:**

1. **perf-record**: Record performance data with call graph information
   ```bash
   make perf-record
   ```
   Creates `perf.data` in the build directory with detailed sampling information.

2. **perf-report**: Interactive analysis of recorded data
   ```bash
   make perf-report
   ```
   Opens an interactive TUI to explore hotspots and call chains.

3. **perf-report-text**: Generate text-based performance report
   ```bash
   make perf-report-text
   ```
   Saves analysis to `build/perf-report.txt` for review or sharing.

4. **perf-stat**: Display hardware counter statistics
   ```bash
   make perf-stat
   ```
   Shows CPU cycles, instructions, cache misses, and other hardware events.

5. **perf-cache**: Analyze cache performance in detail
   ```bash
   make perf-cache
   ```
   Reports cache references, misses, and L1 cache behavior.

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

#### Valgrind Memory Analysis

**Available Targets:**

1. **valgrind-memcheck**: Detect memory leaks and errors
   ```bash
   make valgrind-memcheck
   ```
   Results saved to `build/valgrind/memcheck.log`

2. **valgrind-massif**: Heap memory profiling
   ```bash
   make valgrind-massif
   ```
   Tracks heap memory usage over time. Analyze with:
   ```bash
   ms_print build/valgrind/massif.out
   ```

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

