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

## Building the project

**How do you build the project?**

0. Make sure to be in the project root.

1. Configure the project with cmake by running these commands:

    ```
    mkdir build && cd build
    ccmake ..
    ```
    <details><summary>Configure Options</summary>
    - <b>ENABLE_TESTING</b> enables all tests.<br />
    - <b>ENABLE_BENCHMARK</b> enables benchmarking. (<b>BENCHMARK_DOWNLOAD_DEPENDENCIES</b> should also be enabled)<br />
    - <b>ENABLE_DOXYGEN</b> enables doxygen.<br />
    - <b>COVERAGE</b> enables code coverage reports.<br />
    - <b>ENABLE_VTK_OUTPUT</b> enables output in the vtk format.<br />
    - <b>CMAKE_BUILD_TYPE</b> specifies the build type (Debug, Release,...).<br />
   </details>  

 
2. Now you can interactively change the options (These include Doxygen generation, clang-tidy integration and VTK output).

3. Build the executable with the Makefile by running the command:

    ```
    make -j $(nproc)
    ```

## Running the project

**How do I run the project?**

0. Make sure to be in the project root.

1. Run the executable with the given input file and optional desired delta_T ( -d) and end_T ( -t):

    ```
    ./build/project/MolSim ./input/planets.yaml -d 0.014 -t 1000
    ```
    ```
    ./build/project/MolSim ./input/particles.yaml -d 0.0002 -t 5
    ```

## Testing the project

**How do I test the project?**

0. Make sure to be in the project root.

1. Run the test executable or use ctest:
    ```
    ./build/tests/tests
    ```
    ```
    ctest -V --test-dir ./build/tests 
    ```

## Benchmarking the project

**How do I benchmark the project?**

0. Make sure to be in the project root.

1. Run the benchmark executable:
    ```
    ./build/benchmarks/MolSimBench
    ```

## Optional Tools

**Doxygen**  
If Doxygen support was enabled:

```
make doc_doxygen
```

Creates the doxygen documentation in the doxys_documentation directory.

**Clang-tidy**  
If Clang-tidy support was enabled, clang-tidy will automatically run the checks specified in the .clang-tidy file.  
There is also:

```
make fix
```

Which automatically applies fixes where clang-tidy finds them.

**Test Coverage**  
If coverage was enabled:

```
make coverage_report
```

Creates a coverage report on the coverage directory inside the build directory.
