MolSim Group D
===

Contributors:
- Maxie Helen Bichmann
- Georg Sebastian Eisner
- Henry Jacob Meyran

## Project Structure

```
.
├── cmake
│   └── modules
│       ├── clang-format.cmake
│       ├── clang-tidy.cmake
│       ├── doxygen.cmake
│       └── vtk.cmake
├── CMakeLists.txt
├── Doxyfile
├── include
│   ├── io
│   │   ├── CLIParse.h
│   │   ├── fileReader
│   │   │   └── XVMReader.h
│   │   ├── FileReader.h
│   │   ├── outputWriter
│   │   │   ├── VTKWriter.h
│   │   │   └── XYZWriter.h
│   │   └── OutputWriter.h
│   ├── particles
│   │   ├── container
│   │   │   ├── ContainerRef.h
│   │   │   └── SimpleContainer.h
│   │   ├── ParticleContainer.h
│   │   └── Particle.h
│   ├── physics
│   │   ├── ForceSource.h
│   │   └── GravitationalForce.h
│   └── utils
│       ├── ArrayUtils.h
│       ├── MaxwellBoltzmannDistribution.h
│       ├── Simulation.h
│       └── Vector.h
├── input
│   └── eingabe-sonne.txt
├── README.md
└── src
    ├── io
    │   ├── fileReader
    │   │   └── XVMReader.cpp
    │   └── outputWriter
    │       ├── VTKWriter.cpp
    │       └── XYZWriter.cpp
    ├── MolSim.cpp
    ├── particles
    │   ├── container
    │   │   ├── ContainerRef.cpp
    │   │   └── SimpleContainer.cpp
    │   └── Particle.cpp
    └── physics
        └── GravitationalSource.cpp

19 directories, 32 files
```

## Dependencies
- **Essential:**
    - CMake Version 3.10+
    - Make 4.3+
    - C++ compiler with C++20 support 
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

2. Now you can interactively change the options (These include Doxygen generation clang-tidy integration and VTK output).

3. Build the executable with the Makefile by running the command:
    ```
    make -j $(nproc)
    ```

## Running the project

**How do I run the project?**

0. Make sure to be in the project root.

1. Run the executable with the given input file and optional desired delta T and end T:

    ```
    ./build/MolSim ../input/eingabe-sonne.txt 0.014 1000
    ```
