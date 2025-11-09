MolSim Group D
===

Contributors:

- Maxie Helen Bichmann
- Georg Sebastian Eisner
- Henry Jacob Meyran

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
    ./build/MolSim ./input/eingabe-sonne.txt -d 0.014 -t 1000
    ```

## Optional Tools

**Doxygen:**  
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
