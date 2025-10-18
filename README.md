MolSim Group D
===

Contributors:
- Maxie Helen Bichmann
- Georg Sebastian Eisner
- Henry Jacob Meyran

## Project Structure

tba

## Dependencies

- VTK
- clang-format
- Doxygen
- Build Tools (CMake, make)


## Building the project

How do you build the project?

0. Make sure to be in the project root.

1. Configure the project with cmake by running these commands:
```
mkdir build && cd build
ccmake ..
```

2. Now you can interactively change the options.

3. Build the executable with the Makefile by running the command:
```
make -j $(nproc)
```

## Running the project

How do I run the project?

0. Make sure to be in the project root.

1. Run the executable with the given input file:

```
./build/MolSim ./input/eingabe-sonne.txt
```
