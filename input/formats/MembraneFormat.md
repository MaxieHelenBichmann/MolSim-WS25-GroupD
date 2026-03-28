\page membrane_format Membrane Data Format
# Membrane Data Format

This document describes the YAML format for defining one 2D membrane of particles.

To define a membrane, you create an object with the `format` key set to `"Membrane"`.

## Membrane Object

The membrane block has the following keys:

- `format`: A string describing the format, in this case `"Membrane"`.
- `coordinates`: Position of the corner of the membrane as `[x, y, z]`.
- `velocity`: Initial velocity of all particles as `[vx, vy, vz]`.
- `particleNum`: Number of particles along each dimension as `[nx, ny]`.
- `mass`: The mass of each particle in the membrane.
- `distance`: The spacing between adjacent particles in the membrane.
- `mean_velo`: The mean velocity used for the Maxwell-Boltzmann distribution of the particles' velocities.
- `epsilon`: The Epsilon value of all particles in this membrane.
- `sigma`: The Sigma value of all particles in this membrane.
- `targets`: (Optional) List of particle positions `[x, y]` (grid indices) that should receive targeted forces (types 3/4).

## Example

```yaml
my_membrane:
    format: Membrane
    coordinates: [0.0, 0.0, 0.0]
    velocity: [0.0, 0.0, 0.0]
    particleNum: [40, 8]
    mass: 1.0
    distance: 1.1225
    mean_velo: 0.1
    epsilon: 5.0
    sigma: 1.0
    targets:
      - [0, 0]
      - [39, 7]
```
