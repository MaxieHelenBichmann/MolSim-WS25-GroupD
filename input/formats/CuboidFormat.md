\page cuboid_format Cuboid Data Format
# Cuboid Data Format

This document describes the YAML format for defining one cuboid of particles.

To define a cuboid, you create an object with the `format` key set to `"Cuboid"`.

## Cuboid Object

The cuboid block has the following keys:

- `format`: A string describing the format, in this case `"Cuboid"`.
- `coordinates`: Position of the corner of the cuboid as `[x, y, z]`.
- `velocity`: Initial velocity of all particles as `[vx, vy, vz]`.
- `particleNum`: Number of particles along each dimension as `[nx, ny, nz]`.
- `mass`: The mass of each particle in the cuboid.
- `distance`: The spacing between adjacent particles in the cuboid.
- `mean_velo`: The mean velocity used for the Maxwell-Boltzmann distribution of the particles' velocities.
- `epsilon`: The Epsilon value of all particles in this cuboid.
- `sigma`: The Sigma value of all particles in this cuboid.
- `targets`: (Optional) List of particle positions `[x, y, z]` (grid indices) that should receive targeted forces (types 3/4).

## Example

```yaml
my_cuboid:
    format: Cuboid
    coordinates: [0.0, 0.0, 0.0]
    velocity: [0.0, 0.0, 0.0]
    particleNum: [40, 8, 1]
    mass: 1.0
    distance: 1.1225
    mean_velo: 0.1
    epsilon: 5.0
    sigma: 1.0
    targets:
      - [0, 0, 0]
      - [39, 7, 0]
```
