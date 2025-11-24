\page cuboid_format Cuboid Data Format
# Cuboid Data Format

This document describes the YAML format for defining one cuboid of particles.

To define a cuboid, you create an object with the `format` key set to `"Cuboid"`.

## Cuboid Object

The cuboid block has the following keys:

- `format`: A string describing the format, in this case `"Cuboid"`.
- `coordinates`: A map defining the starting coordinates of the cuboid.
  - `x`: The x-coordinate of the corner of the cuboid.
  - `y`: The y-coordinate of the corner of the cuboid.
  - `z`: The z-coordinate of the corner of the cuboid.
- `velocity`: A map defining the initial velocity of all particles in the cuboid.
  - `vx`: The x-component of the velocity.
  - `vy`: The y-component of the velocity.
  - `vz`: The z-component of the velocity.
- `particleNum`: A map defining the number of particles along each dimension of the cuboid.
  - `nx`: The number of particles in the x-direction.
  - `ny`: The number of particles in the y-direction.
  - `nz`: The number of particles in the z-direction.
- `mass`: The mass of each particle in the cuboid.
- `distance`: The spacing between adjacent particles in the cuboid.
- `mean_velo`: The mean velocity used for the Maxwell-Boltzmann distribution of the particles' velocities.
- `epsilon`: The Epsilon value of all particles in this cuboid
- `sigma`: The Sigma value of all particles in this cuboid
## Example

```yaml
my_cuboid:
    format: Cuboid
    coordinates:
        x: 0.0
        y: 0.0
        z: 0.0
    velocity:
        vx: 0.0
        vy: 0.0
        vz: 0.0
    particleNum:
        nx: 40
        ny: 8
        nz: 1
    mass: 1.0
    distance: 1.1225
    mean_velo: 0.1
    epsilon: 5.0
    sigma: 1.0
```
