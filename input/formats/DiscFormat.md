\page disc_format Disc Data Format
# Disc Data Format

This document describes the YAML format for defining one disc of particles.

To define a disc, you create a top-level block with the `format` key set to `"Disc"`.

## Disc Object

The Disc object has the following keys:

- `format`: A string describing the format, in this case `"Disc"`.
- `coordinates`: A map defining the center coordinates of the disc.
  - `x`: The x-coordinate of the center of the disc.
  - `y`: The y-coordinate of the center of the disc.
  - `z`: The z-coordinate of the center of the disc.
- `velocity`: A map defining the initial velocity of all particles in the disc.
  - `vx`: The x-component of the velocity.
  - `vy`: The y-component of the velocity.
  - `vz`: The z-component of the velocity.
- `radius`: The radius of the disc (number of particles along the radius)
- `mass`: The mass of each particle in the disc.
- `distance`: The spacing between adjacent particles in the disc.
- `mean_velo`: The mean velocity used for the Maxwell-Boltzmann distribution of the particles' velocities.
- `epsilon`: The Epsilon value of all particles in this disc.
- `sigma`: The Sigma value of all particles in this disc.
## Example

```yaml
my_disc:
    format: Disc
    coordinates:
      x: 0.0
      y: 0.0
      z: 0.0
    velocity:
      vx: 0.0
      vy: 0.0
      vz: 0.0
    radius: 5
    mass: 1.0
    distance: 1.1225
    mean_velo: 0.1
    epsilon: 5.0
    sigma: 1.0
```
