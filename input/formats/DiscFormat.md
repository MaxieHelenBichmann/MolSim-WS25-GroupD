\page disc_format Disc Data Format
# Disc Data Format

This document describes the YAML format for defining one disc of particles.

To define a disc, you create a top-level block with the `format` key set to `"Disc"`.

## Disc Object

The Disc object has the following keys:

- `format`: A string describing the format, in this case `"Disc"`.
- `coordinates`: Center coordinates of the disc as `[x, y, z]`.
- `velocity`: Initial velocity of all particles as `[vx, vy, vz]`.
- `radius`: The radius of the disc (number of particles along the radius)
- `mass`: The mass of each particle in the disc.
- `distance`: The spacing between adjacent particles in the disc.
- `mean_velo`: The mean velocity used for the Maxwell-Boltzmann distribution of the particles' velocities.
- `epsilon`: The Epsilon value of all particles in this disc.
- `sigma`: The Sigma value of all particles in this disc.
- `targets`: (Optional) List of particle positions `[x, y, z]` (grid indices) that should receive targeted forces (types 3/4).

## Example

```yaml
my_disc:
    format: Disc
    coordinates: [0.0, 0.0, 0.0]
    velocity: [0.0, 0.0, 0.0]
    radius: 5
    mass: 1.0
    distance: 1.1225
    mean_velo: 0.1
    epsilon: 5.0
    sigma: 1.0
    targets:
      - [0, 0, 0]
      - [4, 4, 0]
```
