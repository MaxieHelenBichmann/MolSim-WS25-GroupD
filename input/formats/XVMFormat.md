\page xvm_format XVM Data Format
# XVM Data Format

This document describes the YAML format for defining a collection of particles.

## General Structure

The YAML file has the following top-level keys:

- `format`: A string describing the format, in this case `"XVM"`.
- `num_particles`: An integer specifying the number of particle objects defined in the file.
- `particles`: A list of particle objects.

## Particle Object

Each object in the `particles` list defines a single particle and has the following keys:

- `coordinates`: A map defining the coordinates of the particle.
  - `x`: The x-coordinate of the particle.
  - `y`: The y-coordinate of the particle.
  - `z`: The z-coordinate of the particle.
- `velocity`: A map defining the initial velocity of the particle.
  - `vx`: The x-component of the velocity.
  - `vy`: The y-component of the velocity.
  - `vz`: The z-component of the velocity.
- `mass`: The mass of the particle.

## Example

```yaml
format: XVM
num_particles: 4
particles:
  - coordinates:
      x: 0.0
      y: 0.0
      z: 0.0
    velocity:
      vx: 0.0
      vy: 0.0
      vz: 0.0
    mass: 1.0

  - coordinates:
      x: 0.0
      y: 1.0
      z: 0.0
    velocity:
      vx: -1.0
      vy: 0.0
      vz: 0.0
    mass: 3.0e-6

  - coordinates:
      x: 0.0
      y: 5.36
      z: 0.0
    velocity:
      vx: -0.425
      vy: 0.0
      vz: 0.0
    mass: 9.55e-4

  - coordinates:
      x: 34.75
      y: 0.0
      z: 0.0
    velocity:
      vx: 0.0
      vy: 0.0296
      vz: 0.0
    mass: 1.0e-14
```
