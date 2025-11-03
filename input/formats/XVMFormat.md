\page xvm_format XVM Data Format
# XVM Data Format

This document describes the YAML format for defining a collection of particles.

To define a set of particles, you create a top-level block with the `format` key set to `"XVM"`.

## XVM Block

The XVM block has the following keys:

- `format`: A string describing the format, in this case `"XVM"`.
- `num_particles`: An integer specifying the number of particle objects defined in this block.
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
my_particles:
    format: XVM
    num_particles: 2
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
```
