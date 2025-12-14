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

Optional keys: 
- `epsilon`: The Epsilon value of the particle.
- `sigma`: The Sigma value of the particle.

## Particle in a Checkpoint File

When creating a checkpoint file, the whole state of each particle has to be defined and has the following **additional** keys:

- `old_coordinates`: A map defining the previous coordinates of the particle.
  - `ox`: The old x-coordinate of the particle.
  - `oy`: The old y-coordinate of the particle.
  - `oz`: The old z-coordinate of the particle.
- `force`: A map defining the force acting on a particle.
  - `fx`: The x-component of the current force.
  - `fy`: The y-component of the current force.
  - `fz`: The z-component of the current force.
- `old_force`: A map defining the previous force acting on a particle.
  - `ofx`: The x-component of the old force.
  - `ofy`: The y-component of the old force.
  - `ofz`: The z-component of the old force.
- `type`: The type of the particle. 

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
