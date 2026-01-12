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

- `coordinates`: Position of the particle as `[x, y, z]`.
- `velocity`: Initial velocity of the particle as `[vx, vy, vz]`.
- `mass`: The mass of the particle. 

Optional keys: 
- `epsilon`: The Epsilon value of the particle.
- `sigma`: The Sigma value of the particle.
- `type`: The type of the particle (default: 0).

## Particle in a Checkpoint File

When creating a checkpoint file, the whole state of each particle has to be defined and has the following **additional** keys:

- `old_coordinates`: Previous coordinates of the particle as `[ox, oy, oz]`.
- `force`: Force acting on the particle as `[fx, fy, fz]`.
- `old_force`: Previous force acting on the particle as `[ofx, ofy, ofz]`.

## Example

```yaml
my_particles:
    format: XVM
    num_particles: 2
    particles:
      - coordinates: [0.0, 0.0, 0.0]
        velocity: [0.0, 0.0, 0.0]
        mass: 1.0

      - coordinates: [0.0, 1.0, 0.0]
        velocity: [-1.0, 0.0, 0.0]
        mass: 3.0e-6
```
