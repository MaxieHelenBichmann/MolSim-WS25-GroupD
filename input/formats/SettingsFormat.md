\page settings_format Settings Data Format
# Settings Data Format

This document describes the YAML format for defining simulation settings.

## General Structure

A settings block is a top-level entry in the YAML file, typically named `settings`.
The settings block **must be the first block** in the YAML file.

```yaml
settings:
  format: Settings
  # ... other settings
```

## Keys

The settings block has the following keys, with every key except format being optional:

- `format`: A string describing the format, which must be `"Settings"`.
- `delta_t`: A floating-point number specifying the time step (`Δt`) for the simulation.
- `end_time`: A floating-point number specifying the end time of the simulation.
- `start_time`: A floating-point number specifying the start time of the simulation.
- `base_name`: A string specifying the base name of the output files of the simulation.
- `force`: A string specifying the force type (either "Lennard Jones" or "Gravitational") of the simulation.
- `frequency`: A whole number specifying the frequency of the output files of the simulation.
- `cutoff`: A floating-point number specifying the cutoff radius of the Linked Cells of the simulation.
- `domain`: A map defining the domain of the simulation (see Domain Configuration below).

## Domain Configuration

The `domain` key contains the following sub-keys:

- `domain_type`: Container type, either `"SIMPLE"` or `"LINKED"` (default: `"LINKED"`).
- `x`: The x-dimension of the domain (default: 1.0).
- `y`: The y-dimension of the domain (default: 1.0).
- `z`: The z-dimension of the domain (default: 1.0).
- `boundaries`: A map defining boundary conditions for each of the 6 domain faces.

### Boundary Conditions

The `boundaries` key can contain up to 6 boundary specifications:

- `left`: Boundary at x = 0 (LEFT, -x direction)
- `right`: Boundary at x = domain_x (RIGHT, +x direction)
- `front`: Boundary at y = 0 (FRONT, -y direction)
- `back`: Boundary at y = domain_y (BACK, +y direction)
- `lower`: Boundary at z = 0 (LOWER, -z direction)
- `upper`: Boundary at z = domain_z (UPPER, +z direction)

Each boundary can have the following keys:

- `type`: The boundary type, either `"OUTFLOW"` or `"REFLECTING"` (default: `"OUTFLOW"`).
- `sigma`: (REFLECTING only) Optional sigma value for ghost particle interactions.
- `epsilon`: (REFLECTING only) Optional epsilon value for ghost particle interactions.

If a boundary is not specified, it defaults to `OUTFLOW`.

## Example

### Minimal Example
```yaml
settings:
    format: Settings
    delta_t: 0.005
    end_time: 500.0
    start_time: 0.0
    base_name: "MD"
    force: "Lennard Jones"
    frequency: 10
    cutoff: 1.
    domain:
      x: 1.
      y: 1.
      z: 1.
```

### Full Example with Boundary Conditions
```yaml
settings:
    format: Settings
    delta_t: 0.0005
    end_time: 20.0
    start_time: 0.0
    base_name: "collision"
    force: "Lennard Jones"
    frequency: 10
    cutoff: 3.0
    container_type: LINKED
    domain:

      x: 180.0
      y: 90.0
      z: 1.0
      boundaries:
        left:
          type: REFLECTING
        right:
          type: REFLECTING
        front:
          type: REFLECTING
        back:
          type: REFLECTING
        upper:
          type: OUTFLOW
        lower:
          type: OUTFLOW
```

### Reflecting Boundary with Custom Parameters
```yaml
settings:
    format: Settings
    delta_t: 0.001
    end_time: 100.0
    domain:
      domain_type: LINKED
      x: 50.0
      y: 50.0
      z: 50.0
      boundaries:
        left:
          type: REFLECTING
          sigma: 1.2
          epsilon: 5.0
        right:
          type: REFLECTING
          sigma: 1.2
          epsilon: 5.0
        front:
          type: OUTFLOW
        back:
          type: OUTFLOW
        upper:
          type: REFLECTING
        lower:
          type: REFLECTING
```
