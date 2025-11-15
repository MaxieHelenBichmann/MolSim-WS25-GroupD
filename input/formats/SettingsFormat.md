\page settings_format Settings Data Format
# Settings Data Format

This document describes the YAML format for defining simulation settings.

## General Structure

A settings block is a top-level entry in the YAML file, typically named `settings`.

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
- `domain`: A map defining the domain of the simulation.
  - `x`: The x-coordinate of the corner of the domain.
  - `y`: The y-coordinate of the corner of the domain.
  - `z`: The z-coordinate of the corner of the domain.





## Example

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
