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

The settings block has the following keys:

- `format`: A string describing the format, which must be `"Settings"`.
- `delta_t`: (Optional) A floating-point number specifying the time step (`Δt`) for the simulation.
- `end_time`: (Optional) A floating-point number specifying the end time of the simulation.
- `start_time`: (Optional) A floating-point number specifying the start time of the simulation.

## Example

```yaml
settings:
    format: Settings
    delta_t: 0.005
    end_time: 500.0
    start_time: 0.0
```
