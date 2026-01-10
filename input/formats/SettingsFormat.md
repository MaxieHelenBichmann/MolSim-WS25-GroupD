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
- `force`: A string specifying the force type, either "Lennard Jones", "Smooth Lennard Jones" or "Gravitational" (default:`"Lennard Jones"`), of the simulation.
- `frequency`: A whole number specifying the frequency of the output files of the simulation.
- `checkpoint`: A whole number specifying the frequency of the checkpoint files of the simulation.
- `cutoff`: A floating-point number specifying the cutoff radius of the Linked Cells of the simulation.
- `smooth`: A floating-point number specifying the smoothing radius of the Smooth Lennard-Jones Force.
- `container`: Container type, either `"SIMPLE"` or `"LINKED"` (default: `"LINKED"`).
- `domain`: A map defining the domain of the simulation (see Domain Configuration below).
- `thermostat`: A map defining the thermostat of the simulation
- `statistics`: A map defining the collection of thermodynamic statistics

## Thermostat Configuration
The thermostat is only activated when the `thermostat` key is provided. It can contain the following subkeys:  
- `initial_temp`: A floating-point number specifying the initial target temperature of the simulation.
- `target_temp`: A floating-point number specifying the normal target temperature of the simulation.
- `n_thermostat`: A whole point number specifying the frequency with which the thermostat gets applied.
- `delta_temp`: A floating-point number specifying the maximum allowed change of temperature with one application of the thermostat.

## Statistics Configuration
The data collection is only activated when the `statistics` key is provided. It can contain the following subkeys:
- `diffusion`:  A whole number specifying the frequency of the data collection regarding diffusion. If key not provided, data will not be collected.
- `rdf`: A whole number specifying the frequency of the data collection regarding RDF. If key not provided, data will not be collected.
- `sample_r`: A floating-point number specifying the sample width of the RDF (default: 1.0).


## Domain Configuration

The `domain` key contains the following sub-keys:
- `x`: The x-dimension of the domain (default: 1.0).
- `y`: The y-dimension of the domain (default: 1.0).
- `z`: The z-dimension of the domain (default: 1.0).
- `g_grav`: The gravitational force present in the simulation domain (default: 0.0).
- `dimensions`: The number of dimensions the domain has. Supported values: 2, 3.
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

- `type`: The boundary type, either `"OUTFLOW"` or `"REFLECTING"` or `"VELOCITYREFLECT"` or `"PERIODIC"` (default: `"OUTFLOW"`).
- `sigma`: (REFLECTING only) Optional sigma value for ghost particle interactions.
- `epsilon`: (REFLECTING only) Optional epsilon value for ghost particle interactions.
- `ghost_on_boundary`: (REFLECTING only) Optional boolean to place ghost particles on the boundary (default: `false`).

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
    container: LINKED
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
    container: LINKED
    domain:
      x: 50.0
      y: 50.0
      z: 50.0
      g_grav: -9.81
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
