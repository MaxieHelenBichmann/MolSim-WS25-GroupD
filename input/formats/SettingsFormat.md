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
- `frequency`: A whole number specifying the frequency of the output files of the simulation.
- `checkpoint`: A whole number specifying the frequency of the checkpoint files of the simulation.
- `cutoff`: A floating-point number specifying the cutoff radius of the Linked Cells of the simulation.
- `smooth`: A floating-point number specifying the smoothing radius of the Smooth Lennard-Jones Force.
- `container`: Container type, either `"SIMPLE"` or `"LINKED"` (default: `"LINKED"`).
- `pairwise_forces`: A sequence containing all desired pairwise forces (forces applied between particle pairs).
  - Available forces: `"GRAVITATIONAL"`, `"LENNARDJONES"`, `"TRUNCLENNARDJONES"`, `"SMOOTHLENNARDJONES"`
  - Multiple forces can be specified and will all be applied
  - Default: `["LENNARDJONES"]`
- `single_forces`: A sequence of force objects defining single particle forces (forces applied to individual particles).
  - Each force object must have a `type` field and force-specific parameters
  - Multiple forces can be specified and will all be applied
  - Available force types:
    - `GRAV`: Uniform gravitational pull
    - `HARMONIC`: Harmonic bonds for membrane structures
  - See Single Force Configuration below for details
- `target_force`: Configuration for targeted force applied to specific particles (see Target Force Configuration below).
- `domain`: A map defining the domain of the simulation (see Domain Configuration below).
- `thermostat`: A map defining the thermostat of the simulation
- `statistics`: A map defining the collection of thermodynamic statistics

## Force Types

### Pairwise Forces
Pairwise forces are applied between pairs of particles within the cutoff radius. Multiple pairwise forces can be active simultaneously.

- **`GRAVITATIONAL`**: Newtonian gravitational attraction between particle pairs. Force magnitude: `F = G * m1 * m2 / r²`
- **`LENNARDJONES`**: Standard Lennard-Jones potential for molecular interactions. Includes both attractive and repulsive components.
- **`TRUNCLENNARDJONES`**: Truncated Lennard-Jones potential cut off at `r = 2^(1/6) * σ`. Only the repulsive part is active. Used for preventing self-penetration in membrane simulations. Only applies to particles with type 2.
- **`SMOOTHLENNARDJONES`**: Smooth Lennard-Jones potential that smoothly transitions to zero at the cutoff radius. Requires `cutoff` and `smooth` parameters. The smoothing radius (`smooth`) must be less than or equal to the cutoff radius (`cutoff`).

### Single Force Configuration

Single forces are applied to individual particles. Each force requires specific parameters:

#### GRAV - Gravitational Pull
Applies uniform gravitational acceleration to all particles.

**Required Parameters:**
- `g_grav`: Gravitational acceleration vector as `[gx, gy, gz]` (default: `[0.0, -9.81, 0.0]`)

**Example:**
```yaml
single_forces:
  - type: GRAV
    g_grav: [0.0, -9.81, 0.0]  # Gravity in -y direction
```

#### HARMONIC - Membrane Bonds
Applies harmonic spring forces between neighboring particles in membrane structures. Only affects particles with type 2 or 4 that have neighbors set.

**Required Parameters:**
- `k`: Stiffness constant (spring constant) in N/m
- `r_0`: Average bond length (equilibrium distance) in meters

**Force Formula:**
- Direct neighbors (up/down/left/right): `F = k/2 * (|r| - r_0) * (r / |r|)`
- Diagonal neighbors: `F = k/2 * (|r| - √2*r_0) * (r / |r|)`

**Example:**
```yaml
single_forces:
  - type: HARMONIC
    k: 300.0      # Stiffness constant
    r_0: 1.2      # Equilibrium bond length
```

## Target Force Configuration

The `target_force` enables forces applied only to specific particles (types 3 and 4). Target particles must be specified in the generator configuration using the `targets` field.

**Parameters:**
- `direction`: Force direction vector as `[dx, dy, dz]`
- `magnitude`: Force magnitude (scalar)
- `max_iterations`: Apply force only for the first N iterations

**Example:**
```yaml
settings:
  target_force:
    direction: [0.0, 1.0, 0.0]  # Pull upward
    magnitude: 10.0
    max_iterations: 1000        # Apply for first 1000 steps
```
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
- `window_size`: A floating-point number specifying the size of the considered window for the RDF (default: 10.0).


## Domain Configuration

The `domain` key contains the following sub-keys:
- `coordinates`: coordinates of the left lower corner of the domain (format [x, y, z]).
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
    frequency: 10
    cutoff: 1.0
    pairwise_forces: [LENNARDJONES]
    single_forces:
      - type: GRAV
        g_grav: [0.0, -9.81, 0.0]
    domain:
      x: 1.0
      y: 1.0
      z: 1.0
```

### Multiple Forces Example
```yaml
settings:
    format: Settings
    delta_t: 0.0005
    end_time: 20.0
    pairwise_forces:
      - LENNARDJONES      # Regular LJ between all particles
      - TRUNCLENNARDJONES # Repulsive-only LJ for membrane particles
    single_forces:
      - type: GRAV
        g_grav: [0.0, -0.981, 0.0]  # Gravity in -y direction
      - type: HARMONIC
        k: 300.0
        r_0: 1.2
    domain:
      x: 50.0
      y: 50.0
      z: 1.0
```

### Membrane Simulation Example
```yaml
settings:
    format: Settings
    delta_t: 0.0001
    end_time: 10.0
    container: LINKED
    cutoff: 3.0
    pairwise_forces:
      - LENNARDJONES
      - TRUNCLENNARDJONES
    single_forces:
      - type: HARMONIC
        k: 300.0        # Spring stiffness
        r_0: 1.1225     # Equilibrium bond length
    target_force:
      direction: [0.0, 1.0, 0.0]
      magnitude: 5.0
      max_iterations: 500
    domain:
      x: 40.0
      y: 40.0
      z: 1.0
```

### Smooth Lennard-Jones Example
```yaml
settings:
    format: Settings
    delta_t: 0.0005
    end_time: 20.0
    container: LINKED
    cutoff: 3.0
    smooth: 2.5       # Smoothing radius (must be <= cutoff)
    pairwise_forces:
      - SMOOTHLENNARDJONES
    domain:
      x: 50.0
      y: 50.0
      z: 1.0
```

### Backward Compatible Example (Deprecated)
```yaml
settings:
    format: Settings
    delta_t: 0.005
    end_time: 500.0
    force: "Lennard Jones"  # Deprecated: converted to pairwise_forces: ["LENNARDJONES"]
    frequency: 10
    domain:
      x: 1.0
      y: 1.0
```

### Full Example with Boundary Conditions
```yaml
settings:
    format: Settings
    delta_t: 0.0005
    end_time: 20.0
    start_time: 0.0
    base_name: "collision"
    frequency: 10
    cutoff: 3.0
    container: LINKED
    pairwise_forces: [LENNARDJONES]
    single_forces:
      - type: GRAV
        g_grav: [0.0, -9.81, 0.0]
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
