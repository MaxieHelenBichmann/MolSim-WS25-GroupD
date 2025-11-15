
\page format_main Formats

This section documents the various input file formats.

## General Structure

The YAML input file is structured as a map of top-level blocks. Each block is a named entry that defines a set of particles or global settings. The parser identifies the content of a block by its `format` key.

A file can contain multiple blocks of different formats.

### Example of a file with multiple blocks

```yaml
# A block for global settings
settings:
  format: Settings
  delta_t: 0.014

# A block defining a cuboid of particles
my_cuboid:
  format: Cuboid
  # ... cuboid data

# A block defining a set of individual particles
my_particles:
  format: XVM
  # ... XVM data
```

The following pages describe the specific keys and values for each format type.

\subpage disc_format
\subpage cuboid_format    
\subpage xvm_format   
\subpage settings_format


