# Benchmarks {#benchmarks}

This page documents the performance benchmarks and complexity analysis for MolSim.

## Complexity Analysis

The simulation supports two particle container types with different algorithmic complexity:

- **LinkedCellContainer**: O(n) complexity using spatial partitioning
- **SimpleContainer (DirectSum)**: O(n²) complexity using direct pairwise computation

### Complexity Benchmark Results

The following plot shows runtime scaling with increasing particle count:

![Complexity Comparison](complexity_plot.png)
