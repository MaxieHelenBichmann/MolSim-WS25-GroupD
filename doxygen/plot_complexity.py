#!/usr/bin/env python3
"""
Plot complexity benchmark data from Google Benchmark output.
Extracts mean aggregates and creates a visualization comparing
LinkedCell O(n) vs DirectSum O(n²) scaling.

Shows runtime per simulation iteration (delta_t step) for better comparison.
Benchmark runs: end_time=1.0, delta_t=0.0005 => 2000 iterations per run.

Usage: python plot_complexity.py [input_file] [output_file]
  Default input:  benchmark_data/complexityData.txt
  Default output: complexity_plot.png
"""

import re
import sys
import matplotlib.pyplot as plt

# Simulation parameters from FullComplexityBenchmarks.cpp
ITERATIONS_PER_RUN = 2000  # end_time=1.0, delta_t=0.0005


def parse_benchmark_data(filepath):
    """Parse benchmark output and extract mean times."""
    with open(filepath, 'r') as f:
        data = f.read()

    # Extract mean values: "Name/particles/repeats:5_mean    TIME ms"
    pattern = r'Simulation/Complexity/(\w+)/(\d+)/repeats:5_mean\s+([\d.]+) ms'
    matches = re.findall(pattern, data)

    linked_cell = {'particles': [], 'time': [], 'time_per_iter': []}
    direct_sum = {'particles': [], 'time': [], 'time_per_iter': []}

    for container, particles, time in matches:
        total_time = float(time)
        time_per_iter = total_time / ITERATIONS_PER_RUN
        if container == 'LinkedCell':
            linked_cell['particles'].append(int(particles))
            linked_cell['time'].append(total_time)
            linked_cell['time_per_iter'].append(time_per_iter)
        else:
            direct_sum['particles'].append(int(particles))
            direct_sum['time'].append(total_time)
            direct_sum['time_per_iter'].append(time_per_iter)

    return linked_cell, direct_sum


def plot_complexity(linked_cell, direct_sum, output_file):
    """Create and save the complexity comparison plot (runtime per iteration)."""
    plt.figure(figsize=(10, 6))

    if linked_cell['particles']:
        plt.plot(linked_cell['particles'], linked_cell['time_per_iter'], 'o-',
                 label='LinkedCell O(n)', linewidth=2, markersize=8)

    if direct_sum['particles']:
        plt.plot(direct_sum['particles'], direct_sum['time_per_iter'], 's-',
                 label='DirectSum O(n²)', linewidth=2, markersize=8)

    plt.xlabel('Particle Count', fontsize=12)
    plt.ylabel('Time per Iteration (ms)', fontsize=12)
    plt.title('Simulation Complexity Comparison (per Iteration)', fontsize=14)
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"Plot saved to {output_file}")


def main():
    input_file = sys.argv[1] if len(sys.argv) > 1 else '../benchmark_data/complexityData.txt'
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'graphics/complexity_plot.png'

    linked_cell, direct_sum = parse_benchmark_data(input_file)

    print(f"LinkedCell data points: {len(linked_cell['particles'])}")
    print(f"DirectSum data points: {len(direct_sum['particles'])}")
    print(f"Iterations per run: {ITERATIONS_PER_RUN}")

    if linked_cell['particles']:
        print("\nLinkedCell:")
        for p, t, t_iter in zip(linked_cell['particles'], linked_cell['time'], linked_cell['time_per_iter']):
            print(f"  {p:>5} particles: {t:>8.0f} ms total, {t_iter:>6.2f} ms/iter")

    if direct_sum['particles']:
        print("\nDirectSum:")
        for p, t, t_iter in zip(direct_sum['particles'], direct_sum['time'], direct_sum['time_per_iter']):
            print(f"  {p:>5} particles: {t:>8.0f} ms total, {t_iter:>6.2f} ms/iter")

    plot_complexity(linked_cell, direct_sum, output_file)


if __name__ == '__main__':
    main()
