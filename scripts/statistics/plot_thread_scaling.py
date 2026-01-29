
#!/usr/bin/env python3
"""
Plot thread scaling benchmark results from Google Benchmark output.
Supports both JSON and console (.out) tabular formats.

Usage:
    python3 scripts/statistics/plot_thread_scaling.py <file1> [file2 ...] [output_prefix]
"""

import json
import sys
import matplotlib.pyplot as plt
import re

PHYSICAL_CORES = 112

def parse_tabular(filepath):
    benchmarks = []
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Skip headers
    start_idx = 0
    for i, line in enumerate(lines):
        if line.startswith('---'):
            start_idx = i + 1
            break
            
    for line in lines[start_idx:]:
        parts = line.split()
        if not parts or not parts[0].startswith('ThreadScaling'):
            continue
            
        # Example: ThreadScaling/Strong/100kParticles/Naive/1/0/repeats:5_mean 359608 ms 359527 ms 5 ...
        name = parts[0]
        try:
            time_val = float(parts[1])
            
            benchmarks.append({
                'name': name,
                'real_time': time_val
            })
        except (ValueError, IndexError):
            continue
    return benchmarks

def parse_and_merge(filepaths):
    all_benchmarks = []
    for filepath in filepaths:
        try:
            if filepath.endswith('.json'):
                with open(filepath, 'r') as f:
                    data = json.load(f)
                    if 'benchmarks' in data:
                        all_benchmarks.extend(data['benchmarks'])
            else:
                # Assume tabular .out format
                all_benchmarks.extend(parse_tabular(filepath))
        except Exception as e:
            print(f"Warning: Could not process {filepath}: {e}")
    return all_benchmarks

def extract_data(benchmarks):
    data = {
        'Strong_100k': {'Naive': {'x': [], 'y': []}, 'Coloring': {'x': [], 'y': []}},
        'Strong_Contest1': {'Naive': {'x': [], 'y': []}, 'Coloring': {'x': [], 'y': []}},
        'ForceCalc': {'Naive': {'x': [], 'y': []}, 'Coloring': {'x': [], 'y': []}},
    }

    for bm in benchmarks:
        name = bm['name']
        if not name.endswith('_mean'):
            continue
            
        base_name = name.replace('_mean', '')
        
        # Determine Type
        b_type = None
        if 'Strong' in base_name:
            if 'Contest_1' in base_name:
                b_type = 'Strong_Contest1'
            else:
                b_type = 'Strong_100k'
        elif 'ForceCalc' in base_name:
            b_type = 'ForceCalc'
        else:
            continue

        # Determine Strategy
        strategy = None
        if 'Naive' in base_name:
            strategy = 'Naive'
        elif 'Coloring' in base_name:
            strategy = 'Coloring'
        else:
            continue

        # Determine Threads
        parts = base_name.split('/')
        threads = -1
        found_ints = []
        for part in reversed(parts):
            if part.startswith('repeats:'):
                continue
            # Handle cases like "1/0" or "threads:1"
            clean_part = re.sub(r'[^0-9]', '', part)
            if clean_part.isdigit():
                found_ints.append(int(clean_part))
        
        if len(found_ints) >= 2:
            threads = found_ints[1]
        else:
            continue

        y_val = bm.get('real_time', 0.0)
        data[b_type][strategy]['x'].append(threads)
        data[b_type][strategy]['y'].append(y_val)

    # Sort and remove duplicates (in case .json and .out are both provided)
    for b_type in data:
        for strat in data[b_type]:
            # Use dict to deduplicate by x, preferring non-zero y values
            unique_points = {}
            for x, y in zip(data[b_type][strat]['x'], data[b_type][strat]['y']):
                if x not in unique_points:
                    unique_points[x] = y
                else:
                    # Prefer non-zero value if existing is zero
                    if unique_points[x] == 0.0 and y != 0.0:
                        unique_points[x] = y
            
            combined = sorted(unique_points.items())
            if combined:
                data[b_type][strat]['x'], data[b_type][strat]['y'] = zip(*combined)
            else:
                data[b_type][strat]['x'], data[b_type][strat]['y'] = [], []
                
    return data

def plot_benchmark(ax, title, data_dict, y_label, is_log_x=True, is_log_y=False):
    ax.set_title(title)
    ax.set_xlabel('Threads')
    ax.set_ylabel(y_label)
    
    has_data = False
    for strategy, color in [('Naive', 'red'), ('Coloring', 'blue')]:
        if data_dict[strategy]['x']:
            ax.plot(data_dict[strategy]['x'], data_dict[strategy]['y'], 
                    marker='o', label=strategy, color=color)
            has_data = True
    
    if has_data:
        ax.legend()
        ax.grid(True, which="both", ls="-", alpha=0.5)
        ax.axvline(x=PHYSICAL_CORES, color='green', linestyle='--', label=f'Physical Cores ({PHYSICAL_CORES})')
        # Add text label for physical cores, slightly shifted right
        ax.text(PHYSICAL_CORES * 1.05, 0.95, f'{PHYSICAL_CORES} Cores', 
                transform=ax.get_xaxis_transform(), 
                color='green', ha='left', va='top', rotation=90)
        
        if is_log_x:
            ax.set_xscale('log', base=2)
            ticks = [1, 2, 4, 8, 16, 32, 64, 128, 256]
            ax.set_xticks(ticks)
            ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
            
        if is_log_y:
            ax.set_yscale('log')

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 plot_thread_scaling.py <file1> [file2 ...] [output_prefix]")
        sys.exit(1)

    args = sys.argv[1:]
    output_prefix = "thread_scaling"
    input_files = []

    # Identify output prefix (last arg if it doesn't have a known extension)
    if not (args[-1].endswith('.json') or args[-1].endswith('.out')):
        if len(args) > 1:
            output_prefix = args[-1]
            input_files = args[:-1]
        else:
            input_files = args
    else:
        input_files = args

    raw_benchmarks = parse_and_merge(input_files)
    if not raw_benchmarks:
        print("Error: No benchmark data found.")
        sys.exit(1)

    data = extract_data(raw_benchmarks)

    # 1. Strong Scaling - 100k
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_benchmark(ax, 'Strong Scaling (100k Particles)', 
                   data['Strong_100k'], 'Time (ms) (Lower is Better)')
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_strong_100k.png")
    print(f"Strong scaling (100k) plot saved to {output_prefix}_strong_100k.png")
    plt.close()

    # 2. Strong Scaling - Contest 1
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_benchmark(ax, 'Strong Scaling (Contest 1)', 
                   data['Strong_Contest1'], 'Time (ms) (Lower is Better)')
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_strong_contest1.png")
    print(f"Strong scaling (Contest 1) plot saved to {output_prefix}_strong_contest1.png")
    plt.close()

    # 3. Force Calculation
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_benchmark(ax, 'Force Calculation Only', 
                   data['ForceCalc'], 'Time (us) (Lower is Better)')
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_force.png")
    print(f"Force calculation plot saved to {output_prefix}_force.png")
    plt.close()

if __name__ == '__main__':
    main()
