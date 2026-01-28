#!/usr/bin/env python3
"""
Plot thread scaling benchmark results from Google Benchmark JSON output.
Generates plots for Strong Scaling, Weak Scaling, and Force Calculation benchmarks.

Usage:
    python3 scripts/statistics/plot_thread_scaling.py <input_json_file> [output_prefix]
"""

import json
import sys
import matplotlib.pyplot as plt
import re

PHYSICAL_CORES = 112

def parse_json(filepath):
    with open(filepath, 'r') as f:
        data = json.load(f)
    return data['benchmarks']

def extract_data(benchmarks):
    data = {
        'Strong': {'Naive': {'x': [], 'y': []}, 'Coloring': {'x': [], 'y': []}},
        'Weak': {'Naive': {'x': [], 'y': []}, 'Coloring': {'x': [], 'y': []}},
        'ForceCalc': {'Naive': {'x': [], 'y': []}, 'Coloring': {'x': [], 'y': []}},
    }

    # Regex to parse benchmark names
    # Example: ThreadScaling/Strong/100kParticles/Naive/1/0
    # Capture groups: Type (Strong/Weak/ForceCalc), Strategy (Naive/Coloring), Threads (arg0)
    # Note: Google Benchmark might format args differently depending on version, 
    # but usually "Name/Arg0/Arg1" or "Name/threads:X"
    
    for bm in benchmarks:
        name = bm['name']
        
        # Skip aggregate rows if present (ends with _mean, _median, etc.)
        if name.endswith('_mean') or name.endswith('_median') or name.endswith('_stddev'):
            # Unless we WANT aggregates. The user's config used ReportAggregatesOnly(true)
            # which produces _mean, _median, _stddev. We should use _mean.
            pass
        else:
            # If aggregates only is on, strict measurements might be missing, 
            # but usually google benchmark outputs the aggregates with special suffixes
            # Let's check if it's an aggregate
            continue

        # Check for our specific benchmarks
        # We expect names like "ThreadScaling/Strong/.../Naive/threads:X/..." or similar
        # Since ReportAggregatesOnly is true, we look for names ending in _mean
        pass 

    # Re-loop handling aggregates
    for bm in benchmarks:
        name = bm['name']
        if not name.endswith('_mean'):
            continue
            
        # Remove suffix
        base_name = name.replace('_mean', '')
        
        # Determine Type
        b_type = None
        if 'Strong' in base_name:
            b_type = 'Strong'
        elif 'Weak' in base_name:
            b_type = 'Weak'
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
        # The args are usually at the end. 
        # For our C++ code: ->Args({threads, strategy_enum})
        # The name usually looks like "BaseName/threads/strategy_enum"
        # Let's try to extract integer segments
        parts = base_name.split('/')
        
        threads = -1
        # Traverse parts in reverse to find the thread count
        # Args are typically appended at the end
        for part in reversed(parts):
            if part.isdigit():
                val = int(part)
                # Threads are usually the first arg, strategy enum (0 or 1) is second
                # But threads are large (1, 2, 4...) and strategy is 0/1.
                # However, strategy is also encoded in the text "Naive"/"Coloring".
                # The thread count is likely the one that varies widely.
                if val >= 1: 
                    # If we found a plausible thread count (1 can be ambiguous with enum, 
                    # but typically thread count comes before strategy arg in Args({t, s}))
                    # Actually, Args({t, s}) -> Name/t/s
                    # So 's' is the last one (0 or 1), 't' is the second to last.
                    pass
        
        # Safer approach: Rely on the split structure matching the registration
        # ThreadScaling/Strong/100kParticles/Naive/1/0_mean
        # parts: [ThreadScaling, Strong, 100kParticles, Naive, 1, 0]
        # threads is parts[-2]
        
        try:
            # We assume the last two parts are threads and strategy integer
            # Example: .../Naive/112/0
            threads = int(parts[-2])
        except (ValueError, IndexError):
            # Fallback for manual check or different formatting
            # Maybe formatted as "threads:1"
            continue

        # Extract Value
        # Strong/Weak use UpdatesPerSec (custom counter)
        # ForceCalc uses real_time (time per iter)
        
        y_val = 0.0
        if b_type == 'ForceCalc':
            # Use real_time (microseconds)
            y_val = bm['real_time']
        else:
            # Use UpdatesPerSec if available
            if 'UpdatesPerSec' in bm:
                y_val = bm['UpdatesPerSec']
            else:
                # Fallback to inverse of real_time if metric missing?
                # But let's assume it's there as per C++ code
                print(f"Warning: UpdatesPerSec missing for {name}, skipping.")
                continue

        data[b_type][strategy]['x'].append(threads)
        data[b_type][strategy]['y'].append(y_val)

    # Sort data by threads
    for b_type in data:
        for strat in data[b_type]:
            combined = sorted(zip(data[b_type][strat]['x'], data[b_type][strat]['y']))
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
        
        # Mark Physical Cores
        ax.axvline(x=PHYSICAL_CORES, color='green', linestyle='--', label=f'Physical Cores ({PHYSICAL_CORES})')
        # Add a text label or legend for the vertical line
        # ax.text(PHYSICAL_CORES, ax.get_ylim()[0], ' Physical Limit', rotation=90, verticalalignment='bottom')
        
        if is_log_x:
            ax.set_xscale('log', base=2)
            # Custom ticks to show 112 clearly
            ticks = sorted(list(set(data_dict['Naive']['x'] + data_dict['Coloring']['x'] + [PHYSICAL_CORES])))
            # Filter ticks to avoid overcrowding if needed, but 10 ticks is fine
            ax.set_xticks(ticks)
            ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
            
        if is_log_y:
            ax.set_yscale('log')

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 plot_thread_scaling.py <input_json> [output_prefix]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_prefix = sys.argv[2] if len(sys.argv) > 2 else "thread_scaling"

    try:
        raw_benchmarks = parse_json(input_file)
    except FileNotFoundError:
        print(f"Error: File {input_file} not found.")
        sys.exit(1)
    except json.JSONDecodeError:
        print(f"Error: Failed to parse JSON from {input_file}.")
        sys.exit(1)

    data = extract_data(raw_benchmarks)

    # 1. Strong Scaling
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_benchmark(ax, 'Strong Scaling (Fixed Problem Size)', 
                   data['Strong'], 'Updates / Sec (Higher is Better)')
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_strong.png")
    print(f"Strong scaling plot saved to {output_prefix}_strong.png")
    plt.close()

    # 2. Weak Scaling
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_benchmark(ax, 'Weak Scaling (Scaled Problem Size)', 
                   data['Weak'], 'Updates / Sec (Higher is Better)')
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_weak.png")
    print(f"Weak scaling plot saved to {output_prefix}_weak.png")
    plt.close()

    # 3. Force Calculation
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_benchmark(ax, 'Force Calculation Only', 
                   data['ForceCalc'], 'Time (us) (Lower is Better)', is_log_y=True)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_force.png")
    print(f"Force calculation plot saved to {output_prefix}_force.png")
    plt.close()

if __name__ == '__main__':
    main()
