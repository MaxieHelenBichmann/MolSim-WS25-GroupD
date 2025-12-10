import sys
import typing
import json 
import pandas as pd
import matplotlib.pyplot as plt


# Usage: python analyse_data.py path/to/data
# Prints run Information and generates plots for each type of Benchmark that was run

def print_context(context_data: typing.Dict[str, typing.Any]): 
    print(f"Name of the host machine: {context_data["host_name"]}")
    print(f"Number of CPU threads available: {context_data["num_cpus"]}")
    print(f"Frequency of the CPU: {context_data["mhz_per_cpu"]}MHz | CPU Scaling enabled: { context_data["cpu_scaling_enabled"]}")
    print("Cache hierarchy: ")
    for cache in context_data["caches"]:
        print(f"L{cache["level"]} {cache["type"]} Cache | Size: {cache["size"]} Byte")
    print(f"Build Type: {context_data["library_build_type"]}")

def plot_benchmark(benchmark_name : str, time_unit : str, data):
    plt.figure(figsize=(10, 6))
    unique_types = data['Type'].unique()
    for run_type in unique_types:
        subset = data[data['Type'] == run_type]
        
        subset = subset.sort_values('Iterations')
        
        plt.plot(
            subset['Iterations'], 
            subset['RealTime_mean'], 
            marker='o', 
            label=run_type  
        )
        
        plt.fill_between(
            subset['Iterations'],
            subset['RealTime_mean'] - subset['RealTime_stddev'],
            subset['RealTime_mean'] + subset['RealTime_stddev'], 
            alpha=0.15,
            color='grey'
        )

    plt.title(f"Benchmark: {benchmark_name}")
    plt.xlabel("Input Size")
    plt.ylabel(f"Time ({time_unit})")
    plt.legend(title="Implementation")
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    
    # Optional: Use Log scale if sizes vary wildly (e.g., 32 to 1,000,000)
    # plt.xscale('log')
    # plt.yscale('log')
    
    plt.tight_layout()
    plt.savefig(f"Benchmark_{benchmark_name.split(" ")[0]}_{benchmark_name.split(" ")[1]}.png")
    plt.close()


def main(args: list[str]):
    filename = args[1]
    with open(file=filename, mode="r", encoding="utf-8") as file:
        file_data = json.load(file)
        print_context(file_data["context"])
        aggregate_data = [b for b in file_data["benchmarks"] if b["run_type"] == "aggregate" and b.get("aggregate_name") not in ["BigO", "RMS"]]
        class_name = ""
        processed_data = []
        for bench in aggregate_data: 
            name_part = bench["run_name"].split(sep='/')
            if class_name == "":
                class_name = name_part[0]
            benchmark_name = name_part[1]
            benchmark_type = name_part[2]
            iterations = int(name_part[3])
            stat_type = bench["aggregate_name"]

            processed_data.append({
                'Benchmark': benchmark_name,
                'Type': benchmark_type,
                'Iterations': iterations,
                'Stat': stat_type,
                'RealTime': bench['real_time'],
                'CpuTime': bench['cpu_time'],
                'TimeUnit': bench['time_unit']
            })
        data_frame = pd.DataFrame(processed_data)

        data_frame_pivot = data_frame.pivot_table(
            index=['Benchmark', 'Type', 'Iterations', 'TimeUnit'],
            columns='Stat',
            values=['RealTime', 'CpuTime']
        )
        data_frame_pivot.columns = [f'{col[0]}_{col[1]}' for col in data_frame_pivot.columns]
        data_frame_pivot = data_frame_pivot.reset_index()

        unique_benchmarks = data_frame_pivot['Benchmark'].unique()


        for bench_name in unique_benchmarks:
            bench_data = data_frame_pivot[data_frame_pivot['Benchmark'] == bench_name]
            unit = bench_data['TimeUnit'].iloc[0]
            plot_benchmark(benchmark_name=(class_name+ " " + bench_name + " Benchmark"), time_unit=unit, data=bench_data)

if __name__ == "__main__":
    main(sys.argv)
