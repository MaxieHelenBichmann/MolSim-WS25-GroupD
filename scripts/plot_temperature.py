#!/usr/bin/env python3
import matplotlib.pyplot as plt
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Plot temperature from simulation log.')
    parser.add_argument('logfile', nargs='?', default='build/1log.txt', help='Path to the log file')
    parser.add_argument('--out', default='temperature_plot.png', help='Output filename for the plot')
    args = parser.parse_args()

    if not os.path.exists(args.logfile):
        print(f"Error: Log file '{args.logfile}' not found.")
        return

    temperatures = []
    with open(args.logfile, 'r') as f:
        for line in f:
            if "Current Temperature:" in line:
                try:
                    temp = float(line.split(":")[1].strip())
                    temperatures.append(temp)
                except (ValueError, IndexError):
                    continue

    if not temperatures:
        print("No temperature data found in the log file.")
        return

    plt.figure(figsize=(10, 6))
    plt.plot(temperatures)
    plt.xlabel('Report Interval')
    plt.ylabel('Temperature')
    plt.title('Simulation Temperature over Time')
    plt.grid(True)
    plt.savefig(args.out)
    print(f"Plot saved to {args.out}")

if __name__ == "__main__":
    main()
