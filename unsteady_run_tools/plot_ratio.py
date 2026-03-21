import matplotlib.pyplot as plt
import re
import sys
import glob
import os

def plot_log(case_folder):
    # 1. Find the .out file inside the folder
    # Matches: case_folder/case_folder.out.12345
    search_pattern = os.path.join(case_folder, f"*.out.*")
    file_list = glob.glob(search_pattern)

    if not file_list:
        print(f"Error: No log file found matching {search_pattern}")
        return
    
    # If multiple job files exist, pick the newest one
    file_path = max(file_list, key=os.path.getmtime)
    print(f"Reading from: {file_path}")

    times = []
    ratios = []
    pattern = re.compile(r"t=([\d\.]+)\s+L1=[\d\.]+\s+ratio=([\d\.]+)")

    try:
        with open(file_path, 'r') as f:
            for line in f:
                match = pattern.search(line)
                if match:
                    times.append(float(match.group(1)))
                    ratios.append(float(match.group(2)))
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    if not times:
        print("No valid data found (simulation might not have started yet).")
        return

    # --- Plotting ---
    plt.figure(figsize=(10, 6))
    plt.plot(times, ratios, marker='o', ms=2, label=f'Residual Ratio ({case_folder})')
    
    plt.yscale('log') # Use log scale for the Y-axis
    plt.xlabel('Time (t)')
    plt.ylabel('Residual Ratio (Log Scale)')
    plt.title(f'Convergence History for {case_folder}')
    plt.grid(True, which="both", ls="-", alpha=0.3)
    plt.legend()
    plt.show()
    
    # Save image in the current directory (outside the case folder)
    #plt.savefig(f"plot_{case_folder}.png", dpi=300)
    #print(f"Success! Plot saved as plot_{case_folder}.png")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_ratio.py <case_folder_name>")
    else:
        plot_log(sys.argv[1])
