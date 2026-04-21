import json
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
from matplotlib import rcParams
import warnings
from matplotlib.lines import Line2D



# Use LaTeX for all text rendering (brings in Computer Modern automatically)
rcParams["text.usetex"] = True
rcParams["font.family"] = "serif"
rcParams["font.serif"] = ["Computer Modern"]

# Optional for clean math mode rendering
rcParams["mathtext.fontset"] = "cm"
# === Increase font sizes ===
rcParams["font.size"] = 16         # base font size
rcParams["axes.titlesize"] = 20    # title
rcParams["axes.labelsize"] = 24    # x and y labels
rcParams["xtick.labelsize"] = 18   # x tick labels
rcParams["ytick.labelsize"] = 18   # y tick labels
rcParams["legend.fontsize"] = 18   # legend text
rcParams["figure.titlesize"] = 24  # figure title
def plot_optimization_heatmap(optimization_results_data, flag): # flag = 1 if LCOH
    # Extract data points, filtering out any entries where essential values are None
    xi_starts = []
    eta_overpotentials = []
    objective_values = []

    for item in optimization_results_data:
        obj_val = item.get('objective_value')
        xi_start = item.get('xi_start')
        eta_overpotential = item.get('eta_overpotential')

        if obj_val is not None and xi_start is not None and eta_overpotential is not None:
            xi_starts.append(xi_start*100)
            eta_overpotentials.append(eta_overpotential*1000000)
            objective_values.append(obj_val)
        else:
            print(f"Warning: Skipping data point due to missing essential values: {item}")

    if not xi_starts:
        print("No valid data points found after filtering for plotting.")
        return

    # Convert lists to NumPy arrays
    points = np.array([xi_starts, eta_overpotentials]).T
    values = np.array(objective_values)

    # Define the grid for interpolation
    # Create a dense grid spanning the range of your parameters
    grid_x = np.linspace(min(xi_starts), max(xi_starts), 100) # 100 points for x-axis
    grid_y = np.linspace(min(eta_overpotentials), max(eta_overpotentials), 100) # 100 points for y-axis
    
    # Use meshgrid to create 2D arrays for the grid
    grid_x_mesh, grid_y_mesh = np.meshgrid(grid_x, grid_y)

    # Perform linear interpolation of the objective_values onto the new grid
    grid_z = griddata(points, values, (grid_x_mesh, grid_y_mesh), method='linear')

    # Create the heatmap
    plt.figure(figsize=(10, 8)) # Set figure size for better readability
    
    # Use pcolormesh for flexibility with non-uniform grids, though imshow also works for regular grids
    # vmin and vmax can be set to control the color scale range
    # cmap can be chosen to represent the data effectively (e.g., 'viridis', 'plasma', 'magma', 'hot')
    if flag:
        c = plt.pcolormesh(grid_x_mesh, grid_y_mesh, grid_z, shading='auto', cmap='viridis')
    else:
        c = plt.pcolormesh(grid_x_mesh, grid_y_mesh, grid_z, shading='auto', cmap='viridis')

    # Add labels and title
    
    custom_xticks = np.linspace(min(xi_starts), max(xi_starts), 7) # 5 ticks evenly spaced
    plt.xticks(custom_xticks)
    
    custom_yticks = np.linspace(min(eta_overpotentials), max(eta_overpotentials), 12) # 7 ticks evenly spaced
    plt.yticks(custom_yticks)
    
    plt.xlabel(r"Electrical Efficiency Start $(\% LHV)$")
    plt.ylabel(r"Overpotential Per Hour $(\mu V/hr)$")
    # plt.title("Optimization Objective Value Heatmap (Smoothed)", fontsize=14)

    # Add a color bar to indicate the objective value
    cbar = plt.colorbar(c)
    if flag:
        cbar.set_label("LCOH")
    else: 
        cbar.set_label("Profit")

    # Improve layout and display the plot
    plt.tight_layout()
    # plt.show()


def plot_optimization_heatmap_unsmoothed(optimization_results_data):
    if not optimization_results_data:
        print("No data provided for plotting the unsmoothed heatmap.")
        return

    # Extract data points, filtering out any entries where essential values are None
    xi_starts = []
    eta_overpotentials = []
    objective_values = []

    for item in optimization_results_data:
        obj_val = item.get('objective_value')
        xi_start = item.get('xi_start')
        eta_overpotential = item.get('eta_overpotential')

        if obj_val is not None and xi_start is not None and eta_overpotential is not None:
            xi_starts.append(xi_start)
            eta_overpotentials.append(eta_overpotential)
            objective_values.append(obj_val)
        else:
            print(f"Warning: Skipping data point due to missing essential values: {item}")

    if not xi_starts:
        print("No valid data points found after filtering for plotting.")
        return

    # Create the scatter plot
    plt.figure(figsize=(10, 8)) # Set figure size
    
    # Use scatter plot where color ('c') is mapped to objective_values
    # s: marker size (e.g., 100 for a visible point)
    # cmap: colormap (e.g., 'viridis')
    scatter = plt.scatter(
        xi_starts,
        eta_overpotentials,
        c=objective_values,
        s=100, # Adjust marker size as needed
        cmap='viridis', # Choose a colormap
        edgecolor='black', # Optional: Add a black edge to points
        alpha=0.8 # Optional: Adjust transparency
    )

    # Add labels and title
    # custom_xticks = np.linspace(min(xi_starts), max(xi_starts), 10) # 5 ticks evenly spaced
    # plt.xticks(custom_xticks, rotation=45)
    # custom_yticks = np.linspace(min(eta_overpotentials), max(eta_overpotentials), 7) # 7 ticks evenly spaced
    # # Or, for specific values: custom_yticks = [0.00005, 0.00006, 0.00007, 0.00008, 0.00009, 0.0001]
    # plt.yticks(custom_yticks)
    plt.xlabel("Xi Start (ξ_start)", fontsize=12)
    plt.ylabel("Eta Overpotential (η_overpotential)", fontsize=12)
    plt.title("Optimization Objective Value (Unsmoothed)", fontsize=14)

    # Add a color bar
    cbar = plt.colorbar(scatter)
    cbar.set_label("Objective Value", fontsize=12)

    # Adjust limits to ensure all points are visible and add some padding
    plt.xlim(min(xi_starts) * 0.98, max(xi_starts) * 1.02)
    plt.ylim(min(eta_overpotentials) * 0.98, max(eta_overpotentials) * 1.02)

    plt.grid(True, linestyle='--', alpha=0.6) # Optional: Add a grid
    plt.tight_layout()
    plt.show()

    
def extract_plotting_data(path):
    """
    Extracts LCOH and profit from one JSON file
    Returns:
        dict or None: A dictionary containing 'profit', 'lcoh', 'xi_start',
                      and 'eta_overpotential' if successful, otherwise None.
                      Returns None if files are not found, cannot be decoded,
                      are missing essential data, or contain non-finite values.
    """

    profit_data = None

    # Load LCOH data
    if not os.path.exists(path):
        warnings.warn(f"LCOH file '{path}' not found. Skipping this data point.")
        return None
    try:
        with open(path, 'r', encoding='utf-8') as f:
            lcoh_json = json.load(f)
        lcoh_data = lcoh_json.get("LCOH")
        profit_data = lcoh_json.get("objective_value")
        # Extract common parameters (xi_start, eta_overpotential) from LCOH file
        electrolyzer_data = lcoh_json.get("el", {})
        xi_start = electrolyzer_data.get("ξ_start")
        eta_overpotential = electrolyzer_data.get("η_overpotential")
    except json.JSONDecodeError as e:
        warnings.warn(f"Error decoding LCOH JSON from {path}: {e}. Skipping.")
        return None
    except Exception as e:
        warnings.warn(f"An unexpected error occurred while processing LCOH file {path}: {e}. Skipping.")
        return None

    return {
        "profit": float(profit_data),
        "lcoh": float(lcoh_data),
        "xi_start": float(xi_start),
        "eta_overpotential": float(eta_overpotential)
    }

def plot_profit_and_lcoh_heatmap_with_contours(plotting_data):
    """
    Generates a heatmap for Profit with overlaid contour lines for LCOH.

    The Profit values define the color gradient of the heatmap.
    The LCOH values are represented by overlaid contour lines.

    Applies custom rcParams for LaTeX rendering and font sizes.

    Args:
        plotting_data (list of dict): A list where each dictionary contains
                                      'profit', 'lcoh', 'xi_start', and 'eta_overpotential'.
                                      This list should be generated by processing your JSON files.
    """
    if not plotting_data:
        print("No valid data provided for plotting the heatmap with contours.")
        return

    # Extract and scale data points
    xi_starts = []
    eta_overpotentials = []
    lcoh_values = []
    profit_values = []

    for item in plotting_data:
        # Values are scaled as per your original heatmap code
        xi_starts.append(item['xi_start'] * 100)
        eta_overpotentials.append(item['eta_overpotential'] * 1000000)
        lcoh_values.append(item['lcoh'])
        profit_values.append(item['profit']/1000000)

    # Convert lists to NumPy arrays for interpolation
    points = np.array([xi_starts, eta_overpotentials]).T
    lcoh_np = np.array(lcoh_values)
    profit_np = np.array(profit_values)

    # Define the grid for interpolation
    grid_x = np.linspace(min(xi_starts), max(xi_starts), 100)
    grid_y = np.linspace(min(eta_overpotentials), max(eta_overpotentials), 100)
    grid_x_mesh, grid_y_mesh = np.meshgrid(grid_x, grid_y)

    # Perform linear interpolation for both Profit and LCOH onto the new grid
    grid_profit = griddata(points, profit_np, (grid_x_mesh, grid_y_mesh), method='linear')
    grid_lcoh = griddata(points, lcoh_np, (grid_x_mesh, grid_y_mesh), method='linear')

    plt.figure(figsize=(10, 8))
    
    # --- Plot Heatmap for Profit ---
    # Choose a sequential colormap where higher values (more profit) are visually emphasized.
    # 'viridis', 'plasma', 'magma', 'inferno' are good options.
    heatmap = plt.pcolormesh(grid_x_mesh, grid_y_mesh, grid_profit, shading='auto', cmap='viridis')
    cbar = plt.colorbar(heatmap)
    cbar.set_label(r"Profit $(\$MM)$") # Label the color bar for Profit

    # --- Overlay Contour Lines for LCOH ---
    # Define `levels` to show meaningful LCOH values. Lower LCOH is better.
    # You might want to choose levels that highlight specific cost targets.
    lcoh_levels = np.linspace(min(lcoh_values), max(lcoh_values), 8) # Example: 6 equally spaced levels
    # Alternatively, specify custom levels:
    # lcoh_levels = [2.5, 3.0, 3.5, 4.0] # Example for $/kg

    contours = plt.contour(grid_x_mesh, grid_y_mesh, grid_lcoh, levels=lcoh_levels, colors='red', linestyles='-') # Use red for contrast
    # Add labels to the contour lines for readability
    # plt.clabel(contours, inline=True, fontsize=20, fmt='%1.2f') 

    # --- Add Labels and Title ---
    plt.xlabel(r"Initial Electrical Efficiency $(\% LHV)$")
    plt.ylabel(r"Operation Degradation Rate $(\mu V/hr)$")
    plt.title(r"$11$-Year NPV Heatmap")

    # --- Customize X and Y Ticks ---
    custom_xticks = np.linspace(min(xi_starts), max(xi_starts), 7)
    plt.xticks(custom_xticks)
    
    custom_yticks = np.linspace(min(eta_overpotentials), max(eta_overpotentials), 11)
    plt.yticks(custom_yticks)
    
    from matplotlib.lines import Line2D
    lcoh_legend_line = Line2D([0], [0], color='red', linestyle='-', label='LCOH Isolines')
    
    labels = plt.clabel(contours, inline=True, fontsize=20, fmt='%1.2f')
    
    for l in labels:
        l.set_fontweight('bold')
        l.set_bbox(dict(facecolor='white', edgecolor='none', alpha=0.8, pad=2))


    # Add the legend to the plot
    plt.legend(handles=[lcoh_legend_line], loc='best', fontsize=plt.rcParams["legend.fontsize"])


    # Improve layout and display the plot
    plt.tight_layout()
    
    plt.savefig("heatmap_lcoh_iso.png", format='png', transparent=True)
    plt.show()


def plot_profit_and_lcoh_heatmap_with_patterns(plotting_data):
    # ... [Keep your existing data extraction and interpolation code here] ...
    if not plotting_data:
        print("No valid data provided for plotting the heatmap with contours.")
        return

    # Extract and scale data points
    xi_starts = []
    eta_overpotentials = []
    lcoh_values = []
    profit_values = []

    for item in plotting_data:
        # Values are scaled as per your original heatmap code
        xi_starts.append(item['xi_start'] * 100)
        eta_overpotentials.append(item['eta_overpotential'] * 1000000)
        lcoh_values.append(item['lcoh'])
        profit_values.append(item['profit']/1000000)

    # Convert lists to NumPy arrays for interpolation
    points = np.array([xi_starts, eta_overpotentials]).T
    lcoh_np = np.array(lcoh_values)
    profit_np = np.array(profit_values)

    # Define the grid for interpolation
    grid_x = np.linspace(min(xi_starts), max(xi_starts), 100)
    grid_y = np.linspace(min(eta_overpotentials), max(eta_overpotentials), 100)
    grid_x_mesh, grid_y_mesh = np.meshgrid(grid_x, grid_y)

    # Perform linear interpolation for both Profit and LCOH onto the new grid
    grid_profit = griddata(points, profit_np, (grid_x_mesh, grid_y_mesh), method='linear')
    grid_lcoh = griddata(points, lcoh_np, (grid_x_mesh, grid_y_mesh), method='linear')


    # Define the grid (reusing your grid_x_mesh, grid_y_mesh, grid_profit, grid_lcoh)
    # ... [Assuming grid_x_mesh, grid_y_mesh, grid_profit, grid_lcoh are processed] ...

    plt.figure(figsize=(12, 9))

    # --- 1. Reverse Colormap (Darker = Less Profit) ---
    # Adding '_r' to any colormap name reverses it. 
    # 'Greys_r' starts dark and goes light.
    profit_levels = np.linspace(np.nanmin(grid_profit), np.nanmax(grid_profit), 8)
    hatches = ['XX', 'x', '//', '/', '..', '.', '', ''] # Reversed hatches to match
    
    profit_contour = plt.contourf(grid_x_mesh, grid_y_mesh, grid_profit, 
                                  levels=profit_levels, 
                                  cmap='Greys_r', 
                                  hatches=hatches, 
                                  alpha=0.8)
    
    cbar = plt.colorbar(profit_contour)
    cbar.set_label(r"Profit $(\$MM)$")

    # --- 2. Enhanced LCOH Lines ---
    lcoh_levels = np.linspace(np.nanmin(grid_lcoh), np.nanmax(grid_lcoh), 6)
    contours = plt.contour(grid_x_mesh, grid_y_mesh, grid_lcoh, 
                           levels=lcoh_levels, 
                           colors='black', 
                           linewidths=2.5) # Thicker lines

    # --- 3. Large Labels with White Bounding Boxes ---
    # We capture the labels in a variable to style them individually
    labels = plt.clabel(contours, inline=True, fontsize=16, fmt='%1.2f', inline_spacing=10)
    
    for l in labels:
        l.set_fontweight('bold')
        l.set_bbox(dict(facecolor='white', edgecolor='none', alpha=0.8, pad=2))

    # --- 4. Final Polish ---
    plt.xlabel(r"Electrical Efficiency Start $(\% LHV)$")
    plt.ylabel(r"Operation Degradation Rate $(\mu V/hr)$")
    plt.title(r"Profit Analysis with High-Visibility LCOH Contours")
    
    # Adding legend for the lines
    custom_lines = [Line2D([0], [0], color='black', lw=2.5)]
    plt.legend(custom_lines, ['LCOH Isolines'], loc='upper right')

    plt.tight_layout()
    plt.savefig("heatmap_mono.png", format='png', transparent=True)
    plt.show()

if __name__ == "__main__":
    folder_path = "data/heatmap"

    efficiencies = ["6", "65", "7", "75"]
    suffixes = ["_0.0", "_1.6", "_3.2", "_4.8", "_6.4", "_8.0", "_10.0"]

    filenames = ["eff_0." + e + "_deg" + s + "_results.JSON" for e in efficiencies for s in suffixes]

    all_results = []
    for f in filenames:
        res = extract_plotting_data(folder_path+f)
        all_results.append(res)

    plot_profit_and_lcoh_heatmap_with_contours(all_results)
    # plot_profit_and_lcoh_heatmap_with_patterns(all_results)