import numpy as np
from numpy.core.numeric import ones
import pandas as pd
import csv as csv
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import PercentFormatter
rcParams.update({'figure.autolayout': True})

# --------------------------------------------------------------------------
# Helper Methods
# --------------------------------------------------------------------------
# Write data to a file
def write_data(data_list, file_name):
    """Helper method that writes each item in a list on a new line to a specified file.

    Args:
        data_list (list): list of data to write
        file_name: full-name of file to write to

    Returns:
        None, writes data to file
    """

    with open(file_name, 'a', newline='') as csvFile:
        for i in range(len(data_list)):
            writer = csv.writer(csvFile)
            writer.writerow(data_list[i])
    csvFile.close()

# Extract relevant parameters from file name
def extract_file_name_header(file_name):
    """Extracts relevant simulation data written from scaffold data structure.

    Args:
        file_name: file name to extract information from

    Returns:
        list of extracted data in the order [Scaffold stiffness, Pore diameter, Cell diameter, Replication rate, Ligand factor, Porosity]
    """

    extracted_header = []
    split_file_name = file_name.split("_")
    print(split_file_name)
    extracted_header.append(split_file_name[1].split("kPa")[0])
    extracted_header.append(split_file_name[2].split("PDµm")[0])
    extracted_header.append(split_file_name[3].split("CDµm")[0])
    extracted_header.append(split_file_name[4].split("RP")[0])
    extracted_header.append(split_file_name[5].split("LF%")[0])
    extracted_header.append(split_file_name[6].split("PS%")[0])

    return extracted_header


# ------------------------------------------------------------------------------------------------
# Extract Average Cell Speeds of Numerous Files
# ------------------------------------------------------------------------------------------------
def extract_file_average_cell_speed(files):
    """Extracts the average cell speeds of all of the cells in a list of data files from the scaffold data structure

    Args:
        file_names: list of files to calculate average cell speed

    Returns:
        None, writes data to the file cell_speeds.csv
    """

    data = []
    data.append(["Stiffness (kPa)", "Pore Diameter (um)", "Cell Diameter (um)", "Replication Rate (s^-1)", "Ligand Factor (%)", "Porosity (%)", "Average Cell Speed (um/hr)", "Standard Deviation"])
    for file in files:
        df = pd.read_csv(file)

        # Find the distance between pores
        distance_between_pores = float(list(df.columns)[-1])

        # Max simulation time
        max_time = max(list(df["Time (hour)"]))

        # Extract all unique cell IDs
        cell_ids = set(df["Cell ID"])

        cell_speeds = []
        for id in cell_ids:
            # Extract all data for particular cell ID
            cell_id_data = df[df["Cell ID"] == id]
            
            # Find the time in which the cell entered the simulation
            cell_time_in =  list(cell_id_data["Time (hour)"])[0]

            # Find the number of moves the cell moves over the entire simulation
            cell_moves = list(cell_id_data["Number of Cell Moves"])[-1]

            if (max_time - cell_time_in) != 0:
                average_speed = (cell_moves * distance_between_pores)/(max_time - cell_time_in)
                cell_speeds.append(average_speed)
            else:
                cell_speeds.append(0)

        # Extract File name header
        header = extract_file_name_header(file[:])

        # Calculate average speed and standard deviation
        average_cell_speed = np.mean(cell_speeds)
        standard_deviation = np.std(cell_speeds)
        header.append(average_cell_speed)
        header.append(standard_deviation)
        
        # Add new data to file
        data.append(header)

    write_data(data, "overall_cell_speeds.csv")


# ------------------------------------------------------------------------------------------------
# Extract Average Cell Speeds for one file
# ------------------------------------------------------------------------------------------------
def extract_each_cell_speed(file):
    """Extract the cell speed of each individual cell from data written by the scaffold data structure

    Args:
        file: file to extract data from

    Returns:
        None, writes data to the file "individual_cell_speeds.csv"
    """

    data = []
    data.append(["Time (hour)", "Cell Generation", "Average Cell Speed (um/hr)", "Standard Deviation"])
    org_df = pd.read_csv(file)
    for gen in range(1,4):
        df = org_df[org_df["Cell Generation"] == gen]

        # Find the distance between pores
        distance_between_pores = float(list(df.columns)[-1])

        # Extract all unique cell IDs
        cell_ids = set(df["Cell ID"])

        for t in range(200, 2200, 200):
            print(t)
            cell_speeds = []
            for id in cell_ids:
                # Extract all data for particular cell ID
                cell_id_data = df[df["Cell ID"] == id]
                
                # Find the time in which the cell entered the simulation
                cell_time_in = list(cell_id_data["Time (hour)"])[0]

                if cell_time_in >= t:
                    continue

                # Find the interested time
                cell_time_interested_data = cell_id_data[cell_id_data["Time (hour)"] == t]


                # Find the number of moves the cell moves over the entire simulation
                cell_moves = list(cell_time_interested_data["Number of Cell Moves"])[0]

                if (t - cell_time_in) != 0:
                    average_speed = (cell_moves * distance_between_pores)/(t - cell_time_in)
                    cell_speeds.append(average_speed)
                else:
                    cell_speeds.append(0)
            data.append([t, gen, np.mean(cell_speeds), np.std(cell_speeds)])

        # Write data to file
        write_data(data, "individual_cell_speeds.csv")

        # Clear data for each generation
        data = []


# ------------------------------------------------------------------------------------------------
# Extract Average Cell Speeds for one file at the end of the simulation
# ------------------------------------------------------------------------------------------------
def extract_cell_speed_end_simulation(file):
    """Extract the cell speed of each individual cell from data written by the scaffold data structure

    Args:
        file: file to extract data from

    Returns:
        None, writes data to the file "end_simulation_cell_speeds.csv"
    """


    write_data([["Time (hour)", "Cell Generation", "Average Cell Speed (um/hr)"]], "end_simulation_cell_speeds.csv")
    df = pd.read_csv(file)

    # Find the distance between pores
    distance_between_pores = float(list(df.columns)[-1])

    # Extract all unique cell IDs
    cell_ids = set(df["Cell ID"])

    cell_speeds = []
    for id in cell_ids:
        # Extract all data for particular cell ID
        cell_id_data = df[df["Cell ID"] == id]
        
        # Find the time in which the cell entered the simulation
        cell_time_in = list(cell_id_data["Time (hour)"])[0]

        # Find the interested time
        cell_time_interested_data = cell_id_data[cell_id_data["Time (hour)"] == 2000]

        # Find the number of moves the cell moves over the entire simulation
        cell_moves = list(cell_time_interested_data["Number of Cell Moves"])[0]

        if (2000 - cell_time_in) != 0:
            average_speed = (cell_moves * distance_between_pores)/(2000 - cell_time_in)
            cell_speeds.append(average_speed)
        else:
            cell_speeds.append(0)

        # Write data to file
        write_data([[2000, list(cell_time_interested_data["Cell Generation"])[0], np.mean(cell_speeds)]], "end_simulation_cell_speeds.csv")


# ------------------------------------------------------------------------------------------------
# Plot Line Graph
# ------------------------------------------------------------------------------------------------
def plot_generation_speeds_chart(file, x_values, y_values, xlabel, ylabel, chart_name):
    """Plots a line graph based on x_values and y_values

    Args:
        file: file to extract data from
        x_values: x values
        y_values: y values
        xlabel: x-axis label
        y-label: y-axis label

    Returns:
        None, writes data to the file cell_speeds_over_time.csv
    """

    df = pd.read_csv(file)

    # Modify the axis for inversion
    ax = plt.gca()
    ax.tick_params(width=2, length=5)

    # Set figure size and dpi
    fig = plt.figure(figsize=[4.65, 5], dpi=2000)
    sp = fig.add_subplot(1, 1, 1)
    # Set "spine" width
    for axis in ['top','bottom','left','right']:
        sp.spines[axis].set_linewidth(2)

    # Plot relevant data
    for i in range(1, 4):
        gen_df = df[df["Cell Generation"] == i]
        x = list(gen_df[x_values])
        y = list(gen_df[y_values])
        sd = list(gen_df["Standard Deviation"])
        ax.errorbar(x,y,yerr=sd, label=str(i), capsize=8)

    # Set legend
    plt.legend(title="Cell Generation", fontsize=16, title_fontsize=18, loc=1)

    # Set tick fonts
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    
    # Set label fonts
    plt.xlabel(xlabel, fontsize=18)
    plt.ylabel(ylabel, fontsize=18)

    # Show figure and save it
    plt.tight_layout()
    plt.savefig(chart_name, bbox_inches="tight")


# ------------------------------------------------------------------------------------------------
# Plot Heat Map
# ------------------------------------------------------------------------------------------------
def heat_map_cell_speeds(file, x_values, y_values, xlabel, ylabel, chart_name):
    """Plots jet-based kaiser interpolated heatmap

    Args:
        file: file to extract data from
        x_values: x values
        y_values: y values
        xlabel: x-axis label
        ylabel: y-axis label
        chart_name: name of file to save to

    Returns:
        None, saves figure to specified chart_name
    """

    df = pd.read_csv(file)

    # -------------------------
    # Modify plot(s) and sizes
    # -------------------------
    # Set figure size and dpi
    fig = plt.figure(figsize=[4.65, 5], dpi=2000)
    sp = fig.add_subplot(1, 1, 1)
    # Set "spine" width
    for axis in ['top','bottom','left','right']:
        sp.spines[axis].set_linewidth(2)

    plot_data = pd.pivot_table(df, values="Average Cell Speed (um/hr)", index=y_values, columns=x_values)
    plot = plt.imshow(plot_data, interpolation='kaiser', cmap='jet')

    # -------------------------
    # Modify the color bar
    # -------------------------
    color_bar = plt.colorbar(plot, fraction=0.046, pad=0.04)
    color_bar.set_label("µm/hr", size=18)
    color_bar.ax.tick_params(labelsize=16)

    # --------------------------------------
    # Modify the axis for inversion
    # --------------------------------------
    ax = plt.gca()
    ax.invert_yaxis()
    ax.tick_params(width=2, length=5)

    # --------------------------------------
    # Set X-axis
    # --------------------------------------
    xtick_labels = list(set(list(df[x_values])))
    xtick_labels.sort()
    xtick_labels.insert(0, 0)
    plt.locator_params(axis='x', nbins=len(xtick_labels))
    ax.set_xticklabels(xtick_labels)
    
    # --------------------------------------
    # Set Y-axis
    # --------------------------------------
    ytick_labels = list(set(list(df[y_values])))
    ytick_labels.sort()
    ytick_labels.insert(0, 0)
    plt.locator_params(axis='y', nbins=len(ytick_labels))
    ax.set_yticklabels(ytick_labels)

    # --------------------------------------
    # Modify X and Y axis label and properties
    # --------------------------------------
    # Set tick fonts
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    
    # Set label fonts
    plt.xlabel(xlabel, fontsize=18)
    plt.ylabel(ylabel, fontsize=18)

    # Show figure and save it
    plt.tight_layout()
    plt.savefig(chart_name, bbox_inches="tight")


# ------------------------------------------------------------------------------------------------
# Plot Histogram
# ------------------------------------------------------------------------------------------------
def plot_histogram(file, xlabel, ylabel, bins, number_of_generations, chart_name):
    """Plots a histogram of cell speeds for all generations in file

    Args:
        file: file to extract data from
        xlabel: x-axis label
        ylabel: y-axis label
        bins: bins to to distribute cell_speeds
        chart_name: name of file to save to

    Returns:
        None, saves figure to specified chart_name
    """

    df = pd.read_csv(file)

    # --------------------------------------
    # Plot each generation
    # --------------------------------------
    # Store X and Y values
    x_values = []
    y_values = []
    
    for gen in range(1, number_of_generations+1):
        cell_gen_df = df[df["Cell Generation"] == gen]
        cell_gen_speeds = list(cell_gen_df["Average Cell Speed (um/hr)"])
        [array, bins, patches] =plt.hist(cell_gen_speeds, weights=np.ones(len(cell_gen_speeds))/len(cell_gen_speeds), bins=bins, label=str(gen), histtype='step')
        plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
        x_values.append(bins)
        y_values.append(array)

    # -------------------------
    # Modify plot(s) and sizes
    # -------------------------
    # Set figure size and dpi
    fig = plt.figure(figsize=[4.65, 5], dpi=2000)
    sp = fig.add_subplot(1, 1, 1)

    # Set "spine" width
    for axis in ['top','bottom','left','right']:
        sp.spines[axis].set_linewidth(2)

    # Modify the axis
    ax = plt.gca()
    ax.tick_params(width=2, length=5)
    
    # -------------------------
    # Plot graph
    # -------------------------
    for i in range(len(x_values)):
        plt.errorbar(x_values[i][1:], y_values[i], label=i+1, capsize=8)
        plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    
    # Set legend
    plt.legend(title="Cell Generation", fontsize=16, title_fontsize=18, loc=1)

    # --------------------------------------
    # Modify X and Y axis label and properties
    # --------------------------------------
    # Set tick fonts
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    
    # Set label fonts
    plt.xlabel(xlabel, fontsize=18)
    plt.ylabel(ylabel, fontsize=18)

    # Show figure and save it
    plt.tight_layout()
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.legend(title="Cell Generation")
    plt.savefig(chart_name, box_inches="tight")


# --------------------------------------
# Method Calls
# --------------------------------------
# Stiffness vs. Pore Diameter - no replication
# stiffness = [2, 4, 6, 8]
# pore_diameter = [60, 80, 100, 120]
# data_file_names = []
# for p in pore_diameter:
#     for s in stiffness:
#         data_file_names.append("Sim_" + str(s) + "kPa_" + str(p) + "PDµm_35CDµm_0RP_100LF%_91PS%.csv")

# extract_file_average_cell_speed(data_file_names)

# heat_map_cell_speeds("overall_cell_speeds.csv", "Stiffness (kPa)", "Pore Diameter (um)", "Stiffness (kPa)", "Pore Diameter (µm)", "Heatmap - Cell Speed by Pore Diameter versus Stiffness - no rep.png")


# Pore Diameter vs. Cell Diameter - with replication
# pore_diameter = [60, 80, 100, 120]
# cell_diameter = [10, 20, 30, 40]
# data_file_names = []
# for p in pore_diameter:
#         for c in cell_diameter:
#             data_file_names.append("Sim_6.3kPa_" + str(p) + "PDµm_" + str(c) + "CDµm_1E-07RP_100LF%_91PS%.csv")

# extract_file_average_cell_speed(data_file_names)

# heat_map_cell_speeds("overall_cell_speeds.csv", "Pore Diameter (um)", "Cell Diameter (um)", "Pore Diameter (µm)", "Cell Diameter (µm)", "Heatmap - Cell Speed by Pore Diameter versus Cell Diameter.png")


# Pore Diameter vs. Porosity - with replication
# pore_diameter = [60, 80, 100, 120]
# porosity = [50, 65, 80 ,95]
# data_file_names = []
# for p in pore_diameter:
#         for psty in porosity:
#             data_file_names.append("Sim_6.3kPa_" + str(p) + "PDµm_35CDµm_1E-07RP_100LF%_" +  str(psty) + "PS%.csv")

# extract_file_average_cell_speed(data_file_names)

# heat_map_cell_speeds("overall_cell_speeds.csv", "Pore Diameter (um)", "Porosity (%)", "Pore Diameter (µm)", "Porosity (%)", "Heatmap - Cell Speed by Pore Diameter versus Porosity.png")


# Pore Diameter vs. Replication Rate
# pore_diameter = [60, 80, 100, 120]
# rep_rate = [1E-7, 5E-7, 1E-6, 5E-6]
# data_file_names = []
# for p in pore_diameter:
#         for r in rep_rate:
#             data_file_names.append("Sim_6.3kPa_" + str(p) + "PDµm_35CDµm_" + str(r) + "RP_100LF%_91PS%.csv")

# extract_file_average_cell_speed(data_file_names)

# heat_map_cell_speeds("overall_cell_speeds.csv", "Pore Diameter (um)", "Replication Rate (s^-1)", "Pore Diameter (µm)", "Replication Rate (s$^{-1}$)", "Heatmap - Cell Speed by Pore Diameter versus Replication Rate.png")


# Stiffness vs. Pore Diameter - with replication
# stiffness = [2, 4, 6, 8]
# pore_diameter = [60, 80, 100, 120]
# data_file_names = []
# for p in pore_diameter:
#     for s in stiffness:
#         data_file_names.append("Sim_" + str(s) + "kPa_" + str(p) + "PDµm_35CDµm_1e-07RP_100LF%_91PS%.csv")

# extract_file_average_cell_speed(data_file_names)

# heat_map_cell_speeds("overall_cell_speeds.csv", "Stiffness (kPa)", "Pore Diameter (um)", "Stiffness (kPa)", "Pore Diameter (µm)", "Heatmap - Cell Speed by Pore Diameter versus Stiffness - with rep.png")


# Stiffness vs. Cell Diameter - with replication
# stiffness = [2, 4, 6, 8]
# cell_diameter = [10, 20, 30, 40]
# data_file_names = []
# for s in stiffness:
#         for c in cell_diameter:
#             data_file_names.append("Sim_" + str(s) + "kPa_100PDµm_" + str(c) + "CDµm_1E-07RP_100LF%_91PS%.csv")

# extract_file_average_cell_speed(data_file_names)

# heat_map_cell_speeds("overall_cell_speeds.csv", "Stiffness (kPa)", "Cell Diameter (um)", "Stiffness (kPa)", "Cell Diameter (µm)", "Heatmap - Cell Speed by Stifffness versus Cell Diameter.png")

# Stiffness vs. Ligand Factor - with replication
# stiffness = [2, 4, 6, 8]
# ligand_factor = [40, 60, 80, 100]
# data_file_names = []
# for s in stiffness:
#         for lf in ligand_factor:
#             data_file_names.append("Sim_" + str(s) + "kPa_100PDµm_35CDµm_1E-07RP_" + str(lf) + "LF%_91PS%.csv")

# extract_file_average_cell_speed(data_file_names)

# heat_map_cell_speeds("overall_cell_speeds.csv", "Stiffness (kPa)", "Ligand Factor (%)", "Stiffness (kPa)", "Ligand Factor (%)", "Heatmap - Cell Speed by Stifffness versus Ligand Factor.png")


# Stiffness vs. Replication Rate
# stiffness = [2, 4, 6, 8]
# rep_rate = [1E-7, 5E-7, 1E-6, 5E-6]
# data_file_names = []
# for s in stiffness:
#         for r in rep_rate:
#             data_file_names.append("Sim_" + str(s) + "kPa_100PDµm_35CDµm_" + str(r) + "RP_100LF%_91PS%.csv")

# extract_file_average_cell_speed(data_file_names)

# heat_map_cell_speeds("overall_cell_speeds.csv", "Stiffness (kPa)", "Replication Rate (s^-1)", "Stiffness (kPa)", "Replication Rate (s$^{-1}$)", "Heatmap - Cell Speed by Stifffness versus Replication Rate.png")


# Cell Speeds Over Time By Generation
# file = "Sim_4kPa_100PDµm_35CDµm_5e-06RP_100LF%_91PS%.csv"
# extract_each_cell_speed(file)

# plot_generation_speeds_chart("individual_cell_speeds.csv", "Time (hour)", "Average Cell Speed (um/hr)", "Time (hour)", "Average Cell Speed (µm/hr)", "Line Graph - Average Speeds of Different Generations Over Time.png")


# Histogram - Cell Speeds
file = "Sim_6.3kPa_100PDµm_35CDµm_1e-05RP_100LF%_91PS%.csv"
extract_cell_speed_end_simulation(file)

plot_histogram("end_simulation_cell_speeds.csv", "Average Cell Speed (µm/hr)", "Percent of Cells (%)", [7, 8, 9, 10, 11, 12], 3, "Histogram - Distribution of Average Cell Speeds by Generation")