import csv
import column_scaffold

# --------------------------------------------------------------------------
# Helper Methods
# --------------------------------------------------------------------------
def write_data(data_list, file_name):
    """
    Writes each row in a data array to a CSV file.

    :param data_list: data list to write to CSV file
    :param file_name: filename to write to
    """
    with open(file_name, 'a', newline='') as csv_file:
        for i in range(len(data_list)):
            writer = csv.writer(csv_file)
            writer.writerow(data_list[i])
    csv_file.close()


# --------------------------------------------------------------------------
# Run Simulation
# --------------------------------------------------------------------------
def simulation(header, file_name, scaffold, seeded_cell_count, ts, time_final, wf,
               rf, rp):
    """
    Runs a time-step dependent simulation of cell migration and replication through a scaffold.

    :param header: list that specifies header information for the CSV
    :param file_name: file name to write to CSV
    :param scaffold: scaffold environment object to simulate
    :param seeded_cell_count: number of cells to seed into the scaffold
    :param ts: specifies the simulation's time-step (hour)
    :param time_final: specifies the simulation's end time (hour)
    :param wf: Frequency of how often the program writes all stored data to the CSV (hours/write to file)
    :param rf: Frequency of how often the program records data (hours/record to data list)
    :param rp: Fixed probability that a cell will replicate (s^-1)
    """
    # Stores data in rows in an expandable list
    data = []

    # Writes the header information to the top of the CSV
    write_data(header, file_name)

    # Generate and seed cells at a specified location within the scaffold
    scaffold.seed_cells_at_top(seeded_cell_count, data)

    # Time starts at time step after seeding
    time = ts

    # Initiation of time step simulation until specified simulation end-time
    while time <= time_final:
        print("Simulation Time: ", time)  # TODO: Remove, testing purposes only

        # Migration and replication of all cells in the designated scaffold
        scaffold.migration_replication(data, time_step, time_final, rf,
                                       rp * 3600 * time_step)

        # Write data to the data list at specified intervals
        if time % wf == 0 or time == time_final:
            # Record data in list
            write_data(data, file_name)
            # Clear out data list after recording
            data = []

        print("Cell Count: " + str(scaffold.get_cell_count()))  # TODO: Remove, testing purposes only

        # Add time step increment
        time += ts  # hours
        time = round(time, 2)


# --------------------------------------------------------------------------
# Simulation and Scaffold Characteristics
# --------------------------------------------------------------------------
# Simulation Parameters
time_step = 0.25  # Time "Step" forward each time to perform migration and/or replication (hour)

# Scaffold Parameters:
dimension = 5000                   # Side length of cubical Scaffold (µm)
porosity = 42.1                     # Measure of void or "empty" volume in scaffold (%)
scaffold_stiffness = 120            # Scaffold Stiffness (MPa)
ligand_factor = 100                 # Percentage of ligands compared to normal (%)
pore_size = 23.9                    # Diameter of pores in the scaffold (µm)
packing_density = 80                # Fraction of void or empty volume within the scaffold that can be occupied by
                                    # cells (%)
pore_layer_count = 30               # The number of layers to represent each pore column

# Cell parameters:
cell_diameter = 35                  # Cell diameter (µm)
initial_cell_count = 60000          # Initial seeding of scaffold with cells
replication_probability = 7.72E-06  # Fixed probability that a cell will replicate (s^-1)

# Simulation Parameters:
simulation_time = 168               # Total simulation time (hour)
writing_frequency = 0.25               # Frequency of how often the program writes all stored data to the
                                    # CSV (hours/write to file)
recording_frequency = 0.25             # Frequency of how often the program records data (hours/record to data list)

# Header Information
csv_file_name = "Sim_" + str(scaffold_stiffness) + "MPa_" + str(pore_size) + "PDµm_" + str(
    cell_diameter) + "CDµm_" + str(replication_probability) + "RP_" + str(ligand_factor) + "LF%_" + str(
    porosity) + "PS%_" + str(pore_layer_count) + "CC" + ".csv"
file_header = [["Time (hour)", "Cell ID", "Cell Generation", "Number of Cell Moves", "X (µm)", "Y (µm)", "Z (µm)"]]

# --------------------------------------------------------------------------
# Generate scaffold environment
# --------------------------------------------------------------------------
# Generates a scaffold based on specified parameters in "Simulation and Scaffold Characteristics"
scaffold_1 = column_scaffold.Column_Scaffold(dimension, porosity / 100, pore_size, packing_density / 100,
                                                cell_diameter, scaffold_stiffness * 10 ** 6, ligand_factor / 100,
                                                pore_layer_count)

# --------------------------------------------------------------------------
# Run Simulation
# --------------------------------------------------------------------------
simulation(file_header, csv_file_name, scaffold_1, initial_cell_count, time_step, simulation_time, writing_frequency,
           recording_frequency, replication_probability)
