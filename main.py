import csv
import random as r
from scaffold import Scaffold

# --------------------------------------------------------------------------
# Helper Method
# --------------------------------------------------------------------------
def write_data(data_list, file_name):
    """TODO"""

    with open(file_name, 'a', newline='') as csvFile:
        for i in range(len(data_list)):
            writer = csv.writer(csvFile)
            writer.writerow(data_list[i])
    csvFile.close()


# --------------------------------------------------------------------------
# Run Simulation
# --------------------------------------------------------------------------
def simulation(header, file_name, scaffold, seeded_cell_count, time_step, time_final, writing_frequency, recording_frequency, replication_probability):
    """TODO"""
    
    # Data Storage
    data = []

    # Data headers to CSV file
    header[0].append(scaffold.pore_array[1] - scaffold.pore_array[0])
    write_data(header, file_name)

    # Generate cells in scaffold
    scaffold.seed_cells(seeded_cell_count, data, 0)

    # Start of time increment stepping in the simulation
    time = time_step
    while time <= time_final:
        print(time)

        # Migration and replication of scaffold
        scaffold.migration_replication(data, time, time_step, time_final, recording_frequency, replication_probability * 3600 * time_step)

        if (time % writing_frequency == 0 or time == time_final):
            # Record data in list
            write_data(data, file_name)
            # Clear out data list after recording
            data = []

        print(scaffold.cell_count)

        # Add time step to increment simulation
        time += time_step # hours
        time = round(time, 1)


# --------------------------------------------------------------------------
# Simulation Characteristics
# --------------------------------------------------------------------------

# Modifies Replication Rate (s^-1)
for rr in [1E-6]:
    # Modifies Cell Diameter (µm)
    for cd in [35]:
        # Modifies Ligand Factor (%)
        for lf in [100]:
            # Modifies Pore Diameter (µm)
            for pd in [100]:
                # Modifies Scaffold Stiffness (kPa)
                for ss in [6.3]:
                    # Modifies Porosity (%)
                    for ps in [91]:

                        # Simulation Parameters
                        time_step = 0.1                     # Amount of time the simulation "steps" forward each time to perform migration and/or replication

                        # Scaffold Parameters:
                        dimension = 1500                    # Side length of cubical Scaffold (µm)
                        porosity = ps                       # Measure of 'empty' volume in scaffold
                        pore_diameter = pd                  # Diameter length of pores in scaffold (µm)
                        scaffold_stiffness = ss             # Scaffold Stiffness (kPa)
                        ligand_factor = lf                  # Ligand Factor (%)
                        packing_density = 80                # Fraction of empty volume that can be occupied by cells (%)

                        # Cell parameters:
                        cell_diameter = cd                                                              # Cell diameter (µm)
                        initial_cell_count = round(0.5 * dimension**3 * (3.25*10**6)/(1*10**12))        # Initial seeding of scaffold with cells
                        replication_probability = rr                                                    # Probability that a cell will replicate (s^-1)

                        # Simulation Parameters:
                        simulation_time = 2000              # Time length of simulation (hour)
                        writing_frequency = 100             # Frequency of how often the program writes all stored data to the CSV (hours per write)
                        recording_frequency = 50            # Frequency of how often the program records data (hours per record)

                        # Header information
                        csv_file_name = "Sim_" + str(scaffold_stiffness) + "kPa_" + str(pore_diameter) + "PDµm_" + str(cell_diameter) + "CDµm_" + str(replication_probability) + "RP_" + str(ligand_factor) + "LF%_" + str(porosity) + "PS%" + ".csv"
                        file_header =[["Time (hour)", "Cell ID", "Cell Generation", "Number of Cell Moves", "Traction Force (N)", "X (um)", "Y (um)", "Z (um)"]]


                        # --------------------------------------------------------------------------
                        # Generate scaffold environment
                        # --------------------------------------------------------------------------
                        env_1 = Scaffold(dimension, porosity/100, pore_diameter, packing_density/100, cell_diameter, scaffold_stiffness*10**3, ligand_factor/100)


                        # --------------------------------------------------------------------------
                        # Run Simulation
                        # --------------------------------------------------------------------------
                        simulation(file_header, csv_file_name, env_1, initial_cell_count, time_step, simulation_time, writing_frequency, recording_frequency, replication_probability)

                        # Test push