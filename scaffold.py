from operator import truediv
import numpy as np
import math as m
import random as r

class Scaffold:
    def __init__(self, dimension, porosity, pore_diameter, packing_density, cell_diameter, scaffold_stiffness, ligand_factor):
        """TODO"""

        # --------------------------------------------------------------------------
        # Tracks the number of cells and cell IDs in the scaffold
        # --------------------------------------------------------------------------
        self.cell_count = 0
        self.cell_ID_counter = 1

        # --------------------------------------------------------------------------
        # Generate scaffold characteristics
        # --------------------------------------------------------------------------
        # Generates various scaffold characteristics from method parameters
        pore_array, pore_number, pore_cell_max, scaffold_cell_max = self.__generate_environment_parameters(dimension, porosity, pore_diameter, packing_density, cell_diameter)

        self.pore_number = pore_number                  # Stores the number of pores in the scaffold
        self.pore_diameter = pore_diameter              # Stores the diameters of pores in the scaffold
        self.cell_diameter = cell_diameter              # Stores the cell diameter of cells that can be seeded in scaffold
        self.pore_array = pore_array                    # Stores an array that represents the coordinates of all side lengths in the scaffold
        self.pore_cell_max = pore_cell_max              # Stores the maximum number of cells one pore in the scaffold can hold
        self.scaffold_cell_max = scaffold_cell_max      # Stores the maximum number of cells the scaffold can hold
        self.scaffold_stiffness = scaffold_stiffness    # Stores the scaffold stiffness (Pa)
        self.ligand_factor = ligand_factor              # Stores the percent ligands normally (%)     

        # --------------------------------------------------------------------------
        # Generate the pores in the scaffold and assigns their positions
        # --------------------------------------------------------------------------
        pore_array_size = len(pore_array)  
        self.scaffold = [[[self.Pore(None, 0, self.pore_cell_max, pore_diameter) for x in range(pore_array_size)] for y in range(pore_array_size)] for z in range(pore_array_size)]
        self.__assign_pore_positions(pore_array_size)

        # --------------------------------------------------------------------------
        # Intiailizes and stores cell objects that are in the scaffold
        # --------------------------------------------------------------------------
        self.cell_objs = None


    def __generate_environment_parameters(self, dimension, porosity, pore_diameter, packing_density, cell_diameter):
        """TODO"""

        # --------------------------------------------------------------------------
        # Calculates volume that is available to pores
        # --------------------------------------------------------------------------
        scaffold_volume = dimension**3                          # Scaffold volume with the assumption of a cubical scaffold (µm^3)
        scaffold_porous_volume = porosity * scaffold_volume     # Amount of volume in the cubical scaffold that is 'porous' (µm^3)
        pore_volume = 4/3*m.pi*(pore_diameter/2)**3             # Calculates the volume in one pore (assuming spherical) (µm^3)

        # --------------------------------------------------------------------------
        # Calculates the number of pores based on the pore volume avilable
        # --------------------------------------------------------------------------
        pore_number = m.floor(scaffold_porous_volume/pore_volume)  # Ensures only that all pores are uniform and have the same volume
        pore_per_side = m.floor(pore_number ** (1/3))              # Calculates the number of pores per side length, also ensures the scaffold is cubical in nature
        pore_number = pore_per_side**3                             # Re-calculates new pore number to maintain cubical scaffold

        # --------------------------------------------------------------------------
        # Create side-length array
        # --------------------------------------------------------------------------
        pore_array = np.linspace(0, dimension, num=pore_per_side).tolist()       # Creates an array that represents the coordinates of all side lengths in the scaffold

        # Calculates the maximum cells in the scaffold and maximum cells in a single pore
        cell_volume = 4/3*np.pi*(cell_diameter/2)**3                        # Calulates the volume of cells in the scaffold (assuming spherical) (µm^3)
        occupiable_pore_volume = packing_density*pore_volume                # Amount of volume in the pore that can be occupied by cells based on a packing density (µm^3)
        cell_max_per_pore =  m.floor(occupiable_pore_volume/cell_volume)     # Maximum number of cells per pore
        max_cells = pore_number*cell_max_per_pore                           # Maximum number of cells that can occupy the scaffold

        # Returns scaffold characteristics
        return pore_array, pore_number, cell_max_per_pore, max_cells


    def __assign_pore_positions(self, pore_array_size):
        """TODO"""

        pX = 0
        pY = 0
        pZ = 0
        while pZ < pore_array_size:
            while pY < pore_array_size:
                while pX < pore_array_size:
                    self.scaffold[pX][pY][pZ].position = [pX, pY, pZ]
                    pX += 1
                pX = 0
                pY += 1
            pY = 0
            pZ += 1


    def seed_cells(self, seeded_cell_count, data, seeded_time):
        """TODO"""

        if (seeded_cell_count > self.scaffold_cell_max):
            raise RuntimeError("The number of cells you are trying to seed cannot be larger than the maximum cell capacity of the scaffold.")

        # Generate seeded cell objects
        self.cell_objs = [self.Cell(self.cell_diameter, 0) for i in range(seeded_cell_count)]

        end_point = len(self.pore_array) - 1
        for cell in self.cell_objs:
            # Assign a unique ID for each seeded cell
            cell.ID = self.cell_ID_counter
            self.cell_ID_counter += 1
            self.cell_count += 1

            # Random seed cells by generating random positions, writes its properties to the data list
            while True:
                new_position = [r.randint(0, end_point), r.randint(0, end_point), r.randint(0, end_point)]
                pore_destination = self.scaffold[new_position[0]][new_position[1]][new_position[2]]
                if (not pore_destination.is_full()):
                    cell.position = new_position
                    pore_destination.cell_number += 1
                    break
            
            data.append(cell.write_data(self.pore_array, seeded_time))


    def migration_replication(self, data, t, ts, tf, rf, rep_prob):
        """TODO"""

        # Ensure unbiased cellular migration and replication
        r.shuffle(self.cell_objs)

        # Migration and replication of parent and daughter cells
        for cell in self.cell_objs:

            # Neighbor capacity checks whether neighboring pores in regards to the pore the cell is in are full. X, Y, and Z in both positive and negative direction are accounted by conditions 1-6
            # while the pore's capacity, in which the current cell is in, is accounted by condition 7
            neighbor_capacity = self.neighbor_conditional_check(cell)[:]

            # If all pores, including the pore the cell is in, are full then no movement and no replication occurs. Only recording of relevant cell properties to the data array is performed.
            if (neighbor_capacity[0] and neighbor_capacity[1] and neighbor_capacity[2] and neighbor_capacity[3] and neighbor_capacity[4] and neighbor_capacity[5] and neighbor_capacity[6]) == True:
                 # Writes cell properties to data list
                if (t % rf == 0 or t == tf):
                    data.append(cell.write_data(self.pore_array, t))
            else:
                if (neighbor_capacity[0] and neighbor_capacity[1] and neighbor_capacity[2] and neighbor_capacity[3] and neighbor_capacity[4] and neighbor_capacity[5]) == True:
                    rp = 0
                    if rep_prob != 0:
                        if round(1/rep_prob) == 1:
                            rp = 1
                        else:
                            rp = r.randint(1, round(1/rep_prob))
                    if rp != 1:
                        if (t % rf == 0 or t == tf):
                            data.append(cell.write_data(self.pore_array, t))
                    # Replication of cell occurs if replication condition is achieved
                    else:
                        new_daughter_cell = self.new_daughter_cell(cell, t)

                        # New cell with an identical pore location to its parent is added to the list of cell objects
                        self.scaffold[cell.position[0]][cell.position[1]][cell.position[2]].cell_number += 1
                        self.cell_objs.append(new_daughter_cell)

                        # Update cell count and cell ID counter
                        self.cell_ID_counter += 1
                        self.cell_count += 1

                        # Writes cell properties to data list
                        if (t % rf == 0 or t == tf):
                            data.append(cell.write_data(self.pore_array, t))
                            data.append(new_daughter_cell.write_data(self.pore_array, t))
                else:
                    # Perform migration on parent cell
                    self.__perform_migration(cell, neighbor_capacity, ts, False)

                    # Perform replication which is dependent on the replication probability
                    rp = 0
                    if rep_prob != 0:
                        if round(1/rep_prob) == 1:
                            rp = 1
                        else:
                            rp = r.randint(1, round(1/rep_prob))
                    if rp != 1:
                        if (t % rf == 0 or t == tf):
                            data.append(cell.write_data(self.pore_array, t))
                    else:
                        # Generate a new daughter cell and perform migration on it
                        new_daughter_cell = self.new_daughter_cell(cell, t)
                        daughter_neighbor_conditions = self.neighbor_conditional_check(new_daughter_cell)
                        self.__perform_migration(new_daughter_cell, daughter_neighbor_conditions, ts, True)

                        # Replicated cell must migrate to exist
                        if new_daughter_cell.moves_made != 0:
                            self.cell_objs.append(new_daughter_cell)

                            # Update cell count and cell ID counter
                            self.cell_ID_counter += 1
                            self.cell_count += 1
                        
                        # Writes cell properties to data list
                        if (t % rf == 0 or t == tf):
                            data.append(cell.write_data(self.pore_array, t))
                            if new_daughter_cell.moves_made != 0:
                                data.append(new_daughter_cell.write_data(self.pore_array, t))
        

    def __perform_migration(self, cell, neighbor_capacity, ts, is_new_daughter_cell):
        
        # --------------------------------------------------------------------------
        # Choose a random direction for migration
        # --------------------------------------------------------------------------
        migration_direction = None
        # Ensure direction is not full or out of bounds
        while True:
            migration_direction = r.randint(1,6)
            if (neighbor_capacity[0] == True and migration_direction == 1):
                continue
            elif (neighbor_capacity[1] == True and migration_direction == 2):
                continue
            elif (neighbor_capacity[2] == True and migration_direction == 3):
                continue
            elif (neighbor_capacity[3] == True and migration_direction == 4):
                continue
            elif (neighbor_capacity[4] == True and migration_direction == 5):
                continue
            elif (neighbor_capacity[5] == True and migration_direction == 6):
                continue
            else:
                break
        
        # --------------------------------------------------------------------------
        # Calculate probability to migrate in chosen direction
        # --------------------------------------------------------------------------
        # Calculate distance 
        cell_velocity = self.__calculate_cell_velocity(cell, migration_direction)           # Cell velocity in chosen direction (µm/s)
        distance_between_pores = self.pore_array[1] - self.pore_array[0]                    # Distance between neighboring pores (equidistant)
        cell_distance_traveled = cell_velocity * ts * 3600                                  # Distance cell can travel (µm)

        # Calculate probability to migrate as a ratio between distance cell can travel and distance between pores
        probability_of_migration = cell_distance_traveled/distance_between_pores * 100      # Probability of migration (%)

        # --------------------------------------------------------------------------
        # Cell Migration
        # --------------------------------------------------------------------------
        # Generates a number between 1 and 100 to determine if cell will migrate
        random_migration_number = r.randint(0, 100000)/1000

        # Migrate in chosen direction based on migration probability
        if probability_of_migration > random_migration_number:
           
            # Update number of moves made
            cell.moves_made += 1
            
            if migration_direction == 1:
                # X > 0
                if not is_new_daughter_cell:
                    self.scaffold[cell.position[0]][cell.position[1]][cell.position[2]].cell_number -= 1
                cell.position[0] += 1
                self.scaffold[cell.position[0]][cell.position[1]][cell.position[2]].cell_number += 1
            elif migration_direction == 2:
                # X < 0
                if not is_new_daughter_cell:
                    self.scaffold[cell.position[0]][cell.position[1]][cell.position[2]].cell_number -= 1
                cell.position[0] -= 1
                self.scaffold[cell.position[0]][cell.position[1]][cell.position[2]].cell_number += 1
            elif migration_direction == 3:
                # Y > 0
                if not is_new_daughter_cell:
                    self.scaffold[cell.position[0]][cell.position[1]][cell.position[2]].cell_number -= 1
                cell.position[1] += 1
                self.scaffold[cell.position[0]][cell.position[1]][cell.position[2]].cell_number += 1
            elif migration_direction == 4:
                # Y < 0
                if not is_new_daughter_cell:
                    self.scaffold[cell.position[0]][cell.position[1]][cell.position[2]].cell_number -= 1
                cell.position[1] -= 1
                self.scaffold[cell.position[0]][cell.position[1]][cell.position[2]].cell_number += 1
            elif migration_direction == 5:
                # Z > 0
                if not is_new_daughter_cell:
                    self.scaffold[cell.position[0]][cell.position[1]][cell.position[2]].cell_number -= 1
                cell.position[2] += 1
                self.scaffold[cell.position[0]][cell.position[1]][cell.position[2]].cell_number += 1
            elif migration_direction == 6:
                # Z < 0
                if not is_new_daughter_cell:
                    self.scaffold[cell.position[0]][cell.position[1]][cell.position[2]].cell_number -= 1
                cell.position[2] -= 1
                self.scaffold[cell.position[0]][cell.position[1]][cell.position[2]].cell_number += 1


    def __calculate_cell_velocity(self, cell, migration_direction):
        """TODO"""

        # Calculate are and generate a normally distributed integrin and ligand distance
        integrin_ligand_separation = np.random.normal(0.04, 0.001, None)
        equilibrium_distance = 0.025
        area = np.pi * (0.06)**2

        # Calculate spring constant
        k = self.scaffold_stiffness * (area/equilibrium_distance)

        # Calculate force
        force = k * (integrin_ligand_separation - equilibrium_distance)

        # Calculate Tau
        tau = 1/(2.5*np.exp(-0.014 * 6.8 * force) + (0.5*10**(-5)) * np.exp(0.3 *0.96 * force)) 

        # --------------------------------------------------------------------------
        # Ligand Number
        # --------------------------------------------------------------------------
        # Calculate ligand number based on the pore the cell plans to migrate to
        x_increment = 0
        y_increment = 0
        z_increment = 0

        if migration_direction == 1:
            x_increment += 1
        elif migration_direction == 2:
            x_increment -= 1
        elif migration_direction == 3:
            y_increment += 1
        elif migration_direction == 4:
            y_increment -= 1
        elif migration_direction == 5:
            z_increment += 1
        elif migration_direction == 6:
            z_increment -= 1


        # *** Area calculations need to be checked and justified ***
        available_area = 0
        next_pore = self.scaffold[cell.position[0] + x_increment][cell.position[1] + y_increment][cell.position[2] + z_increment]
        if next_pore.surface_area < 1E3 * next_pore.calculate_cell_adhesion_area(self.cell_diameter):
            available_area = next_pore.surface_area/(1E3 * next_pore.calculate_cell_adhesion_area(self.cell_diameter))
        else:
            available_area = 1
        
        # Calculate ligand max and ligand number
        ligand_max = available_area*(np.pi*(self.cell_diameter/2)**2 - np.pi*(self.cell_diameter/2 - 4) ** 2) * 500
        ligand_num = self.ligand_factor*ligand_max*tau/(5.7499*10**3)

        # --------------------------------------------------------------------------
        # Calculate velocity
        # --------------------------------------------------------------------------
        c1 = 0.001                                  # TODO
        n = 55                                      # Viscosity (Pa-s)
        c = 6 * np.pi * self.cell_diameter/2        # Constant, dependent on cell shape (µm)

        return ((c1 * self.scaffold_stiffness * ligand_num)/(c * n))  # Cell velocity to planned pore to migrate to (µm/s)


    def neighbor_conditional_check(self, cell_to_check):
        """TODO"""

        # Generate cell conditions
        neighbor_conditions = [False, False, False, False, False, False, False]

        # Evaluate neighbor cell conditions as well as boundary conditions
        boundary_index = len(self.pore_array) - 1
        if (cell_to_check.position[0] + 1) > boundary_index:
            neighbor_conditions[0] = True
        else:
            neighbor_conditions[0] = self.scaffold[cell_to_check.position[0] + 1][cell_to_check.position[1]][cell_to_check.position[2]].is_full()

        if (cell_to_check.position[0] - 1) < 0:
            neighbor_conditions[1] = True
        else:
            neighbor_conditions[1] = self.scaffold[cell_to_check.position[0] - 1][cell_to_check.position[1]][cell_to_check.position[2]].is_full()

        if (cell_to_check.position[1] + 1) > boundary_index:
            neighbor_conditions[2] = True
        else:
            neighbor_conditions[2] = self.scaffold[cell_to_check.position[0]][cell_to_check.position[1] + 1][cell_to_check.position[2]].is_full()
            
        if (cell_to_check.position[1] - 1) < 0:
            neighbor_conditions[3] = True
        else:
            neighbor_conditions[3] = self.scaffold[cell_to_check.position[0]][cell_to_check.position[1] - 1][cell_to_check.position[2]].is_full()

        if (cell_to_check.position[2] + 1) > boundary_index:
            neighbor_conditions[4] = True
        else:
            neighbor_conditions[4] = self.scaffold[cell_to_check.position[0]][cell_to_check.position[1]][cell_to_check.position[2] + 1].is_full()

        if (cell_to_check.position[2] - 1) < 0:
            neighbor_conditions[5] = True
        else:
            neighbor_conditions[5] = self.scaffold[cell_to_check.position[0]][cell_to_check.position[1]][cell_to_check.position[2] - 1].is_full()
        neighbor_conditions[6] = self.scaffold[cell_to_check.position[0]][cell_to_check.position[1]][cell_to_check.position[2]].is_full()

        return neighbor_conditions


    def new_daughter_cell(self, parent_cell, time_in):
        """TODO"""

        daughter_cell = self.Cell(parent_cell.cell_diameter, time_in)
        daughter_cell.position = parent_cell.position[:]
        daughter_cell.ID = self.cell_ID_counter
        daughter_cell.generation = parent_cell.generation + 1

        return daughter_cell


    class Pore:
        def __init__(self, position, cell_number, max_cells, pore_diameter):
            """TODO"""

            self.position = position                                    # Position of pore within scaffold
            self.cell_number = cell_number                              # Number of cells currently occupying pore
            self.max_cells = max_cells                                  # Max number of cells that can occupy the pore
            self.surface_area = 4 * np.pi * (pore_diameter/2)**2        # Calculates the surface area of the pore (assuming spherical)


        def is_full(self):
            """TODO"""

            return self.cell_number >= self.max_cells


        def calculate_cell_adhesion_area(self, cell_diameter):
            """TODO"""

            return self.cell_number * np.pi * (cell_diameter/2) ** 2


    class Cell:
        def __init__(self, cell_diameter, time):
            """TODO"""

            self.position = [0, 0, 0]
            self.ID = 0
            self.generation = 1
            self.moves_made = 0
            self.time = time
            self.cell_diameter = cell_diameter


        # Generates data row regarding relevant cell parameters
        def write_data(self, pore_array, current_time):
            return [current_time, self.ID, self.generation, self.moves_made, round(pore_array[self.position[0]], 3), round(pore_array[self.position[1]], 3), round(pore_array[self.position[2]], 3)]