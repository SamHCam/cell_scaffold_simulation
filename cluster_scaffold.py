import numpy as np
import math as m
import random as r

class Cluster_Scaffold:
    def __init__(self, dimension, porosity, pore_diameter, packing_density, cell_diameter, scaffold_stiffness,
                 ligand_factor, pore_column_layers):
        """Seven parameter constructor that defines a scaffold utilizing individual spherical pores.

        :param dimension: side dimension of cubical scaffold (µm).
        :param porosity: empty or void volume fraction of the scaffold (%).
        :param pore_diameter: diameter of pores occupying scaffold (µm).
        :param packing_density: fraction of void or empty volume within the scaffold that can be occupied by cells (%)
        :param cell_diameter: diameter of cells occupying scaffold (µm).
        :param scaffold_stiffness: the Young's modulus of the scaffold.
        :param ligand_factor: percentage of normal ligand percentage in the scaffold.
        :param pore_column_layers: number of layers to represent a single pore column in the scaffold
        """
        # --------------------------------------------------------------------------
        # Tracks the number of cells in the scaffold, cell IDs in the scaffold, and time elapsed in the scaffold
        # --------------------------------------------------------------------------
        self.__cell_count = 0
        self.__cell_ID_counter = 1
        self.__time = 0

        # --------------------------------------------------------------------------
        # Generate scaffold characteristics
        # --------------------------------------------------------------------------
        # Generates various scaffold characteristics from method parameters
        pore_array, pore_number, pore_cell_max, scaffold_cell_max = \
            self.__generate_pore_scaffold_properties(dimension, porosity, pore_diameter, packing_density, cell_diameter,
                                                     pore_column_layers)

        # Scaffold property member variables
        self.__pore_column_layers = pore_column_layers      # The number of pores in the scaffold
        self.__pore_diameter = pore_diameter                # The diameters of pores in the scaffold
        self.__cell_diameter = cell_diameter                # The cell diameter of cells that can be seeded in
                                                            # scaffold
        self.__pore_array = pore_array                      # An array that represents the coordinates of each
                                                            # pore along a scaffold side length to represent any
                                                            # position in the scaffold
        self.__pore_cell_max = pore_cell_max                # The maximum number of cells one pore segment in the
                                                            # scaffold can hold
        self.__scaffold_cell_max = scaffold_cell_max        # The maximum number of cells the scaffold can hold
        self.__scaffold_stiffness = scaffold_stiffness      # The scaffold stiffness (Pa)
        self.__ligand_factor = ligand_factor                # The percent of normal ligand percentage

        # --------------------------------------------------------------------------
        # Generate the pores in the scaffold and assigns their positions
        # --------------------------------------------------------------------------
        # Number of clusters per scaffold side length
        cluster_array_size = len(pore_array)

        # Generate the scaffold based on the number of clusters
        self.scaffold = [[[self.Pore([x, y, z], self.__pore_cell_max, pore_diameter)
                           for x in range(cluster_array_size)]
                          for y in range(cluster_array_size)]
                         for z in range(pore_column_layers)]

        # --------------------------------------------------------------------------
        # Initializes and stores cell objects that are in the scaffold
        # --------------------------------------------------------------------------
        self.cell_objs = []

    # --------------------------------------------------------------------------
    # Getters and Setters for Scaffold
    # --------------------------------------------------------------------------
    def get_cell_count(self):
        """Returns cell count in the scaffold"""
        return self.__cell_count

    def get_pore_number(self):
        """Returns the number of pores in the scaffold."""
        return self.__pore_number

    def get_cell_diameter(self):
        """Returns the cell diameter (µm) of the cells defined in the scaffold."""
        return self.__cell_diameter

    def get_scaffold_cell_max(self):
        """Returns the maximum number of cells that the scaffold can hold"""
        return self.__scaffold_cell_max

    def get_pore_size(self):
        """Returns the pore diameter (µm) of pore columns defined in the scaffold."""
        return self.__pore_diameter

    def get_scaffold_stiffness(self):
        """Returns the scaffold's stiffness or Young's Modulus (Pa)."""
        return self.__scaffold_stiffness

    # --------------------------------------------------------------------------
    # Private Class Methods
    # --------------------------------------------------------------------------
    @staticmethod
    def __generate_pore_scaffold_properties(dimension, porosity, pore_diameter, packing_density, cell_diameter,
                                            pore_column_layers):
        """Calculates scaffold properties related to the number of pores, pore distribution, and pore cell capacity.

        :param dimension: side dimension of cubical scaffold (µm).
        :param porosity: empty or void volume fraction of the scaffold (%).
        :param pore_diameter: fraction of void or empty volume within the scaffold that can be occupied by cells (%).
        :param packing_density: fraction of void or empty volume within the scaffold that can be occupied by cells (%).
        :param cell_diameter: diameter of cells occupying scaffold (µm).
        :return: distribution of pores for the X, Y, and Z side dimensions, number of pores, maximum cells per
                 pore, maximum cells in scaffold
        """
        # --------------------------------------------------------------------------
        # Calculates volume that is available to pores
        # --------------------------------------------------------------------------
        scaffold_volume = dimension**3                          # Scaffold volume, assuming cube (µm^3)
        scaffold_porous_volume = porosity * scaffold_volume     # Amount of void volume in scaffold (µm^3)

        # --------------------------------------------------------------------------
        # Calculates number of pore columns and column distribution in the scaffold
        # --------------------------------------------------------------------------
        pore_column_cs_area = np.pi * (pore_diameter / 2) ** 2             # Calculates the pore-column cross-sectional
                                                                           # area (µm^2)
        pore_column_volume = pore_column_cs_area * dimension               # Volume of one porous column (µm^3)
        pore_column_count = scaffold_porous_volume / pore_column_volume    # Total potential pore columns in scaffold
        pore_columns_per_side = m.floor(pore_column_count ** (1/3))        # Calculates the number of pores per side and
                                                                           # maintains scaffold's cubical nature

        pore_columns = pore_columns_per_side ** 3   # Recalculate new pore number to maintain a cubical scaffold

        # --------------------------------------------------------------------------
        # Create side-length array that represents the number of columns per scaffold side-length
        # --------------------------------------------------------------------------
        pore_array = \
            np.linspace(0, dimension,
                        num=int(pore_columns_per_side)).tolist()

        # --------------------------------------------------------------------------
        # Calculate the volume of pore_column segments
        # --------------------------------------------------------------------------

        # Calculate the volumes of each pore segment
        pore_segment_volume = pore_column_volume / pore_column_layers   # Volume of one segment of the pore column
        cell_volume = 4 / 3 * np.pi * (cell_diameter / 2) ** 3  # Calculates the volume of a cell in the scaffold
                                                                # assuming spherical (µm^3)

        # Calculate the maximum number of cells in one pore segment
        pore_cell_max = m.floor(pore_segment_volume * packing_density / cell_volume)

        if pore_cell_max < 1:
            raise Exception("Current pore segment volume is too small to hold at least one cell.")

        # Calculate maximum number of cells that can occupy the scaffold
        max_cells = pore_cell_max * pore_columns * pore_column_layers

        # Returns scaffold characteristics
        return pore_array, pore_columns, pore_cell_max, max_cells

# --------------------------------------------------------------------------
# Seeding Methods
# --------------------------------------------------------------------------
    def seed_cells_at_top(self, seeded_cell_count, data, seeded_time):
        """Seeds cells at the center of the top layer of the scaffold"""
        # Generate seeded cell objects
        self.cell_objs = [self.Cell(self.cell_diameter, seeded_time) for i in range(seeded_cell_count)]

        # Find end index and middle index
        end_point = len(self.cluster_array) - 1
        middle_point = end_point/2

        # Seed each cell object
        for cell in self.cell_objs:
            # Assign a unique ID for each seeded cell
            cell.ID = self.cell_ID_counter
            self.cell_ID_counter += 1
            self.cell_count += 1

            # Random seed cells by generating random positions, writes its properties to the data list
            while True:
                new_position = [r.randint(int(middle_point - end_point/2), int(middle_point + end_point/2)),
                                r.randint(int(middle_point - end_point/2), int(middle_point + end_point/2)),
                                r.randint(int(end_point - 2), int(end_point))]
                pore_destination = self.scaffold[new_position[0]][new_position[1]][new_position[2]]
                if (not pore_destination.is_full()):
                    cell.x = new_position[0]
                    cell.y = new_position[1]
                    cell.z = new_position[2]
                    pore_destination.cell_number += 1
                    break

# --------------------------------------------------------------------------
# Public Methods
# --------------------------------------------------------------------------
    def migration_replication(self, data, ts, tf, rf, rep_rate):
        """Performs a migration and replication cycle on all cells in the scaffold at a specified replication rate.
        :param data: Stores data in rows in an expandable list
        :param ts: specifies the simulation's time-step (hour)
        :param tf:
        :param rf:
        :param rep_rate:
        """

        # Get current time elapsed
        t = self.__time

        # Ensure unbiased cellular migration and replication
        r.shuffle(self.cell_objs)

        # Migration and replication of parent and daughter cells
        for cell in self.cell_objs:

            # Get current status of neighboring pore capacity
            n_cap = self.__neighbor_conditional_check(cell)[:]

            # If all pores, including the pore the cell is in, are full then no movement and no replication occurs
            if n_cap[0] and n_cap[1] and n_cap[2] and n_cap[3] and n_cap[4] and n_cap[5] and n_cap[6]:
                # Writes cell properties to data list
                if t % rf == 0 or t == tf:
                    data.append(cell.write_data(self.__pore_array, t))
            else:
                # Only the pore the cell is currently in is not full
                if n_cap[0] and n_cap[1] and n_cap[2] and n_cap[3] and n_cap[4] and n_cap[5]:
                    # Generate a random replicate rate number based on specified replication rate
                    rp = self.generate_replication_rate_number(rep_rate)

                    # Replicate only if replication number matches
                    if rp != 1:
                        # Writes cell properties to data list
                        if t % rf == 0 or t == tf:
                            data.append(cell.write_data(self.__pore_array, t))
                    # Replication of cell occurs if replication condition is achieved
                    else:
                        # Create a new daughter cell
                        new_daughter_cell = self.__new_daughter_cell(cell, t)

                        # New cell with an identical pore location to its parent is added to the list of cell objects
                        self.scaffold[cell.x][cell.y][cell.z].cell_number += 1
                        self.cell_objs.append(new_daughter_cell)

                        # Increment cell ID and cell count
                        self.increment_cell_tracking()

                        # Writes cell properties to data list
                        if t % rf == 0 or t == tf:
                            data.append(cell.write_data(self.__pore_array, t))
                            data.append(new_daughter_cell.write_data(self.__pore_array, t))

                # Neighboring pores are not all full
                else:
                    # Perform migration on parent cell
                    self.__perform_migration(cell, n_cap, ts, False)

                    # Generate a random replicate rate number based on specified replication rate
                    rp = self.generate_replication_rate_number(rep_rate)

                    # Replicate only if replication number matches
                    if rp != 1:
                        # Writes cell properties to data list
                        if t % rf == 0 or t == tf:
                            data.append(cell.write_data(self.__pore_array, t))
                    else:
                        # Generate a new daughter cell and perform migration on it
                        new_daughter_cell = self.__new_daughter_cell(cell, t)
                        daughter_neighbor_conditions = self.__neighbor_conditional_check(new_daughter_cell)
                        self.__perform_migration(new_daughter_cell, daughter_neighbor_conditions, ts, True)

                        # Replicated cell must migrate to exist
                        if new_daughter_cell.moves_made != 0:
                            self.cell_objs.append(new_daughter_cell)

                            # Increment cell ID and cell count
                            self.increment_cell_tracking()

                        # Writes cell properties to data list
                        if t % rf == 0 or t == tf:
                            data.append(cell.write_data(self.__cluster_array, t))
                            # Write daughter cell properties to data list only if it exists
                            if new_daughter_cell.moves_made != 0:
                                data.append(new_daughter_cell.write_data(self.__cluster_array, t))

    def increment_cell_tracking(self):
        """Update cell count and cell ID counter"""
        self.cell_ID_counter += 1
        self.cell_count += 1

    @staticmethod
    def __generate_replication_rate_number(rep_rate):
        """generates a replication rate number that determines if replication occurs

        :param rep_rate: replication rate (s^-1)
        :return: replicate rate number
        """
        rp = 0
        # Replicate rate is not set to zero
        if rep_rate != 0:
            # Migration rate is capped
            if round(1 / rep_rate) == 1:
                print("Warning: replication rate is capped, consider increasing your time step.")
                rp = 1
            else:
                # Generate a random replication number
                rp = r.randint(1, round(1 / rep_rate))
        return rp

    def __perform_migration(self, cell, neighbor_capacity, ts, is_new_daughter_cell):
        """Performs cell migration on the specified cell"""
        # --------------------------------------------------------------------------
        # Choose a random direction for migration
        # --------------------------------------------------------------------------
        migration_direction = None
        # Ensure direction is not full or out of bounds
        while True:
            migration_direction = r.randint(1, 6)
            if cell.z < 0:
                if (neighbor_capacity[0] == True) and (migration_direction == 1):
                    continue
                elif (neighbor_capacity[1] == True) and (migration_direction == 2):
                    continue
                elif (neighbor_capacity[2] == True) and (migration_direction == 3):
                    continue
                elif (neighbor_capacity[3] == True) and (migration_direction == 4):
                    continue
            else:
                if (neighbor_capacity[4] == True) and (migration_direction == 5):
                    continue
                elif (neighbor_capacity[5] == True) and (migration_direction == 6):
                    continue
                else:
                    break

        # --------------------------------------------------------------------------
        # Calculate probability to migrate in chosen direction
        # --------------------------------------------------------------------------
        # Calculate distance
        cell_velocity = self.__calculate_cell_velocity(cell, migration_direction)       # Cell velocity in chosen direction (µm/s)
        distance_between_clusters = self.__pore_array[1] - self.__pore_array[0]                 # Distance between neighboring pores (µm, equidistant)
        cell_distance_traveled = cell_velocity * ts * 3600                                      # Distance cell can travel (µm)

        # Calculate probability to migrate as a ratio between distance cell can travel and distance between pores
        probability_of_migration = cell_distance_traveled/distance_between_clusters * 100       # Probability of migration (%)

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
        """Calculates the cell velocity (speed) towards one of six non-full neighboring pores. Returns the cell velocity."""

        # Calculate area and generate a normally distributed integrin and ligand distance
        integrin_ligand_separation = np.random.normal(0.04, 0.001, None)
        equilibrium_distance = 0.025
        area = np.pi * 0.06 ** 2

        # Calculate spring constant
        k = self.__scaffold_stiffness * (area/equilibrium_distance)

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


        # Calculate available surface area in a given pore based on the current capacity of the pore
        available_area = 0
        next_pore = \
            self.scaffold[cell.position[0] + x_increment][cell.position[1] + y_increment][cell.position[2] + z_increment]
        if next_pore.surface_area < 1E3 * next_pore.calculate_cell_adhesion_area(self.__cell_diameter):
            available_area = next_pore.surface_area/(1E3 * next_pore.calculate_cell_adhesion_area(self.__cell_diameter))
        else:
            available_area = 1

        # Calculate ligand max and ligand number
        ligand_max = available_area*(np.pi*(self.cell_diameter/2)**2 - np.pi*(self.cell_diameter/2 - 8) ** 2) * 500
        ligand_num = self.ligand_factor*ligand_max*tau/1E6

        # --------------------------------------------------------------------------
        # Calculate velocity
        # --------------------------------------------------------------------------
        c1 = 0.001                                  # Cell shape factor
        n = 55                                      # Viscosity (Pa-s)
        c = 6 * np.pi * self.cell_diameter/2        # Constant, dependent on cell shape (µm)

        f_traction = c1 * self.__scaffold_stiffness * ligand_num

        return f_traction/(c * n)  # Cell velocity to planned pore to migrate to (µm/s)


    def __neighbor_conditional_check(self, cell_to_check):
        """Checks whether the cell's current pore's neighboring cells are full. Returns a list of booleans indicating which pores are full."""

        # Generate cell conditions
        neighbor_conditions = [False, False, False, False, False, False, False]

        # Evaluate neighbor cell conditions as well as boundary conditions
        x_y_bounds = len(self.cluster_array) - 1
        z_bounds = self.__pore_column_layers - 1
        if (cell_to_check.x + 1) > x_y_bounds:
            neighbor_conditions[0] = True
        else:
            neighbor_conditions[0] = self.scaffold[cell_to_check.x + 1][cell_to_check.y][cell_to_check.z].is_full()

        if (cell_to_check.x - 1) < 0:
            neighbor_conditions[1] = True
        else:
            neighbor_conditions[1] = self.scaffold[cell_to_check.x - 1][cell_to_check.y][cell_to_check.z].is_full()

        if (cell_to_check.y + 1) > x_y_bounds:
            neighbor_conditions[2] = True
        else:
            neighbor_conditions[2] = self.scaffold[cell_to_check.x][cell_to_check.y + 1][cell_to_check.z].is_full()

        if (cell_to_check.y - 1) < 0:
            neighbor_conditions[3] = True
        else:
            neighbor_conditions[3] = self.scaffold[cell_to_check.x][cell_to_check.y - 1][cell_to_check.z].is_full()

        if (cell_to_check.y + 1) > z_bounds:
            neighbor_conditions[4] = True
        else:
            neighbor_conditions[4] = self.scaffold[cell_to_check.x][cell_to_check.y][cell_to_check.z + 1].is_full()

        if (cell_to_check.z - 1) < 0:
            neighbor_conditions[5] = True
        else:
            neighbor_conditions[5] = self.scaffold[cell_to_check.x][cell_to_check.y][cell_to_check.z - 1].is_full()
        # Check if pore cell is currently in is full
        neighbor_conditions[6] = self.scaffold[cell_to_check.x][cell_to_check.y][cell_to_check.z].is_full()

        return neighbor_conditions


    def __new_daughter_cell(self, parent_cell, time_in):
        """Generate a new daughter cell based on the parent cell. Returns the new daughter cell object"""

        daughter_cell = self.Cell(parent_cell.cell_diameter, time_in)
        daughter_cell.position = parent_cell.position[:]
        daughter_cell.ID = self.__cell_ID_counter
        daughter_cell.generation = parent_cell.generation + 1

        return daughter_cell


    class Pore:
        def __init__(self, position, max_cells, pore_diameter):
            """Four parameter constructor for a pore object
            :param position:
            :param cell_number:
            :param max_cells:
            :param pore_diameter:
            """

            self.x = position[0]                                        # X-position of pore within scaffold
            self.y = position[1]                                        # Y-position of pore within scaffold
            self.z = position[2]                                        # Z-position of pore within scaffold
            self.cell_number = 0                                        # Number of cells currently occupying pore
            self.max_cells = max_cells                                  # Max number of cells that can occupy the pore
            self.surface_area = 4 * np.pi * (pore_diameter/2)**2        # Calculates the surface area of the pore (assuming spherical pores)


        def is_full(self):
            """Checks whether the pore is currently full. Returns true if it is."""

            return self.cell_number >= self.max_cells


        def calculate_cell_adhesion_area(self, cell_diameter):
            """Calculates the current area of the pore taken up by cells"""

            return self.cell_number * np.pi * (cell_diameter/2) ** 2

    class Cell:
        def __init__(self, cell_diameter, time_in):
            """Two parameter cell constructor that takes in the cell's diameter and its time of entry into the simulation."""

            self.x = 0
            self.y = 0
            self.z = 0
            self.ID = 0
            self.generation = 1
            self.moves_made = 0
            self.time_in = time_in
            self.cell_diameter = cell_diameter
            self.f_trac = 0


        # Generates data row regarding relevant cell parameters
        def write_data(self, cluster_array, current_time):
            """Returns a data array detailing various properties of the cell including the current simulation time, its unique cell ID, the number of moves it has made,
            the force of traction, and its X, Y, and Z position in the scaffold."""
            # return [current_time, self.ID, self.generation, self.moves_made, self.f_trac, self.position[0], self.position[1], self.position[2]]
            return [current_time, self.ID, self.generation, self.moves_made, self.f_trac, round(cluster_array[self.position[0]], 3), round(cluster_array[self.position[1]], 3), round(cluster_array[self.position[2]], 3)]