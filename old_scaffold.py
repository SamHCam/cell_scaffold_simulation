from operator import truediv
import numpy as np
import math as m
import random as r

class Scaffold:
    def __init__(self, dimension, porosity, pore_diameter, packing_density, cell_diameter, scaffold_stiffness,
                 ligand_factor):
        """Seven parameter constructor that defines a scaffold utilizing individual spherical pores.

        :param dimension: side dimension of cubical scaffold (µm).
        :param porosity: empty or void volume fraction of the scaffold (%).
        :param pore_diameter: diameter of pores occupying scaffold (µm).
        :param packing_density: fraction of void or empty volume within the scaffold that can be occupied by cells (%)
        :param cell_diameter: diameter of cells occupying scaffold (µm).
        :param scaffold_stiffness: the Young's modulus of the scaffold.
        :param ligand_factor: percentage of normal ligand percentage in the scaffold.
        """
        # --------------------------------------------------------------------------
        # Tracks the number of cells in the scaffold, cell IDs in the scaffold, and the time elapsed in the scaffold
        # --------------------------------------------------------------------------
        self.__cell_count = 0
        self.__cell_ID_counter = 1
        self.__time = 0

        # --------------------------------------------------------------------------
        # Generate scaffold characteristics
        # --------------------------------------------------------------------------
        # Generates various scaffold characteristics from method parameters
        cluster_array, pore_number, cluster_cell_max, scaffold_cell_max, pores_per_cluster = \
            self.__generate_pore_scaffold_properties(dimension, porosity, pore_diameter, packing_density, cell_diameter)


        # Member Variables
        self.__pore_number = pore_number                    # The number of pores in the scaffold
        self.__pore_diameter = pore_diameter                # The diameters of pores in the scaffold
        self.__cell_diameter = cell_diameter                # The cell diameter of cells that can be seeded in
                                                            # scaffold
        self.__cluster_array = cluster_array                # An array that represents the coordinates of each
                                                            # pore cluster along a scaffold side length to represent any
                                                            # position in the scaffold
        self.__cluster_cell_max = cluster_cell_max          # The maximum number of cells one pore cluster in the
                                                            # scaffold can hold
        self.__scaffold_cell_max = scaffold_cell_max        # The maximum number of cells the scaffold can hold
        self.__scaffold_stiffness = scaffold_stiffness      # The scaffold stiffness (Pa)
        self.__ligand_factor = ligand_factor                # The percent of normal ligand percentage

        # --------------------------------------------------------------------------
        # Generate the pores in the scaffold and assigns their positions
        # --------------------------------------------------------------------------
        # Number of clusters per scaffold side length
        cluster_array_size = len(cluster_array)

        # Generate the scaffold based on the number of clusters
        self.scaffold = [[[self.Pore_cluster(None, 0, self.__cluster_cell_max, pore_diameter, pores_per_cluster)
                           for x in range(cluster_array_size)]
                          for y in range(cluster_array_size)]
                         for z in range(cluster_array_size)]

        # Assign each pore cluster a unique position within the scaffold
        self.__assign_positions(cluster_array_size)

        # --------------------------------------------------------------------------
        # Initializes and stores cell objects that are in the scaffold
        # --------------------------------------------------------------------------
        self.cell_objs = []


    def __init__(self, dimension, porosity, pore_diameter, packing_density, cell_diameter, scaffold_stiffness,
                 ligand_factor, pore_cluster_count):
        """Eight parameter constructor that defines a scaffold utilizing clusters of spherical pores rather than
        individual pores.

        :param dimension: side dimension of cubical scaffold (µm).
        :param porosity: empty or void volume fraction of the scaffold (%).
        :param pore_diameter: diameter of pores occupying scaffold (µm).
        :param packing_density: fraction of void or empty volume within the scaffold that can be occupied by cells (%)
        :param cell_diameter: diameter of cells occupying scaffold (µm).
        :param scaffold_stiffness: the Young's modulus of the scaffold.
        :param ligand_factor: percentage of normal ligand percentage in the scaffold.
        :param pore_cluster_count: number of clusters that define all pores within the scaffold
        """
        # --------------------------------------------------------------------------
        # Tracks the number of cells in the scaffold, cell IDs in the scaffold
        # --------------------------------------------------------------------------
        self.__cell_count = 0
        self.__cell_ID_counter = 1

        # --------------------------------------------------------------------------
        # Generate scaffold characteristics
        # --------------------------------------------------------------------------
        # Generates various scaffold characteristics from method parameters
        cluster_array, pore_number, cluster_cell_max, scaffold_cell_max, pores_per_cluster = \
            self.__generate_environment_properties(dimension, porosity, pore_diameter, packing_density, cell_diameter,
                                                   pore_cluster_count)

        # Member Variables
        self.__pore_number = pore_number                    # The number of pores in the scaffold
        self.__cluster_count = pore_cluster_count           # The number of clusters in the scaffold
        self.__pore_diameter = pore_diameter                # The diameters of pores in the scaffold
        self.__cell_diameter = cell_diameter                # The cell diameter of cells that can be seeded in
                                                            # scaffold
        self.__cluster_array = cluster_array                # An array that represents the coordinates of each
                                                            # pore cluster along a scaffold side length to represent any
                                                            # position in the scaffold
        self.__cluster_cell_max = cluster_cell_max          # The maximum number of cells one pore cluster in the
                                                            # scaffold can hold
        self.__scaffold_cell_max = scaffold_cell_max        # The maximum number of cells the scaffold can hold
        self.__scaffold_stiffness = scaffold_stiffness      # The scaffold stiffness (Pa)
        self.__ligand_factor = ligand_factor                # The percent of normal ligand percentage

        # --------------------------------------------------------------------------
        # Generate the pores in the scaffold and assigns their positions
        # --------------------------------------------------------------------------
        # Number of clusters per scaffold side length
        cluster_array_size = len(cluster_array)

        # Generate the scaffold based on the number of clusters
        self.scaffold = [[[self.Pore_cluster(None, 0, self.__cluster_cell_max, pore_diameter, pores_per_cluster)
                           for x in range(cluster_array_size)]
                          for y in range(cluster_array_size)]
                         for z in range(cluster_array_size)]

        # Assign each pore cluster a unique position within the scaffold
        self.__assign_positions(cluster_array_size)

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

    def get_cluster_count(self):
        """Returns the number of pore clusters within the scaffold."""
        return self.__cluster_count

    def get_cell_diameter(self):
        """Returns the cell diameter (µm) of the cells defined in the scaffold."""
        return self.__cell_diameter

    def get_pore_size(self):
        """Returns the pore diameter (µm) of the cells defined in the scaffold."""
        return self.__pore_diameter

    def get_scaffold_stiffness(self):
        """Returns the scaffold's stiffness or Young's Modulus (Pa)."""
        return self.__scaffold_stiffness

    def get_time(self):
        """Returns the total time elapsed in the scaffold"""
        return self.__time

    # --------------------------------------------------------------------------
    # Private Class Methods
    # --------------------------------------------------------------------------
    def __generate_pore_scaffold_properties(self, dimension, porosity, pore_diameter, packing_density, cell_diameter):
        """Calculates scaffold properties related to the number of pores, pore distribution, and pore cell capacity

        :param dimension: side dimension of cubical scaffold (µm).
        :param porosity: empty or void volume fraction of the scaffold (%).
        :param pore_diameter: fraction of void or empty volume within the scaffold that can be occupied by cells (%).
        :param packing_density: fraction of void or empty volume within the scaffold that can be occupied by cells (%).
        :param cell_diameter: diameter of cells occupying scaffold (µm).
        :return: distribution of pores for the X, Y, and Z side dimensions, number of pores, maximum cells per
                 pore cluster, maximum cells in scaffold, pores per cluster.
        """
        # --------------------------------------------------------------------------
        # Calculates volume that is available to pores
        # --------------------------------------------------------------------------
        scaffold_volume = dimension**3                          # Scaffold volume, assuming cube (µm^3)
        scaffold_porous_volume = porosity * scaffold_volume     # Amount of void volume in scaffold (µm^3)
        pore_volume = 4/3*m.pi*(pore_diameter/2)**3             # Pore volume, assuming spherical (µm^3)

        # --------------------------------------------------------------------------
        # Calculates the number of pores based on the pore volume available
        # --------------------------------------------------------------------------
        pore_number = m.floor(scaffold_porous_volume/pore_volume)   # Calculates number of pores possible in scaffold
        pores_per_side = m.floor(pore_number ** (1/3))              # Re-calculates an appropriate number of pores
                                                                    # per side such that the scaffold is cubical
        pore_number = pores_per_side**3                             # Re-calculates new pore number to maintain
                                                                    # a cubical scaffold

        # --------------------------------------------------------------------------
        # Create side-length array
        # --------------------------------------------------------------------------
        pore_array = \
            np.linspace(0, dimension,
                        num=int(pores_per_side).tolist())       # Creates an array that represents
                                                                # the positions of clusters along
                                                                # one side-length of the scaffold.

        # Calculates the maximum cells in the scaffold and maximum cells in a single pore
        cell_volume = 4/3*np.pi*(cell_diameter/2)**3        # Calculates the volume of a cell in the scaffold
                                                            # assuming spherical (µm^3)
        occupied_pore_volume = \
            packing_density*pore_volume                     # Amount of volume in the pore that can be
                                                            # occupied by cells based on a packing density (µm^3)
        pore_cell_max =  \
            m.floor(occupied_pore_volume / cell_volume)  # Maximum number of cells per pore cluster
        max_cells = pore_cell_max * pore_number   # Maximum number of cells that can occupy the scaffold

        # Returns scaffold characteristics
        return pore_array, pore_number, pore_cell_max, pore_number, max_cells

    def __generate_cluster_scaffold_properties(self, dimension, porosity, pore_diameter, packing_density, cell_diameter,
                                          pore_cluster_count):
        """Calculates scaffold properties related to the number of pores, pore distribution, and pore cell capacity

        :param dimension: side dimension of cubical scaffold (µm).
        :param porosity: empty or void volume fraction of the scaffold (%).
        :param pore_diameter: fraction of void or empty volume within the scaffold that can be occupied by cells (%).
        :param packing_density: fraction of void or empty volume within the scaffold that can be occupied by cells (%).
        :param cell_diameter: diameter of cells occupying scaffold (µm).
        :param pore_cluster_count: number of clusters that define all pores within the scaffold.
        :return: distribution of pores for the X, Y, and Z side dimensions, number of pores, maximum cells per
                 pore cluster, maximum cells in scaffold, pores per cluster.
        """
        # --------------------------------------------------------------------------
        # Calculates volume that is available to pores
        # --------------------------------------------------------------------------
        scaffold_volume = dimension**3                          # Scaffold volume, assuming cube (µm^3)
        scaffold_porous_volume = porosity * scaffold_volume     # Amount of void volume in scaffold (µm^3)
        pore_volume = 4/3*m.pi*(pore_diameter/2)**3             # Pore volume, assuming spherical (µm^3)

        # --------------------------------------------------------------------------
        # Calculates the number of pores based on the pore volume available
        # --------------------------------------------------------------------------
        pore_number = m.floor(scaffold_porous_volume/pore_volume)   # Calculates number of pores possible in scaffold
        pores_per_side = m.floor(pore_number ** (1/3))              # Re-calculates an appropriate number of pores
                                                                    # per side such that the scaffold is cubical
        pore_number = pores_per_side**3                             # Re-calculates new pore number to maintain
                                                                    # a cubical scaffold
        pores_per_cluster = \
            m.floor(pore_number/pore_cluster_count)                 # Calculate the number of pores contained
                                                                    # in a single pore cluster

        # --------------------------------------------------------------------------
        # Create side-length array
        # --------------------------------------------------------------------------
        cluster_array = \
            np.linspace(0, dimension,
                        num=int(pore_cluster_count**(1/3))).tolist()    # Creates an array that represents
                                                                        # the positions of clusters along
                                                                        # one side-length of the scaffold.

        # Calculates the maximum cells in the scaffold and maximum cells in a single pore
        cell_volume = 4/3*np.pi*(cell_diameter/2)**3        # Calculates the volume of cells in the scaffold
                                                            # (assuming spherical) (µm^3)
        occupied_cluster_volume = \
            packing_density*pore_volume*pores_per_cluster   # Amount of volume in the pore that can be
                                                            # occupied by cells based on a packing density (µm^3)
        cluster_cell_max =  \
            m.floor(occupied_cluster_volume/cell_volume)  # Maximum number of cells per pore cluster
        max_cells = cluster_cell_max * pore_cluster_count   # Maximum number of cells that can occupy the scaffold

        # Returns scaffold characteristics
        return cluster_array, pore_number, cluster_cell_max, max_cells, pores_per_cluster


    def __assign_positions(self, array_size):
        """Assigns each of the clusters in the scaffold to their respective locations.

        :param array_size:An array that represents the coordinates of each pore cluster along a scaffold side length
        """
        p_x = 0
        p_y = 0
        p_z = 0
        while p_z < array_size:
            while p_y < array_size:
                while p_x < array_size:
                    self.scaffold[p_x][p_y][p_z].position = [p_x, p_y, p_z]
                    p_x += 1
                p_x = 0
                p_y += 1
            p_y = 0
            p_z += 1

# --------------------------------------------------------------------------
# Seeding Methods
# --------------------------------------------------------------------------
    def seed_cells_randomly(self, seeded_cell_count, data, seeded_time):
        """Randomly seed a specified number of cells throughout all pores in the scaffold."""

        # Throw if attempting to seed more cells than allowed in the scaffold
        if (seeded_cell_count > self.scaffold_cell_max):
            raise RuntimeError("The number of cells you are trying to seed cannot be larger than the maximum cell capacity of the scaffold.")

        # Generate seeded cell objects
        self.cell_objs = [self.Cell(self.cell_diameter, 0) for i in range(seeded_cell_count)]

        end_point = len(self.cluster_array) - 1
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

            # Write cell data
            data.append(cell.write_data(self.cluster_array, seeded_time))


    def seed_cells_at_top(self, seeded_cell_count, data, seeded_time):
        """Seeds cells at the center of the top layer of the scaffold"""
        # Generate seeded cell objects
        self.cell_objs = [self.Cell(self.cell_diameter, 0) for i in range(seeded_cell_count)]

        end_point = len(self.cluster_array) - 1
        for cell in self.cell_objs:
            # Assign a unique ID for each seeded cell
            cell.ID = self.cell_ID_counter
            self.cell_ID_counter += 1
            self.cell_count += 1

            # Random seed cells by generating random positions, writes its properties to the data list
            while True:
                new_position = [r.randint(int(end_point/2), int(end_point - end_point/2)), r.randint(int(end_point/2), int(end_point - end_point/2)), r.randint(int(end_point/2), int(end_point - end_point/2))]
                pore_destination = self.scaffold[new_position[0]][new_position[1]][new_position[2]]
                if (not pore_destination.is_full()):
                    cell.position = new_position
                    pore_destination.cell_number += 1
                    break


# --------------------------------------------------------------------------
# Public Methods
# --------------------------------------------------------------------------
    def scaffold_migration_replication(self, data, ts, tf, rf, rep_rate):
        """Performs a migration and replication cycle on all cells in the scaffold at a specified replication rate.

        :param data: data list to write cell data to
        :param ts:
        :param tf:
        :param rf:
        :param rep_rate:
        """
        # Ensure unbiased cellular migration and replication
        r.shuffle(self.cell_objs)

        # Migration and replication of parent and daughter cells
        for cell in self.cell_objs:

            # Neighbor capacity checks whether neighboring pores in regards to the pore the cell is in are full. X, Y, and Z in both positive and negative direction are accounted by conditions 1-6
            # while the pore's capacity, in which the current cell is in, is accounted by condition 7
            neighbor_capacity = self.__neighbor_conditional_check(cell)[:]

            # If all pores, including the pore the cell is in, are full then no movement and no replication occurs. Only recording of relevant cell properties to the data array is performed.
            if (neighbor_capacity[0] and neighbor_capacity[1] and neighbor_capacity[2] and neighbor_capacity[3] and neighbor_capacity[4] and neighbor_capacity[5] and neighbor_capacity[6]) == True:
                 # Writes cell properties to data list
                if (t % rf == 0 or t == tf):
                    data.append(cell.write_data(self.cluster_array, t))
            else:
                if (neighbor_capacity[0] and neighbor_capacity[1] and neighbor_capacity[2] and neighbor_capacity[3] and neighbor_capacity[4] and neighbor_capacity[5]) == True:
                    rp = 0
                    if rep_rate != 0:
                        if round(1/rep_rate) == 1:
                            rp = 1
                        else:
                            rp = r.randint(1, round(1/rep_rate))
                    if rp != 1:
                        if (t % rf == 0 or t == tf):
                            data.append(cell.write_data(self.cluster_array, t))
                    # Replication of cell occurs if replication condition is achieved
                    else:
                        new_daughter_cell = self.__new_daughter_cell(cell, t)

                        # New cell with an identical pore location to its parent is added to the list of cell objects
                        self.scaffold[cell.position[0]][cell.position[1]][cell.position[2]].cell_number += 1
                        self.cell_objs.append(new_daughter_cell)

                        # Update cell count and cell ID counter
                        self.cell_ID_counter += 1
                        self.cell_count += 1

                        # Writes cell properties to data list
                        if (t % rf == 0 or t == tf):
                            data.append(cell.write_data(self.cluster_array, t))
                            data.append(new_daughter_cell.write_data(self.cluster_array, t))
                else:
                    # Perform migration on parent cell
                    self.__perform_migration(cell, neighbor_capacity, ts, False)

                    # Perform replication which is dependent on the replication probability
                    rp = 0
                    if rep_rate != 0:
                        if round(1/rep_rate) == 1:
                            rp = 1
                        else:
                            rp = r.randint(1, round(1/rep_rate))
                    if rp != 1:
                        if (t % rf == 0 or t == tf):
                            data.append(cell.write_data(self.cluster_array, t))
                    else:
                        # Generate a new daughter cell and perform migration on it
                        new_daughter_cell = self.__new_daughter_cell(cell, t)
                        daughter_neighbor_conditions = self.__neighbor_conditional_check(new_daughter_cell)
                        self.__perform_migration(new_daughter_cell, daughter_neighbor_conditions, ts, True)

                        # Replicated cell must migrate to exist
                        if new_daughter_cell.moves_made != 0:
                            self.cell_objs.append(new_daughter_cell)

                            # Update cell count and cell ID counter
                            self.cell_ID_counter += 1
                            self.cell_count += 1

                        # Writes cell properties to data list
                        if (t % rf == 0 or t == tf):
                            data.append(cell.write_data(self.cluster_array, t))
                            if new_daughter_cell.moves_made != 0:
                                data.append(new_daughter_cell.write_data(self.cluster_array, t))


    def __perform_migration(self, cell, neighbor_capacity, ts, is_new_daughter_cell):
        """Performs cell migration on the specified cell"""
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
        f_trac, cell_velocity = self.__calculate_cell_velocity(cell, migration_direction)           # Cell velocity in chosen direction (µm/s)
        distance_between_clusters = self.cluster_array[1] - self.cluster_array[0]                      # Distance between neighboring pores (equidistant)
        cell_distance_traveled = cell_velocity * ts * 3600                                          # Distance cell can travel (µm)

        # Calculate probability to migrate as a ratio between distance cell can travel and distance between pores
        probability_of_migration = cell_distance_traveled/distance_between_clusters * 100      # Probability of migration (%)

        # --------------------------------------------------------------------------
        # Cell Migration
        # --------------------------------------------------------------------------
        # Generates a number between 1 and 100 to determine if cell will migrate
        random_migration_number = r.randint(0, 100000)/1000

        # Migrate in chosen direction based on migration probability
        if probability_of_migration > random_migration_number:
            # Update cell's traction force
            cell.f_trac = f_trac

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
        area = np.pi * (0.06)**2


        stiffness_to_use = self.scaffold_stiffness
        if self.scaffold_stiffness >= 1E6:
            stiffness_to_use = 1E3

        # Calculate spring constant
        k = stiffness_to_use * (area/equilibrium_distance)

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
        next_pore = self.scaffold[cell.position[0] + x_increment][cell.position[1] + y_increment][cell.position[2] + z_increment]
        if next_pore.surface_area < 1E3 * next_pore.calculate_cell_adhesion_area(self.cell_diameter):
            available_area = next_pore.surface_area/(1E3 * next_pore.calculate_cell_adhesion_area(self.cell_diameter))
        else:
            available_area = 1

        # Calculate ligand max and ligand number
        ligand_max = available_area*(np.pi*(self.cell_diameter/2)**2 - np.pi*(self.cell_diameter/2 - 8) ** 2) * 500
        ligand_num = self.ligand_factor*ligand_max*tau/(1E6)

        # --------------------------------------------------------------------------
        # Calculate velocity
        # --------------------------------------------------------------------------
        c1 = 0.001                                  # TODO
        n = 55                                      # Viscosity (Pa-s)
        c = 6 * np.pi * self.cell_diameter/2        # Constant, dependent on cell shape (µm)

        f_traction = c1 * stiffness_to_use * ligand_num

        return f_traction, 0.005#(f_traction/(c * n))  # Cell velocity to planned pore to migrate to (µm/s)


    def __neighbor_conditional_check(self, cell_to_check):
        """Checks whether the cell's current pore's neighboring cells are full. Returns a list of booleans indicating which pores are full."""

        # Generate cell conditions
        neighbor_conditions = [False, False, False, False, False, False, False]

        # Evaluate neighbor cell conditions as well as boundary conditions
        boundary_index = len(self.cluster_array) - 1
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


    def __new_daughter_cell(self, parent_cell, time_in):
        """Generate a new daughter cell based on the parent cell. Returns the new daughter cell object"""

        daughter_cell = self.Cell(parent_cell.cell_diameter, time_in)
        daughter_cell.position = parent_cell.position[:]
        daughter_cell.ID = self.cell_ID_counter
        daughter_cell.generation = parent_cell.generation + 1

        return daughter_cell


    class Pore:
        def __init__(self, position, cell_number, max_cells, pore_diameter):
            """Four parameter pore constructor that takes in its position in the scaffold, the number of cells currently occupying the cell, the maximum
            cells it can hold, and its pore diameter."""

            self.position = position                                    # Position of pore within scaffold
            self.cell_number = cell_number                              # Number of cells currently occupying pore
            self.max_cells = max_cells                                  # Max number of cells that can occupy the pore
            self.surface_area = 4 * np.pi * (pore_diameter/2)**2        # Calculates the surface area of the pore (assuming spherical pores)


        def is_full(self):
            """Checks whether the pore is currently full. Returns true if it is."""

            return self.cell_number >= self.max_cells


        def calculate_cell_adhesion_area(self, cell_diameter):
            """Calculates the current area of the pore taken up by cells"""

            return self.cell_number * np.pi * (cell_diameter/2) ** 2

    class Pore_cluster:
        def __init__(self, position, cell_number, max_cells, pore_diameter, pores_per_cluster):
            """Four parameter constructor representing a cluster of pores in the scaffold."""

            self.position = position                                                        # Position of pore within scaffold
            self.cell_number = cell_number                                                  # Number of cells currently occupying pore
            self.max_cells = max_cells                                                      # Max number of cells that can occupy the pore
            self.surface_area = 4 * np.pi * (pore_diameter/2)**2 * pores_per_cluster        # Calculates the total surface area of the pore cluster (assuming spherical pores)


        def is_full(self):
            """Checks whether the pore cluster is currently full. Returns true if it is."""

            return self.cell_number >= self.max_cells


        def calculate_cell_adhesion_area(self, cell_diameter):
            """Calculates the current area of the pore cluster taken up by cells"""

            return self.cell_number * np.pi * (cell_diameter/2) ** 2



    class Cell:
        def __init__(self, cell_diameter, time):
            """Two parameter cell constructor that takes in the cell's diameter and its time of entry into the simulation."""

            self.position = [0, 0, 0]
            self.ID = 0
            self.generation = 1
            self.moves_made = 0
            self.time = time
            self.cell_diameter = cell_diameter
            self.f_trac = 0


        # Generates data row regarding relevant cell parameters
        def write_data(self, cluster_array, current_time):
            """Returns a data array detailing various properties of the cell including the current simulation time, its unique cell ID, the number of moves it has made,
            the force of traction, and its X, Y, and Z position in the scaffold."""
            # return [current_time, self.ID, self.generation, self.moves_made, self.f_trac, self.position[0], self.position[1], self.position[2]]
            return [current_time, self.ID, self.generation, self.moves_made, self.f_trac, round(cluster_array[self.position[0]], 3), round(cluster_array[self.position[1]], 3), round(cluster_array[self.position[2]], 3)]