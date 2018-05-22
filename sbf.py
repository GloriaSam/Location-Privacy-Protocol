"""
    Spatial Bloom Filter Python Library (libSBF-python)
    Copyright (C) 2017  Luca Calderoni, Dario Maio,
                        University of Bologna
    Copyright (C) 2017  Paolo Palmieri,
                        University College Cork

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
 """


import hashlib
import numpy as np
import base64
import csv
from sys import byteorder
from datetime import datetime


class sbf:
    """ Implements the Spatial Bloom Filter data structure.

    This class implements the Spatial Bloom Filter (SBF) data structure, which
    is an extention of the classical Bloom filter concept. An SBF is a compact,
    hash-based data structure representing a number of disjoint sets (called
    areas). SBFs can be queried to ascertain membership of an element in a set,
    with a known false-positive probability, and a 0 false-negative probability.
    For more informatio and the relevant scientific literature, refer to the
    project website: http://sbf.csr.unibo.it

    Attributes:
        bit_mapping     actual size of the filter (as in number of cells): for
                        instance, bit_mapping = 10 states that the SBF will be
                        composed by 2^10 cells. As such, the size can only be
                        defined as a power of 2. This value is bounded by the
                        MAX_BIT_MAPPING constant.
        hash_family     specifies the hash function to be used. Currently
                        available: MD4, MD5, SHA1
        num_hash        number of digests to be produced (by running the hash
                        function specified above using different salts) in the
                        insertion and check phases.
        num_areas       number of areas over which to build the filter.
        hash_salt_path  path to the file where to read from/write to the hash
                        salts (which are used to produce different digests).
                        If the file exists, reads one salt per line.
                        If the file doesn't exist, the salts are randomly
                        generated during the filter creation phase
    """

    # The maximum string (as a char array) length in bytes of each element'
    # given as input to be mapped in the SBF'
    MAX_INPUT_SIZE = 128
    # This value defines the maximum size (as in number of cells) of the SBF:
    # MAX_BIT_MAPPING = 32 states that the SBF will be composed at most by
    # 2^32 cells. The value is the number of bits used for SBF indexing.
    MAX_BIT_MAPPING = 32
    # Utility byte value of the above MAX_BIT_MAPPING
    MAX_BYTE_MAPPING = MAX_BIT_MAPPING/8
    # The maximum number of allowed areas. This way, we limit the memory size'
    # (which is the memory size of each cell) to 2 bytes'
    MAX_AREA_NUMBER = 65535
    # The maximum number of allowed digests
    MAX_HASH_NUMBER = 1024
    # The available hash families
    HASH_FAMILIES = ['md4', 'MD4', 'md5', 'MD5', 'sha', 'SHA', 'sha1', 'SHA1']


    def __init__(self, bit_mapping, hash_family, num_hashes, num_areas, hash_salt_path):
        """Inits the SBF class.

        The constructor performs argument validation, class variables
        initialization and loads or create the hash salts file.

        Raises:
            AttributeError: The arguments are out of bounds.
        """

        self.bit_mapping = bit_mapping
        self.hash_family = hash_family
        self.num_hashes = num_hashes
        self.num_areas = num_areas
        self.hash_salt_path = hash_salt_path

        # Argument validation
        if (self.bit_mapping <= 0) or (self.bit_mapping > self.MAX_BIT_MAPPING):
            raise AttributeError("Invalid bit mapping.")
        if (not self.hash_family in self.HASH_FAMILIES):
            raise AttributeError("Invalid hash family.")
        if (self.num_hashes <= 0) or (self.num_hashes > self.MAX_HASH_NUMBER):
            raise AttributeError("Invalid number of hash runs.")
        if (self.num_areas <= 0) or (self.num_areas > self.MAX_AREA_NUMBER):
            raise AttributeError("Invalid number of areas.")
        if (self.hash_salt_path == ''):
            raise AttributeError("Invalid hash salt path.")

        # Defines the number of bytes required for each cell depending on AREA_number
        # In order to reduce the memory footprint of the filter, we use 1 byte
        # for a number of areas <= 255, 2 bytes for a up to MAX_AREA_NUMBER
        self.cell_size = 0
        if (self.num_areas <= 255):
            self.cell_size = 1
        elif (self.num_areas > 255):
            self.cell_size = 2

        # Initializes the hash_salts list
        self.hash_salts = []
        # Tries to load the hash salts from the specified file
        try:
            with open(self.hash_salt_path) as self.salt_file:
                self.hash_salts = self.load_hash_salt(self.salt_file)
        except IOError:
            # File couldn't be opened, thus creates it as well as the salts
            with open(self.hash_salt_path, 'w') as self.salt_file:
                self.hash_salts = self.create_hash_salt(self.salt_file)

        # Defines the number of cells in the filter
        self.num_cells = pow(2,self.bit_mapping)

        # Defines the total size in bytes of the filter
        #self.size = self.cell_size*self.num_cells

        # Allocates the right amount of memory for the array by specifying the dtype,
        # and initializes the cells to 0
        if (self.cell_size == 1):
            self.filter = np.array(np.repeat(0,self.num_cells), dtype=np.uint8)
        elif (self.cell_size == 2):
            self.filter = np.array(np.repeat(0,self.num_cells), dtype=np.uint16)

        # Filter statistics initialization:
        # number of elements in the filter
        self.members = 0
        # number of collisions in the filter
        self.collisions = 0
        # number of members for each area
        self.area_members = [0]*(self.num_areas + 1)
        # number of cells for each area
        self.area_cells = [0]*(self.num_areas + 1)
        # number of self collisions for each area
        self.area_self_collisions = [0]*(self.num_areas + 1)
        # false positive probability for each area
        # (initialized in the method compute_area_fpp)
        #self.area_fpp = [0]*(self.num_areas + 1)
        # list of file from which elements have been inserted
        self.insert_file_list = []


    def __del__(self):
        """Destructs the SBF object"""
        pass


    def __enter__(self):
        """ Enters the SBF class (to be used with the 'with' statement)."""
        try:
            return self
        except StopIteration:
            raise RuntimeError("Error during objct generation")


    def __exit__(self, type, value, traceback):
        """ Exits the SBF class (to be used with the 'with' statement)."""
        pass


    def create_hash_salt(self, salt_file):
        """ Creates and saves to disk the hash salts.

        Creates a hash salt for each hash (the number of hashes is num_hashes).
        Each input element will be combined with the salt via XOR, by the insert
        and check methods. The length of salts is MAX_INPUT_SIZE bytes.
        Hashes are stored encoded in base64.

        Args:
            salt_file:  the file where to store the hash salts.

        Returns:
            The hash salts (in a list).
        """

        self.salt_file = salt_file
        self.hash_salts = []

        for self.i in range(0, self.num_hashes):
            self.hash_salts.append(np.random.bytes(self.MAX_INPUT_SIZE))

        for self.i in self.hash_salts:
            self.salt_file.write(base64.b64encode(self.i).decode('UTF-8') + "\n")

        self.salt_file.close()
        del self.salt_file

        return self.hash_salts


    def load_hash_salt(self, salt_file):
        """ Load from disk the hash salts.

        Loads from the file in input the hash salts, one salt per line.
        Hashes are stored encoded in base64, and need to be decoded.

        Args:
            salt_file:  the file where to load the hash salts from.

        Returns:
            The hash salts (in a list).
        """

        self.salt_file = salt_file
        self.hash_salts = []

        self.salt_file_lines = self.salt_file.readlines()

        for self.line in self.salt_file_lines:
            self.hash_salts.append(base64.b64decode(self.line.strip()))

        del self.salt_file_lines
        del self.line
        self.salt_file.close()
        del self.salt_file

        return self.hash_salts


    def insert(self, element, area):
        """ Inserts an element into the SBF filter.

        Maps a single element (a string) into the SBF. For each hash function,
        internal method set_cell is called, passing elements coupled with the area
        labels. The elements MUST be passed following the ascending-order of area
        labels. If this is not the case, the self-collision calculation (done by
        set_cell) will likely be wrong.

        Args:
            element:    element to be mapped (as a string).
            area:       the area label (int).
        """

        self.element = element
        self.area = area

        # Computes the hash digest of the input 'num_hashes' times; each
        # iteration combines the input string with a different hash salt
        for self.i in range(0, self.num_hashes):

            # XOR byte by byte (as char by char) of the element string with the salt
            self.buffer = bytes(''.join([ chr(ord(a) ^ b) for (a,b) in zip(self.element, self.hash_salts[self.i]) ]), 'latin-1')

            # Initializes the hash function according to the hash family
            # it is faster to do this manually than to use
            # self.m = hashlib.new(self.hash_family)
            # although it does not work for MD4
            if (self.hash_family in ['md5', 'MD5']):
                self.m = hashlib.md5()
            elif (self.hash_family in ['sha', 'SHA', 'sha1', 'SHA1']):
                self.m = hashlib.sha1()
            else:
                self.m = hashlib.new(self.hash_family)

            self.m.update(self.buffer)

            # We allow a maximum SBF mapping of 32 bit (resulting in 2^32 cells).
            # Thus, the hash digest is truncated after the first 4 bytes.
            self.digest = self.m.digest()[:4]

            #self.index = int.from_bytes(self.digest, byteorder=byteorder) % (pow(2,self.bit_mapping))
            self.index = int(int.from_bytes(self.digest, byteorder=byteorder) / (pow(2,(self.MAX_BIT_MAPPING - self.bit_mapping))))

            self.set_cell(self.index, self.area);

        self.members += 1
        self.area_members[self.area] += 1

        del self.buffer
        del self.digest
        del self.m
        del self.i
        del self.element


    def insert_from_file(self, dataset_path, dataset_delimiter = ','):
        """ Inserts in the SBF filter a set of elements from a file.

        Insert the elements from a dataset CSV file (dataset_path). The CSV has
        one area-element couple per line (in this order), seprated by the value
        dataset_delimiter, which defaults to ',').

        Args:
            dataset_path:       the path to the CSV file containing the elements
                                to be inserted into the filter.
            dataset_delimiter:  the delimiter character used in the CSV (optional,
                                defaults to ',').
        """

        self.dataset_path = dataset_path
        self.dataset_delimiter = dataset_delimiter

        with open(self.dataset_path, 'r') as self.dataset_file:
            self.dataset_reader = csv.reader(self.dataset_file, delimiter = self.dataset_delimiter)
            for self.row in self.dataset_reader:
                self.insert(self.row[1], int(self.row[0]))

        self.insert_file_list.append(self.dataset_path)

        del self.dataset_delimiter
        del self.dataset_path
        self.dataset_file.close()
        del self.dataset_file
        del self.dataset_reader
        del self.row


    def set_cell(self, index, area):
        """ Sets a cell in the filter to the specified area label.

        Sets the cell to the specified input (the area label). This method is called
        by insert with the cell index, and the area label. It manages the two
        different possible cell sizes (one or two bytes) automatically set during
        filter construction.

        Args:
            index:  the index value poitning to the filter cell to be set.
            area:   the area label to set the cell to.

        Raises:
            AttributeError: The area label is out of bounds.
        """

        self.index = index
        self.area = area

        if (self.area <= 0) or (self.area > self.MAX_AREA_NUMBER):
            raise AttributeError("Invalid area number.")

        # Collisions handling
        self.cell_value = self.filter[self.index]
        if (self.cell_value == 0):
            # Sets cell value
            self.filter[self.index] = self.area
            self.area_cells[self.area] += 1
        elif (self.cell_value < self.area):
            # Sets cell value
            self.filter[self.index] = self.area
            self.collisions += 1
            self.area_cells[self.area] += 1
            self.area_cells[self.cell_value] -= 1
        elif (self.cell_value == self.area):
            self.collisions += 1
            self.area_self_collisions[self.area] += 1
        else:
            # This condition should never be reached as long as elements are
            # processed in ascending order of area label. Self-collisions may
            # contain several miscalculations if elements are not passed
            # following the ascending order of area labels.
            self.collisions += 1

        del self.index
        del self.cell_value


    def check(self, element):
        """ Checks if an element is a member of the areas mapped into the SBF.

        Verifies weather the input element belongs to one of the mapped areas.
        Returns the area label (i.e. the identifier of the set) if the element
        belongs to an area, 0 otherwise.

        Args:
            element:    the element (string) to be checked against the filter.

        Returns:
            The area the element belongs to, or 0 if the element is not a
            member of any area.
        """

        self.element = element
        self.area = 0
        self.current_area = 0

        # Computes the hash digest of the input 'num_hashes' times; each
        # iteration combines the input string with a different hash salt
        for self.i in range(0, self.num_hashes):

            # XOR byte by byte (as char by char) of the element string with the salt
            self.buffer = bytes(''.join([ chr(ord(a) ^ b) for (a,b) in zip(self.element, self.hash_salts[self.i]) ]), 'latin-1')

            # Initializes the hash function according to the hash family
            # it is faster to do this manually than to use
            # self.m = hashlib.new(self.hash_family)
            # although it does not work for MD4
            if (self.hash_family in ['md5', 'MD5']):
                self.m = hashlib.md5()
            elif (self.hash_family in ['sha', 'SHA', 'sha1', 'SHA1']):
                self.m = hashlib.sha1()
            else:
                self.m = hashlib.new(self.hash_family)

            self.m.update(self.buffer)

            # We allow a maximum SBF mapping of 32 bit (resulting in 2^32 cells).
            # Thus, the hash digest is truncated after the first 4 bytes.
            self.digest = self.m.digest()[:4]

            #self.index = int.from_bytes(self.digest, byteorder=byteorder) % (pow(2,self.bit_mapping))
            self.index = int(int.from_bytes(self.digest, byteorder=byteorder) / (pow(2,(self.MAX_BIT_MAPPING - self.bit_mapping))))

            self.current_area = self.filter[self.index]

            # If one hash points to an empty cell, the element does not belong
            # to any set.
            if (self.current_area == 0):
                return 0
            # Otherwise, stores the lower area label, among those which were returned
            elif (self.area == 0):
                self.area = self.current_area
            elif (self.current_area < self.area):
                self.area = self.current_area

        del self.buffer
        del self.digest
        del self.m
        del self.i
        del self.current_area
        del self.index
        del self.element

        return self.area


    def check_from_file(self, dataset_path, dataset_delimiter = ','):
        """ Checks against the SBF filter a set of elements from a file.

        Checks an existing filter using the elements from a dataset CSV file (path).
        The CSV has one area-element couple per line (in this order), seprated by the
        value dataset_delimiter, which defaults to ','). If the CSV file only has one
        value, it assumes the value is the element, and checks it against the filter.

        It stores the check results in two variables:
            check_results:  a list with the check results: if the CSV specifies
                            both areas and elements, the list contains, for each
                            element: a boolean value (returned area matches the
                            CSV area), the expected area, and the returned area.
                            If the CSV specifies only elements, returns the numeric
                            value of the returned area.
            check_results_stats: a list with the results statistics: if the CSV
                            specifies both areas and elements, the list contains
                            the number of matching elements (True) and the number
                            of non-matching elements (False, not in any area). If
                            the CSV specifies only elements, the list contains,
                            for each area, the number of elements recognized as
                            members of the area.

        TODO: returns the variables above.

        Args:
            dataset_path:       the path to the CSV file containing the elements
                                to be checked against the filter.
            dataset_delimiter:  the delimiter character used in the CSV (optional,
                                defaults to ',').

        Raises:
            IOError:    The CSV file is malformed(either more than two values per
            line, or an inconsistent number of values)(either more than two values
            per line, or an inconsistent number of values).
        """

        self.dataset_path = dataset_path
        self.dataset_delimiter = dataset_delimiter
        self.check_list = []
        self.check_results = []

        # Creates a list with the contents of the file for better handling
        with open(self.dataset_path, 'r') as self.dataset_file:
            self.dataset_reader = csv.reader(self.dataset_file, delimiter = self.dataset_delimiter)
            for self.row in self.dataset_reader:
                self.check_list.append(self.row)

        # CSV specifies only elements (one value per line)
        if (len(self.check_list) == sum(len(self.x) for self.x in self.check_list)):
            self.check_results_stats = [[self.x, 0] for self.x in range(self.num_areas+1)]
            for self.check_index, self.check_value in enumerate(self.check_list):
                self.check_val = self.check(self.check_value[0])
                self.check_results.append(self.check_val)
                self.check_results_stats[self.check_val][1] += 1
        # CSV specifies both areas and elements (two values per line)
        elif ((2*len(self.check_list)) == (sum(len(self.x) for self.x in self.check_list))):
            self.check_results_stats = [[True, 0], [False, 0]]
            for self.check_index, self.check_value in enumerate(self.check_list):
                self.check_val = self.check(self.check_value[1])
                self.check_bool = (str(self.check_val) == self.check_value[0])
                self.check_results.append([self.check_bool, self.check_value[0], self.check_val])
                if (self.check_bool):
                    self.check_results_stats[0][1] += 1
                else:
                    self.check_results_stats[1][1] += 1
            del self.check_val
            del self.check_bool
        # CSV specifies either more than two values per line, or an inconsistent number of values
        else:
            raise IOError("CSV file is malformed (either more than two values per line, or an inconsistent number of values).")

        del self.dataset_path
        del self.dataset_delimiter
        self.dataset_file.close()
        del self.dataset_file
        del self.dataset_reader
        del self.row
        del self.check_list
        del self.check_index
        del self.check_value
        del self.x


    def print_filter(self, mode, precision = 5):
        """ Prints the filter and related statistics to the standart output.

        Prints the filter and related statistics (according to the specified
        operation mode) to the standart output.

        Args:
            mode:       If 0, prints the SBF statistics only; if 1, prints the
                        SBF statistics as well as the full SBF content
            precision:  Sets the precision (number of decimal places) to use when
                        printing float values.

        Raises:
            AttributeError: The mode argument is invalid (not in [0,1]).
        """

        self.mode = mode
        self.precision = precision

        if (self.mode not in [0,1]):
                    raise AttributeError("Invalid mode.")

        self.expected_area_cells()
        self.compute_area_fpp()
        self.compute_apriori_area_fpp()
        self.compute_apriori_area_isep()
        self.compute_area_isep()

        print('\nSpatial Bloom Filter stats:')
        print('---------------------------')
        print('Hash details:')
        print(' - Hash family: ' + str(self.hash_family))
        print(' - Number of hash runs: ' + str(self.num_hashes))
        print('Filter details:')
        print(' - Number of cells: ' + str(self.num_cells))
        print(' - Size in bytes: ' + str(self.cell_size * self.num_cells))
        print(' - Filter sparsity: ' + str('{:.{prec}f}'.format(round(self.filter_sparsity(), self.precision), prec=self.precision)))
        print(' - Filter a-priori fpp: ' + str('{:.{prec}f}'.format(round(self.filter_apriori_fpp(), self.precision), prec=self.precision)))
        print(' - Filter fpp: ' + str('{:.{prec}f}'.format(round(self.filter_fpp(), self.precision), prec=self.precision)))
        print(' - Filter a-priori safeness probability: ' + str('{:.{prec}f}'.format(round(self.safeness, self.precision), prec=self.precision)))
        print(' - Number of mapped elements: ' + str(self.members))
        print(' - Number of hash collisions: ' + str(self.collisions))

        if (self.mode == 1):
            print('\nFilter:')
            print('-------')
            for self.i in range(0, self.num_cells):
                # For readability purposes, we print a line break after 16 cells
                if (self.i % 16 == 0) and (self.i != 0):
                    print('')
                print(str(self.filter[self.i]).rjust(3) + " ", end='')
            del self.i

        print('\n')
        print('Area properties:')
        print('----------------')
        for self.j in range(1, self.num_areas + 1):
            self.potential_elements = (self.area_members[self.j] * self.num_hashes) - self.area_self_collisions[self.j]
            print('Area ' + str(self.j).rjust(len(str(self.num_areas))) + ': ' \
                  + str(self.area_members[self.j]) + ' members, ' \
                  + str(round(self.area_expected_cells[self.j])) + ' expected cells, ' \
                  + str(self.area_cells[self.j]) + ' cells out of ' \
                  + str(self.potential_elements) + ' potential (' \
                  + str(self.area_self_collisions[self.j]) + ' self-collisions)')

        print('\nEmersion, FPP and ISEP:\n')
        for self.j in range(1, self.num_areas + 1):
            print('Area ' + str(self.j).rjust(len(str(self.num_areas))) + \
                  ': expected emersion ' + str('{:.{prec}f}'.format(round(self.expected_area_emersion(self.j), self.precision), prec=self.precision)) + \
                  ', emersion ' + str('{:.{prec}f}'.format(round(self.area_emersion(self.j), self.precision), prec=self.precision)) + \
                  ', a-priori fpp ' + str('{:.{prec}f}'.format(round(self.area_apriori_fpp[self.j], self.precision), prec=self.precision)) + \
                  ', fpp ' + str('{:.{prec}f}'.format(round(self.area_fpp[self.j], self.precision), prec=self.precision)) + \
                  ', a-priori isep ' + str('{:.{prec}f}'.format(round(self.area_apriori_isep[self.j], self.precision), prec=self.precision)) + \
                  ', expected ise ' + str('{:.{prec}f}'.format(round((self.area_apriori_isep[self.j] * self.area_members[self.j]), self.precision), prec=self.precision)) + \
                  ', isep ' + str('{:.{prec}f}'.format(round(self.area_isep[self.j], self.precision), prec=self.precision)) + \
                  ', a-priori safep ' + str('{:.{prec}f}'.format(round(self.area_apriori_safep[self.j], self.precision), prec=self.precision)))

        del self.j
        del self.mode
        del self.potential_elements


    def save_filter(self, mode, filter_path = '', precision = 5):
        """ Saves the filter and related statistics onto a CSV file.

        Saves to disk the filter or its statistics (according to the specified
        operation mode) to the specified path.

        Args:
            filter_path:    the path to the file where to store the filter
                            information.
            mode:           If 0, writes the SBF metadata (CSV: key,value);
                            if 1, writes the SBF cells (CSV: value).
            precision:      Sets the precision (number of decimal places) to use
                            when printing float values.

        Raises:
            AttributeError: The mode argument is invalid (not in [0,1]).
            OSError:        The file cannot be created.
        """

        self.mode = mode
        self.filter_path = filter_path
        self.precision = precision

        if (self.filter_path == ''):
            self.filter_path = 'sbf-stats-' + datetime.now().strftime("%Y%m%d-%H%M%S") + '.csv'

        # Tries to load the hash salts from the specified file
        try:
            with open(self.filter_path, 'w') as self.filter_file:

                if (self.mode == 0):

                    self.compute_area_fpp()
                    self.compute_apriori_area_fpp()
                    self.compute_apriori_area_isep()
                    self.compute_area_isep()
                    self.expected_area_cells()

                    self.filter_file.write("hash_family" + ";" + self.hash_family + "\n")
                    self.filter_file.write("hash_number" + ";" + str(self.num_hashes) + "\n")
                    self.filter_file.write("area_number" + ";" + str(self.num_areas) + "\n")
                    self.filter_file.write("bit_mapping" + ";" + str(self.bit_mapping) + "\n")
                    self.filter_file.write("cells_number" + ";" + str(self.num_cells) + "\n")
                    self.filter_file.write("cell_size" + ";" + str(self.cell_size) + "\n")
                    self.filter_file.write("byte_size" + ";" + str(self.cell_size*self.num_cells) + "\n")
                    self.filter_file.write("members" + ";" + str(self.members) + "\n")
                    self.filter_file.write("collisions" + ";" + str(self.collisions) + "\n")
                    self.filter_file.write("sparsity" + ";" + str('{:.{prec}f}'.format(round(self.filter_sparsity(), self.precision), prec=self.precision)) + "\n")
                    self.filter_file.write("a-priori fpp" + ";" + str('{:.{prec}f}'.format(round(self.filter_apriori_fpp(), self.precision), prec=self.precision)) + "\n")
                    self.filter_file.write("fpp" + ";" + str('{:.{prec}f}'.format(round(self.filter_fpp(), self.precision), prec=self.precision)) + "\n")
                    self.filter_file.write("a-priori safeness probability" + ";" + str('{:.{prec}f}'.format(round(self.safeness, self.precision), prec=self.precision)) + "\n")
                    self.filter_file.write("area;members;expected cells;self-collisions;cells;expected emersion;emersion;a-priori fpp;fpp;a-priori isep;expected ise;isep;a-priori safep\n")

                    # area-related parameters:
                    # area, members, expected cells, self-collisions, cells, expected emersion, emersion, apriori_fpp, fpp, apriori_isep, expected ise, isep
                    for self.j in range(1, self.num_areas+1):
                        self.filter_file.write(str(self.j) + ";" + \
                                               str(self.area_members[self.j]) + ";" + \
                                               str(round(self.area_expected_cells[self.j])) + ";" + \
                                               str(self.area_self_collisions[self.j]) + ";" + \
                                               str(self.area_cells[self.j]) + ";" + \
                                               str('{:.{prec}f}'.format(round(self.expected_area_emersion(self.j), self.precision), prec=self.precision)) + ";" + \
                                               str('{:.{prec}f}'.format(round(self.area_emersion(self.j), self.precision), prec=self.precision)) + ";" + \
                                               str('{:.{prec}f}'.format(round(self.area_apriori_fpp[self.j], self.precision), prec=self.precision)) + ";" + \
                                               str('{:.{prec}f}'.format(round(self.area_fpp[self.j], self.precision), prec=self.precision)) + ";" + \
                                               str('{:.{prec}f}'.format(round(self.area_apriori_isep[self.j], self.precision), prec=self.precision)) + ";" + \
                                               str('{:.{prec}f}'.format(round((self.area_apriori_isep[self.j] * self.area_members[self.j]), self.precision), prec=self.precision)) + ";" + \
                                               str('{:.{prec}f}'.format(round(self.area_isep[self.j], self.precision), prec=self.precision)) + ";" + \
                                               str('{:.{prec}f}'.format(round(self.area_apriori_safep[self.j], self.precision), prec=self.precision)) + "\n")

                    del self.j

                elif(self.mode == 1):
                    for self.i in range(0, self.num_cells):
                        self.filter_file.write(str(self.filter[self.i]) + "\n")
                    del self.i

                else:
                    raise AttributeError("Invalid mode.")

                self.filter_file.close()
                del self.filter_file

        except IOError:
            # File couldn't be opened
            raise OSError(2, "File cannot be opened or written into", self.filter_path)

        del self.filter_path
        del self.mode


    def compute_area_fpp(self):
        """ Computes a-posteriori false positives probability for each area.

        Computes the a-posteriori false positives probability (fpp) for each
        area. This method needs to be called manually after the last insert
        operation in order to generate the statistics correctly. Alternatively,
        it is called by both the print_filter and the save_filter (when saving
        the filter statistics).

        Returns:
            The list of a-posteriori false positives probability (fpp) for the areas.
        """

        self.area_fpp = [0]*(self.num_areas + 1)

        for self.i in range(self.num_areas, 0, -1):

            self.c = 0

            for self.j in range(self.i, self.num_areas+1):
                self.c += self.area_cells[self.j]

            self.p = self.c / self.num_cells
            self.area_fpp[self.i] = pow(self.p, self.num_hashes)

            for self.j in range(self.i, self.num_areas):
                self.area_fpp[self.i] -= self.area_fpp[self.j + 1]

            if (self.area_fpp[self.i] < 0):
                self.area_fpp[self.i] = 0

        del self.j
        del self.c
        del self.p
        del self.i

        return self.area_fpp


    def compute_apriori_area_fpp(self):
        """ Computes a-priori false positives probability for each area.

        Computes the a-priori false positives probability (fpp) for each
        area. This method needs to be called manually after the last insert
        operation in order to generate the statistics correctly. Alternatively,
        it is called by both the print_filter and the save_filter (when saving
        the filter statistics).

        Returns:
            The list of a-priori false positives probability (fpp) for the areas.
        """

        self.area_apriori_fpp = [0]*(self.num_areas + 1)

        for self.i in range(self.num_areas, 0, -1):

            self.c = 0
            self.p = 0

            for self.j in range(self.i, self.num_areas+1):
                self.c += self.area_members[self.j]

            self.p = 1 - (1  / self.num_cells)

            self.p = 1 - pow(self.p, (self.num_hashes * self.c))

            self.p = pow(self.p, self.num_hashes)

            self.area_apriori_fpp[self.i] = self.p

            for self.j in range(self.i, self.num_areas):
                self.area_apriori_fpp[self.i] -= self.area_apriori_fpp[self.j + 1]

            if (self.area_apriori_fpp[self.i] < 0):
                self.area_apriori_fpp[self.i] = 0

        del self.j
        del self.c
        del self.p
        del self.i

        return self.area_apriori_fpp


    def expected_area_cells(self):
        """ Computes the expected number of cells for each area

        Computes the expected number of cells for each area. The expected value
        depends only on the filter properties (size, number of hash functions,
        number of sets and elements, etc.) and not an actual instance of the
        filter.

        Returns:
            The list of expected number of cells for the areas.
        """

        self.area_expected_cells = [0]*(self.num_areas + 1)

        for self.i in range(self.num_areas, 0, -1):

            self.nfill = 0

            for self.j in range(self.i+1, self.num_areas+1):
                self.nfill += self.area_members[self.j]

            self.p1 = 1 - (1  / self.num_cells)

            self.p2 = pow(self.p1, (self.num_hashes * self.nfill))

            self.p1 = 1 - pow(self.p1, (self.num_hashes * self.area_members[self.i]))

            self.p1 = self.num_cells * self.p1 * self.p2

            self.area_expected_cells[self.i] = self.p1

        del self.nfill
        del self.p1
        del self.p2
        del self.i
        del self.j

        return self.area_expected_cells


    def compute_apriori_area_isep(self):
        """ Computes a-priori inter-set error probability for each area.

        Computes the a-priori inter-set error probability (isep) for each
        area. This method needs to be called manually after the last insert
        operation in order to generate the statistics correctly. Alternatively,
        it is called by both the print_filter and the save_filter (when saving
        the filter statistics).

        Returns:
            The list of a-priori inter-set error probability (isep) for the areas.
        """

        self.area_apriori_isep = [0]*(self.num_areas + 1)
        self.area_apriori_safep = [0]*(self.num_areas + 1)
        self.safeness = 1

        for self.i in range(self.num_areas, 0, -1):

            self.nfill = 0
            self.p1 = 0
            self.p2 = 0

            for self.j in range(self.i+1, self.num_areas+1):
                self.nfill += self.area_members[self.j]

            self.p1 = 1 - (1  / self.num_cells)

            self.p1 = 1 - pow(self.p1, (self.num_hashes * self.nfill))

            self.p1 = pow(self.p1, self.num_hashes)

            self.p2 = 1 - self.p1
            self.p2 = pow(self.p2, self.area_members[self.i])

            self.safeness *= self.p2

            self.area_apriori_isep[self.i] = self.p1
            self.area_apriori_safep[self.i] = self.p2

        del self.nfill
        del self.p1
        del self.p2
        del self.i
        del self.j

        return self.area_apriori_isep


    def compute_area_isep(self):
        """ Computes a-posteriori inter-set error probability for each area.

        Computes the a-posteriori inter-set error probability (isep) for each
        area. This method needs to be called manually after the last insert
        operation in order to generate the statistics correctly. Alternatively,
        it is called by both the print_filter and the save_filter (when saving
        the filter statistics).

        Returns:
            The list of a-posteriori inter-set error probability (isep) for the areas.
        """

        self.area_isep = [0]*(self.num_areas + 1)

        for self.i in range(self.num_areas, 0, -1):

            self.p = 1 - self.area_emersion(self.i)
            self.p = pow(self.p, self.num_hashes)

            self.area_isep[self.i] = self.p

        del self.p

        return self.area_isep


    def filter_sparsity(self):
        """ Returns the sparsity of the SBF.

        Returns the sparsity (ratio of non zero cells to the total number of cells)
        for the entire SBF filter.

        Returns:
            The filter sparsity.
        """

        self.sparsity_sum = 0

        for self.i in range(1, self.num_areas + 1):
            self.sparsity_sum += self.area_cells[self.i]

        return 1 - (self.sparsity_sum / self.num_cells)


    def filter_fpp(self):
        """ Computes a-posteriori false positives probability for the filter.

        Computes the a-posteriori false positive probability over the entire
        filter (i.e. not area-specific).

        Returns:
            The filter a-posteriori false positives probability (fpp).
        """

        self.c = 0

        # Counts non-zero cells
        for self.i in range(1, self.num_areas + 1):
            self.c += self.area_cells[self.i]

        self.p = self.c / self.num_cells

        return pow(self.p, self.num_hashes)


    def filter_apriori_fpp(self):
        """ Computes a-priori false positives probability for the filter.

        Computes the a-priori false positive probability over the entire
        filter (i.e. not area-specific).

        Returns:
            The filter a-priori false positives probability (fpp).
        """

        self.p = 1 - (1 / self.num_cells)

        self.p = 1 - pow(self.p, (self.num_hashes * self.members))

        self.p = pow(self.p, self.num_hashes)

        return self.p


    def area_emersion(self, area):
        """ Computes the emersion value for an area.

        Computes the emersion value for the input area. The emersion indicates how
        much the values in the cells that maps elements of an area store the area
        label value (as opposed to a higher value). It is the ratio between the
        cells actually set to the area label in the filter and the cells that would
        be theoreticall set to the area label if no collision from higher areas
        happen.

        Args:
            area:   the area for which to calculate the emersion value.

        Returns:
            The emersion value (float).
        """

        self.area = area

        if (self.area_members[self.area] == 0):
            return -1
        else:
            return (self.area_cells[self.area] / ((self.area_members[self.area] * self.num_hashes) - self.area_self_collisions[self.area]))


    def expected_area_emersion(self, area):
        """ Computes the expected emersion value for an area.

        Computes the expected emersion value for the input area. The expected
        value depends only on the filter properties (size, number of hash functions,
        number of sets and elements, etc.) and not an actual instance of the filter.
        The emersion indicates how much the values in the cells that maps elements
        of an area store the area label value (as opposed to a higher value). It
        is the ratio between the cells actually set to the area label in the filter
        and the cells that would be theoreticall set to the area label if no
        collision from higher areas happen.

        Args:
            area:   the area for which to calculate the emersion value.

        Returns:
            The expected emersion value (float).
        """

        self.area = area
        self.nfill = 0

        for self.i in range(self.area+1, self.num_areas+1):
                self.nfill += self.area_members[self.i]

        self.p = 1 - (1  / self.num_cells)

        self.p = pow(self.p, (self.num_hashes * self.nfill))

        del self.nfill

        return self.p
