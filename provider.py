#!/usr/bin/env python

import math
import util
import socket
import numpy as np

from sbf import sbf
from phe import paillier
from joblib import Parallel, delayed
from random import randint
from random import shuffle
from collections import Counter

class SbfProvider(object):
   
    def __init__(self):
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.bind(('localhost', 2323))
        self.sock.listen(5)

    def close_stream(self):
        """This method will close the socket stream."""
        self.sock.close()
        
    def create_sbf(self):
        """Inits the provider's sbf.

        The function asks to the user:
            - name of the areas data file
            - desired false positive percentage (fpp)
            - desired hash family
        It determines the number of areas and the number of elements, then
        it calculates the best size, and number of hash run for the sbf
        creation.

        Returns:
            sbf_vector (sbf): provider's sbf
            self.bit_mapping (int): sbf vector size
            hash_family (str): hash function used for sbf
            self.hash_runs (int): number of the hash runs
        """
        file_path = input("Enter the path for the verification test\n").rstrip()

        self.num_elements = 0
        areas = []
        with open(file_path) as f:
            for line in f:
                self.num_elements += 1
                areas.append(line.split(",")[0])
            self.num_areas = len(Counter(areas))

        max_fpp = float(input(
            "Enter the desired false positive probability (fpp) [less then 0.001]\n").rstrip())
        cells = math.ceil(-self.num_elements * math.log(max_fpp)) / (math.log(2) ** 2)
        self.bit_mapping = math.ceil(math.log(cells, 2))
        self.hash_runs = math.ceil((cells / self.num_elements) * math.log(2))
        self.hash_family = input(
            "Enter the desired hash family [sha1, md4, md5]\n").rstrip()

        print("\nNumber of elements: ", self.num_elements)
        print("Number of areas: ", self.num_areas)
        print("Resulting array size: 2^", self.bit_mapping)
        print("Resulting number of hash runs: ", self.hash_runs)

        # TODO: How should I set the salt?
        self.sbf_vector = sbf(self.bit_mapping, self.hash_family, self.hash_runs, self.num_areas, "salt")
        self.sbf_vector.insert_from_file(file_path)
        self.sbf_vector.print_filter(0)

        return self.sbf_vector

    def obfuscate_zeros(self):
        """ Hide the server's AoI to the client, and it hide the client's position to the server

        1. count the no zero values inside the sbf (NZ)
        2. count the number of elements per the N areas Ei = (E1, E2 ... EN)
        3. groups the no zero values, creating M sets of size randomly between min(Ei) and max(Ei) 
        4. put everything in a list 
        5. shuffle the list
        """
        unique, counts = np.unique(self.sbf_vector.filter, return_counts=True)
        self.e_i = dict(zip(unique, counts))
        M = self.e_i[0]
        del self.e_i[0]
        max_e_i = max(self.e_i) 
        min_e_i = min(self.e_i) 
        tot = M
        self.obfuscated_zeros = []
        while (tot > 0):
            size_zero_set = randint(min_e_i, max_e_i)
            if size_zero_set < tot: 
                size_zero_set = tot
            r_value = self.public_key.get_random_lt_n()
            for y in range(size_zero_set):
               self.obfuscated_zeros.append(self.public_key.raw_encrypt(0, r_value)) 
            tot -= size_zero_set 

        shuffle(self.obfuscated_zeros)
        return

    def encrypt_sbf(self):
        """Encrypt sbf with pailler publick key.

        Args:
            sbf_vector
            self.public_key (PaillerPublicKey): pailler public key class
        Returns:
            enc_sbf_vector (list): encrypted sbf

        """
        self.public_key, self.private_key = paillier.generate_paillier_keypair(None, 1024)
        self.obfuscate_zeros()
        
        e_i_r_value = self.e_i 
        for k,v in self.e_i.items():
            e_i_r_value[k] = self.public_key.get_random_lt_n()

        obf_index = 0
        self.enc_sbf_vector = []
        for element in self.sbf_vector.filter:
            if element == 0:
                self.enc_sbf_vector.append(self.obfuscated_zeros[obf_index])
                obf_index += 1
            else:
                self.enc_sbf_vector.append(self.public_key.raw_encrypt(int(element), e_i_r_value[element]))
        
        print("\nFILTER ENCRYPTED") 
        return self.enc_sbf_vector

    def check_user_position(self, num_nozero, mask):
        """Chek if users is inside any area, and in case return it.

        Args:
            enc_client_sbf (list): client's encrypted and shuffled sbf
            num_nozero (int): number of no zero req_types in the client's computed sbf
        Returns:
            client_area (int): client's located area
        """
        if (num_nozero <= 0):
            print("User it outside any area\n")
            return
        dec_client_sbf = Parallel(n_jobs=4)(
            delayed(self.private_key.raw_decrypt)(int(x)) for x in mask if x != 0)
        client_area = min(dec_client_sbf)
        return client_area

    def loop(self):
        """Create socket send sbf datas then waith until recives client datas.
        
        receive S as send sbf data
        receive C as check client position 

        Args:
            enc_sbf_vector (list)
            public_key (PaillierPublicKey)
            bit_mapping (int)
            hash_family (str)
            hash_runs (int)

            mask (list): encrypted and shuffles client's sbf
            non_zero (int): number of no zero values in the client's computed sbf
            """
        client_socket, addr = self.sock.accept()
        while True:
            msg, req_type  = util.packet_recv(client_socket)
            if req_type is 'S':
                print("Received sbf request")
                self.create_sbf() 
                self.encrypt_sbf()
                data = {
                    "enc_sbf_vector": self.enc_sbf_vector,
                    "bit_mapping": self.bit_mapping,
                    "hash_family": self.hash_family,
                    "hash_runs": self.hash_runs
                }
                util.packet_send(client_socket, 'S', data)
            elif req_type is 'C':
                print("Received check request")
                client_area = self.check_user_position(msg['non_zero'], msg['mask'])
                if (client_area <= 0):
                    print(addr, "it outside any area")
                else:
                    print(addr, "is inside area ", client_area)
        client_socket.close()

def main():
    try:
        provider = SbfProvider()
        provider.loop()
    except socket.error:
        print('socket error')
    finally:
        provider.close_stream()

if __name__ == '__main__':
    main()
