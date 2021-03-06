#!/usr/bin/env python

import math
import util
import numpy as np
import socket

from sbf import sbf
from phe import paillier
from joblib import Parallel, delayed
from random import randint
from random import shuffle
from collections import Counter

class SbfProvider():
    """Provider Class of the location privacy protocol.
       The protocol is realized in the loop() funcion of this class.
       
       1. The provider generates a pair of asimetric key using Pailler's homomorphic cryptosystem
       2. The provider sends to the user the enrcyption of the precomputed sbf: Enc(b#)
       3. The provider receives the user's shuffled e# = Enc(b#) * b_u# and the
           number of non-zero values in b_u#: z. If z is > of the number of
           non-zero in Dec(e#) the user's position is outside any monitored area.
           Otherwise the smallest non-zero value in Dec(e#) is the area's index
           identifying user's position.  
    """
    def __init__(self, HOST="localhost", PORT=2324):
        """Init socket stream."""

        if (PORT <= 0) or (PORT >65535):
            raise ValueError("Invalid port")
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.bind((HOST, PORT))
        self.sock.listen(5)

    def close_stream(self):
        """Close socket stream."""
        self.sock.close()
        
    def create_sbf(self):
        """Init the provider's sbf.

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
        file_path = input("Enter the path for the verification test [vatican-city-0001.csv]\n").rstrip()

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

        self.hash_salt_path = 'salt'
        self.sbf_vector = sbf(self.bit_mapping, self.hash_family, self.hash_runs, self.num_areas, self.hash_salt_path)
        with open(self.hash_salt_path,'w') as salt_file:
            self.sbf_vector.create_hash_salt(salt_file)
        self.sbf_vector.insert_from_file(file_path)
        self.sbf_vector.print_filter(0)

        return self.sbf_vector

    """
    def inter_salt(self):
        self.sbf_vector.compute_apriori_area_isep()
        while self.sbf_vector.safep
    """
    

    def obfuscate_zeros(self):
        """ Hide the server's AoI to the user, and it hide the user's position to the server

        1. count the no zero values inside the sbf (NZ)
        2. count the number of elements per the N areas Ei = (E1, E2 ... EN)
        3. groups the no zero values, creating M sets of size randomly between min(Ei) and max(Ei) 
        4. put everything in a list 
        5. shuffle the list
        """
        unique, counts = np.unique(self.sbf_vector.filter, return_counts=True)
        # e_i -> key = area_index, value = number_of_elements
        self.e_i = dict(zip(unique, counts))
        M = self.e_i[0]
        del self.e_i[0]
        max_e_i = max(self.e_i.values()) 
        min_e_i = min(self.e_i.values()) 
        tot = M
        self.obfuscated_zeros = []
        while (tot > 0):
            size_zero_set = randint(min_e_i, max_e_i)
            if size_zero_set > tot: 
                size_zero_set = tot
            r_value = self.public_key.get_random_lt_n()
            self.obfuscated_zeros += size_zero_set * [self.public_key.raw_encrypt(0, r_value)]
            tot -= size_zero_set 

        shuffle(self.obfuscated_zeros)
        return

    def encrypt_sbf(self):
        """Encrypt sbf with pailler publick key.

        Returns:
            sbf_vector_enc (list): encrypted sbf
        """
        self.public_key, self.private_key = paillier.generate_paillier_keypair(None, 1024)
        self.obfuscate_zeros()
        
        e_i_enc = self.e_i 
        for k,v in self.e_i.items():
            e_i_r_value = self.public_key.get_random_lt_n()
            e_i_enc[k] = self.public_key.raw_encrypt(int(k), e_i_r_value)

        obf_index = 0
        self.sbf_vector_enc = []
        for e in self.sbf_vector.filter:
            if e == 0:
                self.sbf_vector_enc.append(self.obfuscated_zeros[obf_index])
                obf_index += 1
            else:
                self.sbf_vector_enc.append(e_i_enc[e])
        
        print("\nFILTER ENCRYPTED\n") 
        return self.sbf_vector_enc

    def check_user_position(self, user_non_zero, mask):
        """Chek if users is inside any area, and in case return it.

        Args:
            enc_user_sbf (list): user's encrypted and shuffled sbf
            num_nozero (int): number of no zero values in the user's computed sbf
        Returns:
            user_area (int): user's located area
        
        """
        user_sbf_dec = Parallel(n_jobs=4)(
            delayed(self.private_key.raw_decrypt)(int(x)) for x in mask if x != 0)
        non_zero_dec = np.count_nonzero(user_sbf_dec) 
        if (non_zero_dec < user_non_zero):
            #user outside of any area
            user_area = -1
        else:
            user_area = min(user_sbf_dec)
        return user_area

    def loop(self):
        """Create socket send sbf datas then waith until recives user datas.
        
        receive S as send sbf data
        receive C as check user position 

        Args:
            sbf_vector_enc (list)
            public_key (PaillierPublicKey)
            bit_mapping (int)
            hash_family (str)
            hash_runs (int)
            salt (str)

            mask (list): encrypted and shuffles user's sbf
            non_zero (int): number of no zero values in the user's computed sbf
            """
        user_socket, addr = self.sock.accept()
        while True:
            msg, req_type  = util.packet_recv(user_socket)
            if req_type is 'S':
                print("Received sbf request")
                self.create_sbf()
                self.encrypt_sbf()
                data = {
                    "sbf_vector_enc": self.sbf_vector_enc,
                    "bit_mapping": self.bit_mapping,
                    "hash_family": self.hash_family,
                    "hash_runs": self.hash_runs,
                    "salt": self.hash_salt_path
                }
                util.packet_send(user_socket, 'S', data)
            elif req_type is 'C':
                print("Received check request")
                user_area = self.check_user_position(msg['non_zero'], msg['mask'])
                if (user_area <= 0):
                    print(addr, "user it outside any area")
                else:
                    print(addr, "user is inside area:", user_area)
        user_socket.close()

def main():
    p = SbfProvider()
    try:
        p.loop()
    except socket.error:
        print('socket error')
    finally:
        p.close_stream()

if __name__ == '__main__':
    main()
