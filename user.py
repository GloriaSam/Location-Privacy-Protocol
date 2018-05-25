#!/usr/bin/env python

import socket
import numpy as np
import util
from sbf import sbf

class SbfUser(object):
    """User Class of the location privacy protcol.
       The protocol is realized in the loop() function of this class.
       
       1. The user receives the following sbf info:
            - the encrypted and precomputed SBF Enc(b#)  
            - the set of k hash functions H
            - the value m 
            - the conventional grid E

       2. The user determines his geographic position (long. lat.) at regular
           time intervals or whenever a specific application requires it, then he has to
           selects the corresponding grid region e_u in E. (In this implementation
           we use, as convetional grid, longitude and latitute trimmed at the
           fourth decimal) 
           He creates the sbf b_u#, using the received sbf info, over e_u and
           counts the numbers of non zero values therein. 
           
       3. The user computes e# = Enc(b#)*b_u# using the homomorphic properties
           of the Pailler cryptosystem. Then applies a random permutation over the
           e#. After that he send:
            - non zero values
            - e#
    """
    def __init__(self, HOST="localhost", PORT=2324):
        """This method creates the socket stream."""
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.connect((HOST, PORT))

    def close_stream(self):
        """This method closes the socket stream."""
        self.sock.close()
        
    def generate_sbf(self, position, sbf_data):
        """Generates client sbf over the user's position.
             
        Args:
            position (str): longitude latitude data of the user's postion
            sbf_data (dict): data for the sbf creation

        Returns: 
            tuple: the user's sbf (list) and the number of non zero values
            therein (int)
        """
        user_sbf = sbf(sbf_data['bit_mapping'], sbf_data['hash_family'], sbf_data['hash_runs'], 1, sbf_data['salt'])
        user_sbf.insert(position, 1)
        non_zero = np.count_nonzero(user_sbf.filter == 1)
        return user_sbf, non_zero

    def mask_message(self, user_sbf_vector, provider_sbf_vector_enc):
        """Computes element-wise multiplication:
            mask = provider_sbf_vector_enc * user_sbf_vector.
            
        Args:
            user_sbf_vector (list) the user's sbf
            provider_sbf_vector_enc (list) the encrypted provider's sbf

        Returns:
            mask (list) 
            """
        mask = []
        for i, enc in enumerate(provider_sbf_vector_enc):
            mask.append(user_sbf_vector[i] * enc)
        return mask
        
    def loop(self):
        """Main user's method. It realize the comunication"""
        util.packet_send(self.sock, 'S')
        sbf_data, _ = util.packet_recv(self.sock)
        print('Received sbf data')
        
        while True:
            position = input("Enter your position (longitude, latitude) as [xx.xxxx#yy.yyyy]\n").rstrip()
            #TODO:input check
            user_sbf, non_zero = self.generate_sbf(position, sbf_data)
            mask = self.mask_message(user_sbf.filter, sbf_data['sbf_vector_enc'])
            data = {
                    "non_zero": non_zero,
                    "mask": mask
            }
            util.packet_send(self.sock, 'C', data)
        self.close_stream()

def main():
    try:
        u = SbfUser()
        u.loop()
    except socket.error:
        print('socket error')
    finally:
        u.close_stream()

if __name__ == '__main__':
    main()
