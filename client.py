#!/usr/bin/env python

import sys
import socket
import numpy as np
import util
from sbf import sbf
from phe import PaillierPublicKey

class SbfClient(object):
    """

    """
    def __init__(self, HOST="localhost", PORT=2323):
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.connect((HOST, PORT))

    def close_stream(self):
        """This method will close the socket stream."""
        self.sock.close()
        
    def generate_sbf(self, position, sbf_data):
        """
        Generate client sbf, with his position string
        """
        user_sbf = sbf(sbf_data['bit_mapping'], sbf_data['hash_family'], sbf_data['hash_runs'], 1,"salt")
        user_sbf.insert(position, 1)
        non_zero = np.count_nonzero(user_sbf.filter == 1)
        return user_sbf, non_zero

    def mask_message(self, user_sbf_vector, enc_provider_sbf_vector):
        mask = []
        for i,enc in enumerate(enc_provider_sbf_vector):
            mask.append(user_sbf_vector[i] * enc)
        return mask
        
    def loop(self):
        """This is the main method of the client, it will enter
        in a receive/send loop."""
        
        util.packet_send(self.sock, 'S')
        sbf_data, _ = util.packet_recv(self.sock)
        print('Received sbf data')
        
        while True:
            position = input("Enter your position (longitude, latitude) as [xx.xxxx#yy.yyyy]\n").rstrip()
            #TODO: input check
            user_sbf, non_zero = self.generate_sbf(position, sbf_data)
            mask = self.mask_message(user_sbf.filter, sbf_data['enc_sbf_vector'])
            data = {
                    "non_zero": non_zero,
                    "mask": mask
            }
            util.packet_send(self.sock, 'C', data)
        self.close_stream()

def main():
    try:
        cli = SbfClient()
        cli.loop()
    except socket.error:
        print('socket error')
    finally:
        cli.close_stream()

if __name__ == '__main__':
    main()
