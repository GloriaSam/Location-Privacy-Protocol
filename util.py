#!/usr/bin/env python

import json
import socket

BUFF_SIZE = 4096

def packet_send(sock, msg_type, data=""):
    if data is not '':
        json_data = json.dumps(data)
        size = len(json_data)
        len_payload = str(len(json_data))
        format_len = len_payload.rjust(16)
        tosend = msg_type[0] + format_len + json_data
    else:
        tosend = msg_type[0] + '0'.rjust(16)
    sock.sendall(tosend.encode('ascii'))

def packet_recv(conn):
    header = conn.recv(17).decode('ascii')
    msg_type = header[0]
    msg_len = int(header[1:])
    msg = ''
    if msg_len != 0:
        json_msg = b''
        while msg_len > BUFF_SIZE:
            json_msg += conn.recv(BUFF_SIZE)
            msg_len -= BUFF_SIZE
        json_msg += conn.recv(msg_len)
        json_msg_str = json_msg.decode('ascii')
        msg = json.loads(json_msg_str)
    return msg, header[0]
