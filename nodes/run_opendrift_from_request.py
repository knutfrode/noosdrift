#!/usr/bin/env python

import sys
import opendrift
from noosdrift.seed import noos_seed


def run_opendrift_simulation_request(request_json_file):

    f = open(request_json_file)
    if f is None:
        raise ValueError('File does not exist')

    print f

    print('exiting')
    sys.exit(2)
    print('exited')


if __name__ == '__main__':
    
    try:
        json_request_file = sys.argv[1]
    except:
        raise ValueError('Usage: %s <json_request_file>' %
                         sys.argv[0])

    run_opendrift_simulation_request(json_request_file)
