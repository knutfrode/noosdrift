#!/usr/bin/env python

import sys
import json
from opendrift.readers import reader_netCDF_CF_generic
from noosdrift.seed.noos_seed import noos_seed, noos_seed_from_json


def run_opendrift_simulation_request(request_json_file):

    # Import JSON request file
    f = open(request_json_file)
    if f is None:
        raise ValueError('File does not exist')

    file = open(request_json_file)
    j = json.load(file)

    # Check if oil or leeway simulation
    if j['drifter']['drifter_type'] == 'oil':
        from opendrift.models.openoil3D import OpenOil3D as Model
    elif j['drifter']['drifter_type'] == 'object':
        from opendrift.models.leeway import Leeway as Model
    o = Model()  # initiate simulation object

    # Seed elements at requested time and positions
    lons,lats,times = noos_seed_from_json(j['initial_condition'])
    o.seed_elements(lons, lats, time=times)

    # Prepare readers for forcing data
    current_source = j['model_set_up']['ocean_forcing']
    wind_source = j['model_set_up']['wind_forcing']

    if current_source == 'cmems_nws7':
        current_URL = '/path/to/nws_file.nc'
    elif current_source == 'norkyst':
        current_URL = '/path/to/norkyst_file.nc'
    else:
        raise ValueError('Current source %s not available' %
                         current_source)

    current_reader = reader_netCDF_CF_generic.Reader(current_URL)

    # Check that readers cover the requested time and area
    if not current_reader.covers_position(centerlon, centerlat):
        print('Current reader does not cover seeding location')
        status_code = 1
    if not wind_reader.covers_position(centerlon, centerlat):
        print('Wind reader does not cover seeding location')
        status_code = 1

    # Run simulation, and save output to netCDF file
    outname = j['simulation_description']['output_path'] + \
                'noosdrift_%s_opendrift_%s_%s' % (
                j['simulation_description']['request_id'],
                current_source, wind_source)

    try:
        o.run(outfile=outname + '.nc')
        status_code = 0  # simulation is successful
    except:
        status_code = 8  # something went wrong

    # Save status to JSON file
    output_JSON = {'simulation_result' :
        {'status_code' : status_code,
         'result_file_name' : outname + '.json'}}
    with open(outname + '.json', 'w') as outfile:
        json.dump(output_JSON, outfile)


if __name__ == '__main__':
    
    try:
        json_request_file = sys.argv[1]
        run_opendrift_simulation_request(json_request_file)
    except:
        raise ValueError('Usage: %s <json_request_file>' %
                         sys.argv[0])
