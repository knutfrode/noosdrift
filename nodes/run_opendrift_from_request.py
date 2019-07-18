#!/usr/bin/env python

import sys
from datetime import datetime
from collections import OrderedDict
import json
from opendrift.readers import reader_netCDF_CF_generic
from noosdrift.seed.noos_seed import noos_seed, noos_seed_from_json


def run_opendrift_simulation_request(request_json_file):

    status_codes = {
        0: 'Model simulation successfully completed (no error)',
        1: 'ERROR: initial position out of model domain',
        2: 'ERROR: initial position on land',
        3: 'ERROR: Simulation start and/or end time are not in the forcing availability period [today-4 days, today+4 days]',
        4: 'ERROR : Release time of Lagrangian particle out of the simulation start time and end time windows',
        5: 'ERROR : Drifter type unknown or not available in the model',
        6: 'ERROR: Model cannot handle the requested set-up -> set-up has been adapted',
        7: 'ERROR: any other error in the model pre-processing',
        8: 'ERROR: any error in the model processing',
        9: 'ERROR: any error in the model post-processing : preparation of the model output'
    }

    current_sources = {
        'cmems_nws7': 'sample_forcing_data/nws_current_14Jul2019.nc',
        'norkyst': 'sample_forcing_data/norkyst_current_14Jul2019.nc'}
    wind_sources = {
        'ecmwf': 'sample_forcing_data/ec_wind_14Jul2019.nc',
        'arome': 'sample_forcing_data/arome_wind_14Jul2019.nc'}

    oil_type_mapping = {  # from Noos-name to OpenDrift oil name
        'gasoline': '*GENERIC GASOLINE',
        'kerosene': 'FUEL OIL NO.1 (JET FUEL A)',
        'light crude oil': '*GENERIC LIGHT CRUDE',
        'diesel oil': '*GENERIC DIESEL',
        'heavy crude oil': '*GENERIC HEAVY CRUDE',
        'fuel oil no. 6': '*GENERIC FUEL OIL No. 6',
        'heavy fuel oil': '*GENERIC HEAVY FUEL OIL'}

    def return_status_code(status_code, comment=None,
                           resultfile=None):
        print('Returning with status code %i: %s' %
                (status_code, status_codes[status_code]))
        j['simulation_result'] = {'status_code': status_code}
        if resultfile is not None:
            j['simulation_result']['result_file_name'] = resultfile
        if comment is not None:
            j['simulation_result']['comment'] = comment
            print(comment)
        with open('out.json', 'w') as outfile:
            json.dump(j, outfile, indent=2)
        sys.exit(status_code)

    # Import JSON request file
    f = open(request_json_file)
    if f is None:
        raise ValueError('File does not exist')

    file = open(request_json_file)
    j = json.load(file, object_pairs_hook=OrderedDict)

    # Read some parameters
    end_time = datetime.strptime(
        j['simulation_description']['simulation_end_time'],
        '%Y-%m-%dT%H:%M:%Sz')

    # Prepare readers for forcing data
    current_source = j['model_set_up']['ocean_forcing']
    if current_source not in current_sources:
        return_status_code(9,
            comment='Current source %s not available' %
                current_source)
    else:
        current_URL = current_sources[current_source]
    wind_source = j['model_set_up']['wind_forcing']
    if wind_source not in wind_sources:
        return_status_code(9,
            comment='Wind source %s not available' %
                wind_source)
    else:
        wind_URL = wind_sources[wind_source]

    current_reader = reader_netCDF_CF_generic.Reader(current_URL)
    wind_reader = reader_netCDF_CF_generic.Reader(wind_URL)

    # Check that readers cover the requested time and area
    lons,lats,times = noos_seed_from_json(j['initial_condition'])
    if len(current_reader.covers_positions(
            lons.mean(), lats.mean())[0])==0:
        return_status_code(1,
            'Current reader does not cover seeding location')
    if len(wind_reader.covers_positions(
            lons.mean(), lats.mean())[0])==0:
        return_status_code(1,
            'Wind reader does not cover seeding location')

    # Check if oil or leeway simulation
    seed_kwargs = {}
    drifter_name = j['drifter']['drifter_name']
    if j['drifter']['drifter_type'] == 'oil':
        from opendrift.models.openoil3D import OpenOil3D
        seed_kwargs['oiltype'] = oil_type_mapping[drifter_name]
        o = OpenOil3D(weathering_model='noaa')
    elif j['drifter']['drifter_type'] == 'object':
        seed_kwargs['objectType'] = drifter_name
        from opendrift.models.leeway import Leeway
        o = Leeway()

    # Adding readers
    o.add_reader([current_reader, wind_reader])

    # Seed elements at requested time and positions
    o.seed_elements(lons, lats, time=times, **seed_kwargs)

    # Run simulation, and save output to netCDF file
    resultfile = j['simulation_description']['output_path'] + \
                    '/noosdrift_%s_opendrift_%s_%s.nc' % (
                    j['simulation_description']['request_id'],
                    current_source, wind_source)

    try:
        o.run(outfile=resultfile, end_time=end_time)
        # simulation is successful
        return_status_code(0, resultfile=resultfile)
    except SystemExit:
        sys.exit(0)
    except:
        return_status_code(8)  # simulation failed

    # Save status to JSON file
    output_JSON = {'simulation_result' :
        {'status_code' : status_code,
         'result_file_name' : outname + '.json'}}
    with open(outname + '.json', 'w') as outfile:
        json.dump(output_JSON, outfile)


if __name__ == '__main__':
    
    try:
        json_request_file = sys.argv[1]
    except Exception as e:
        print(e)
        raise ValueError('Usage: %s <json_request_file>' %
                         sys.argv[0])
    run_opendrift_simulation_request(json_request_file)
