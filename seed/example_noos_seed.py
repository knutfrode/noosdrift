#!/usr/bin/env python

from noos_seed import noos_seed_from_json

# Examples from specification_request_package_v2.0_20190208

# From file
#############
lon, lat, time = noos_seed_from_json('seed_example4.2.json', plot=True)
lon, lat, time = noos_seed_from_json('seed_example4.5.json', plot=True)

# From dict
#############
json_example4_1 = {
    'geometry': 'point',
    'lon': 2.159,
    'lat': 52.1,
    'radius': 1000,
    'number': 2500,
    'time': '2018-10-17T13:07:12z'
}
lon, lat, time = noos_seed_from_json(json_example4_1, plot=True)

json_example_cone = {
    'geometry': 'point',
    'lon': [2.159, 2.3],
    'lat': [52.1, 52.12],
    'radius': [10, 200],
    'number': 2500,
    'time': ['2018-10-17T12:00:00z', '2018-10-17T18:00:00z']
}
lon, lat, time = noos_seed_from_json(json_example_cone, plot=True)
