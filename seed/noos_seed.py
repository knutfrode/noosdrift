#!/usr/bin/env python

# Written by Knut-Frode Dagestad (knutfd@met.no) Feb 2019
# for the Noos-Drift project

from datetime import datetime, timedelta
import numpy as np
import pyproj
import json


def noos_seed(lon, lat, time=None, z=None, number=1000, radius=0,
              plot=False, randgen=0):
    '''Provide arrays of lon/lat from shorthand seed notation.

    lon: scalar or array/list with 2 elements
    lat: scalar or array/list with 2 elements
    time: datetime object, or list with 2 dt objects [start, end]
    number: integer, the number of particles
    radius: scalar or array/list with 2 elements. The radius in m
            around the (lon,lat) position, within which the elements
            are seeded. Not absolute bound, but gaussian std.
    randgen: A fixed value will always give the same random numbers.
             Use None to get different numbers each time.
    '''

    if randgen is not None:
        np.random.seed(randgen)
        
    #Seeding - horizontal coordinates
    lon = np.atleast_1d(lon)
    lat = np.atleast_1d(lat)
    radius = np.atleast_1d(radius)
    if len(lon) == 1:
        lon = np.repeat(lon, 2)
    if len(lat) == 1:
        lat = np.repeat(lat, 2)
    if len(radius) == 1:
        radius = np.repeat(radius, 2)

    geod = pyproj.Geod(ellps='WGS84')
    # Centerpoints of "cone":
    conelonlats = geod.npts(lon[0], lat[0], lon[1], lat[1],
                            number, radians=False)
    lon, lat = list(zip(*conelonlats))
    radius = np.linspace(radius[0], radius[1], number)

    # Perturbations
    x = np.random.randn(number)*radius.T
    y = np.random.randn(number)*radius.T
    azimuth = np.degrees(np.arctan2(x, y))
    dist = np.sqrt(x*x+y*y)
    ones = np.ones(number)
    
    lons, lats, az = geod.fwd(lon*ones, lat*ones,
                              azimuth, dist, radians=False)
    
    #Seeding vertical coordinates
    depths = None
    if z is None:
       z = 0.0
    z = np.atleast_1d(z)
    if len(z) == 1:
       z = np.repeat(z, 2)
    depths = np.linspace(z[0], z[1], number)

   
    #Seeding time coordinates
    times = None   
    if time is not None:
        if isinstance(time, datetime):
            time = [time, time]    
        #compute timegaps between 2 particles
        dt = (time[1]-time[0]).total_seconds()/float(number)
        #generate a datetime object for each release time
        dts = dt*np.arange(number)
        times = [(time[0] + timedelta(seconds=i)) for i in dts]
    
    
    if plot is True:
        import matplotlib.pyplot as plt
        plt.plot(lons, lats, '.')
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        if time[0] == time[-1]:
            plt.title(str(time[0]))
        else:
            plt.title('%s - %s' % (time[0], time[-1]))
        plt.show() 

    return lons, lats, depths, times # Arrays of length *number*


def noos_seed_from_json(json_source=None, plot=False):
    '''Return arrays of lon, lat and time based on json obsect
    
    json_source may be either filename or a dictionary
    '''

    if isinstance(json_source, str):
        file = open(json_source)
        seed = json.load(file)
    elif isinstance(json_source, dict):
        seed = json_source
    if 'initial_condition' in seed:
        seed = seed['initial_condition']

    t = seed['time']
    if isinstance(t, list) and len(t) == 2:
        time = [datetime.strptime(t[0], '%Y-%m-%dT%H:%M:%SZ'),
                datetime.strptime(t[1], '%Y-%m-%dT%H:%M:%SZ')]
    elif isinstance(t, list) and len(t) == 1:
        time = [datetime.strptime(t[0], '%Y-%m-%dT%H:%M:%SZ'),
                datetime.strptime(t[0], '%Y-%m-%dT%H:%M:%SZ')]
    else:
        time = datetime.strptime(t, '%Y-%m-%dT%H:%M:%SZ')

    try :
      z = seed['z']
    except:    
      z = None


    lons, lats, depths, times = noos_seed(seed['lon'], seed['lat'],
                                 number=seed['number'],
                                 radius=seed['radius'],
                                 z=z,
                                 time=time,
                                 plot=plot)

    return lons, lats, depths, times
