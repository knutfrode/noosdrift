#!/usr/bin/env python

# Multi-Model-Ensemble (MME) analysis tool 
# developed within the Noos-Drift project
# by Knut-Frode Dagestad (MET Norway) 
# Sept 2019

import sys
import argparse
import os
import glob
from math import radians, cos, sin, asin, sqrt
import json
import pyproj
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import patches
from netCDF4 import Dataset, num2date
from sklearn.cluster import MeanShift

colors = ['r', 'g', 'b', 'm', 'y', 'k']


def process_folder(inputfolder, outputfolder, noosID='requestID'):
    '''Perform MME analysis of all simulations stored in a folder'''

    simulation_files = glob.glob(inputfolder + '/*.nc')
    s = SimulationCollection(simulation_files)
    s.noosID = noosID

    print('Performing MME-analysis')
    s.mme_analysis(outfile='mme.json')

    plot_mme_analysis(filename='mme.json', simulationcollection=s)

    # Produce JSON file for each simulation
    print('Writing point JSON files for each simulation')
    for sim in s.simulations:
        sim.write_point_geojson(filename = outputfolder +
            '/noosdrift_%s_%s_%s_%s.json' % (
            noosID, sim.model, sim.current, sim.wind))

def plot_mme_analysis(filename, simulationcollection=None):
    '''Import and plot the contents of a MME output JSON file'''

    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        import matplotlib.patches as mpatches
    except:
        raise ImportError('Please install Cartopy to make plots')

    j = json.load(open(filename, 'r'))
    lonmin = j['coverage']['lonmin']
    lonmax = j['coverage']['lonmax']
    latmin = j['coverage']['latmin']
    latmax = j['coverage']['latmax']
    centerlon = j['coverage']['centerlon']
    centerlat = j['coverage']['centerlat']

    times = np.arange(0, len(j['features']))
    for ti in times:
        f = j['features'][ti]
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, projection=ccrs.Mercator(
                    central_longitude=centerlon))
        ax.gridlines(draw_labels=True)
        buffer = .1
        ax.set_extent([lonmin-buffer, lonmax+buffer,
                       latmin-buffer, latmax+buffer],
                      crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.NaturalEarthFeature(
            'physical', 'land', '110m', edgecolor='black',
            facecolor=cfeature.COLORS['land']))
        for num, e in enumerate(f['ellipses']):
            el = f['ellipses'][e]
            lon = np.array(el['centerlon'])
            lat = np.array(el['centerlat'])
            minor_axis = el['ellipsis_minor_axis']
            major_axis = el['ellipsis_major_axis']
            angle = el['ellipsis_major_axis_azimuth_angle']
            xy = ccrs.Mercator().transform_points(
                ccrs.PlateCarree(), lon, lat)[0]
            ax.add_patch(mpatches.Ellipse(
                xy=[xy[0], xy[1]],
                height=major_axis*10, width=minor_axis*10,
                angle=-angle, fill=False,
                color='k', alpha=1,
                transform=ccrs.Mercator(), zorder=30))
            ax.plot(lon, lat, '*k', zorder=50,
                    transform=ccrs.Geodetic())

        for cluster in f['clusters']:
            c = f['clusters'][cluster]
            if len(c['members'])>1:
                # Plot each cluster, if more than one member
                xy = ccrs.Mercator().transform_points(
                        ccrs.PlateCarree(),
                        np.array(c['centerlon']),
                        np.array(c['centerlat']))[0]
                ax.add_patch(mpatches.Circle(
                    xy=xy, radius= 3*(c['longest_ellipsis_axis'] + 
                        np.mean(c['distance_from_cluster_centre'])),
                    fill=False, color='r', lw=2,
                    transform=ccrs.Mercator(), zorder=50))

        if simulationcollection is not None:
            for s in simulationcollection.simulations:
                ax.plot(s.lon[:,ti], s.lat[:,ti], '.',
                        transform=ccrs.Geodetic())
            
        plt.show()

def haversine(lon1, lat1, lon2, lat2):
    ''' Calculate the great circle distance between two points'''
    
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371000 # Radius of earth in m
    return c * r

def get_ellipse(lons, lats):
    '''Calculate best-fit ellipse for a cloud of lon,lat positions'''
    lons = lons[lons.mask==0]
    lats = lats[lats.mask==0]
    centerlon = np.mean(lons)
    centerlat = np.mean(lats)
    localproj = pyproj.Proj(
                    '+proj=stere +lat_0=%f +lat_ts=%f +lon_0=%s' %
                    (centerlat, centerlat, centerlon))
    x, y = localproj(lons, lats)
    centerx = x.mean()
    centery = y.mean()
    xy = np.row_stack((x, y))
    center = xy.mean(axis=-1)
    eigvals, eigvecs = np.linalg.eig(np.cov(xy))

    if eigvals[1] > eigvals[0]:
        major = 1
    else:
        major = 0
    major_axis = np.sqrt(eigvals)[major]
    minor_axis = np.sqrt(eigvals)[1-major]
    angle = 270 - np.degrees(np.math.atan2(eigvecs[1][major],
                                     eigvecs[1][1-major]))
    if major == 1:
        angle = -angle

    # Make azimuth angle between -180 and +180
    angle = angle % 360
    angle = (angle + 360) % 360
    if angle > 180:
        angle -= 360

    return major_axis, minor_axis, angle

# To encode JSON
class MyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.floating):
            return np.round(float(obj), 4)
        else:
            return super(MyEncoder, self).default(obj)


class Simulation():
    '''Containing data for a single simulation'''

    def __init__(self, filename):
        print('Importing: ' + filename)

        self.filename = os.path.basename(filename)
        infile = Dataset(filename, 'r')
        # Time
        self.times = num2date(infile.variables['time'][:],
                              infile.variables['time'].units)
        self.start_time = self.times[0]
        self.end_time = self.times[-1]
        self.time_step = self.times[1]-self.times[0]
        # Particles
        self.num_elements = len(infile.dimensions['trajectory'])
        self.num_timesteps = len(infile.dimensions['time'])

        self.lon = infile.variables['lon'][:]
        self.lat = infile.variables['lat'][:]
        self.lonmin = self.lon.min()
        self.lonmax = self.lon.max()
        self.latmin = self.lat.min()
        self.latmax = self.lat.max()
        self.centerlon = np.mean(self.lon, axis=0)
        self.centerlat = np.mean(self.lat, axis=0)

        #attributes = infile.ncattrs()
        attributes = infile.__dict__

        self.status = infile.variables['status'][:]
        self.status_categories = infile.variables['status'].flag_meanings.split()
        #if 'stranded' not in self.status_categories:
        #    self.num_stranded = np.zeros(len(self.time_step))
        #else:
        #    stranded_index = self.status_categories.index('stranded')*1
        #    self.num_stranded = np.sum(self.status==stranded_index, axis=0)
        infile.close()

        # Determining current and wind source from filename
        # Would be better to have attributes in netCDF files
        self.model = self.filename.split('_')[0].capitalize()
        self.current = 'NWS'
        if 'without_nws' in self.filename:
            self.current = 'TIDAL'
        if 'with_ibi' in self.filename:
            self.current = 'IBER'
        if 'with_psy4' in self.filename:
            self.current = 'MERC'
        self.wind = 'ECMWF'
        # Override with names from netCDF file, if existing
        if 'current_name' in attributes:
            self.current = attributes['current_name']
        if 'wind_name' in attributes:
            self.wind = attributes['wind_name']
        self.label = '{0:<12}{1:<8}{2:<8}'.format(
                    self.model, self.current, self.wind)

        geod = pyproj.Geod(ellps='WGS84')
        res = [geod.inv(self.centerlon[0], self.centerlat[0],
            lon, lat) for lon, lat in zip(self.centerlon, self.centerlat)]
        self.azimuth, backaz, self.distance = zip(*res)
        self.azimuth = np.array(self.azimuth)
        self.azimuth[0] = 0
        self.distance = np.array(self.distance)

    def get_ellipses(self):
        '''Get best-fit-ellipsis for each timestep of a simulation'''
        self.major_axis = np.ones(self.num_timesteps)
        self.minor_axis = np.ones(self.num_timesteps)
        self.angle = np.ones(self.num_timesteps)
        for i in range(self.num_timesteps):
            self.major_axis[i], self.minor_axis[i], self.angle[i] = \
                get_ellipse(self.lon[:,i], self.lat[:,i])
        self.area = np.pi*self.major_axis*self.minor_axis/4

    def plot_timestep(self, ax, i):
        '''Plot timestep of a simulation'''
        x = self.x[:,i]
        y = self.y[:,i]
        ax.plot(x, y, '.', color=self.color,
                alpha=.1, markeredgewidth=0)
        ellipse = patches.Ellipse(
            (x.mean(), y.mean()),
            3*self.minor_axis[i], 3*self.major_axis[i],
            angle=-self.angle[i], linewidth=2, fill=False,
            zorder=10, color=self.color)
        ax.add_patch(ellipse)
        ax.plot([self.x[:,0].mean(), x.mean()],
                [self.y[:,0].mean(), y.mean()], color=self.color,
                label=self.label)
        ax.legend()

    def json_summary(self):
        '''Return some common JSON properties'''
        s = {'TimeStep': '%sH' % (self.time_step.seconds/3600.),
             'StartTime': self.times[0].isoformat('T')+'Z',
             'EndTime': self.times[-1].isoformat('T')+'Z',
             'number_of_times': len(self.times)}
        return s

    def write_point_geojson(self, filename=None):

        pg = {'type': 'FeatureCollection',
              'properties': self.json_summary(),
              'features': []}

        # Temporarily add hardoded forcing names
        pg['properties']['modelName'] = self.model
        try:
            pg['properties']['requestID'] = self.noosID
        except:
            pg['properties']['requestID'] = 'requestID'
        pg['properties']['wind_forcing'] = self.wind
        pg['properties']['ocean_forcing'] = self.current

        for i in range(len(self.times)):
            lon = self.lon[:,i]
            lat = self.lat[:,i]
            localproj = pyproj.Proj(
                '+proj=stere +lat_0=%f +lat_ts=%f +lon_0=%s' %
                (lat.mean(), lat.mean(), lon.mean()))
            x, y = localproj(lon, lat, inverse=False)
            major_axis, minor_axis, angle = \
                get_ellipse(lon, lat)

            # individual points
            coords = [ [lo,la] for lo,la in
                       zip(self.lon[:,i], self.lat[:,i]) ]

            pg['features'].append({
                'time': self.times[i].isoformat('T')+'Z',
                'latitude_of_center': np.mean(lat),
                'longitude_of_center': np.mean(lon),
                'ellipsis_major_axis': np.round(major_axis, 2),
                'ellipsis_minor_axis': np.round(minor_axis, 2),
                'ellipsis_major_axis_azimuth_angle': np.round(angle, 2),
                'distance_of_center_from_start': np.round(self.distance[i], 2),
                'azimuth_direction_of_center_from_start':
                    np.round(self.azimuth[i], 2),
                'geometry': {'type': 'MultiPoint',
                             'coordinates': coords}})

        if filename is not None:
            with open(filename, 'w') as outfile:
                json.dump(pg, outfile, cls=MyEncoder, indent=2)

        major_axis = [f['ellipsis_major_axis']
                        for f in pg['features']]

        return pg, major_axis, minor_axis

    def __repr__(self):
        return 'Simulation: ' + self.filename


class SimulationCollection():
    '''Contains a collection of individual trajectory simulations'''

    def __init__(self, *simulations):

        simulations = [Simulation(s) if isinstance(s, str) else s
                        for s in simulations[0]]
        self.simulations = simulations

        self.lonmin = min([s.lonmin for s in self.simulations])
        self.lonmax = max([s.lonmax for s in self.simulations])
        self.latmin = min([s.latmin for s in self.simulations])
        self.latmax = max([s.latmax for s in self.simulations])
        self.centerlon = np.mean([s.lon.mean()
            for s in self.simulations])
        self.centerlat = np.mean([s.lat.mean()
            for s in self.simulations])
        self.proj = pyproj.Proj(
            '+proj=stere +lat_0=%f +lat_ts=%f +lon_0=%s' %
            (self.centerlat, self.centerlat, self.centerlon))
        for i, s in enumerate(self.simulations):
            s.color = colors[i]
            s.x, s.y = self.proj(s.lon, s.lat)
            s.x = np.ma.masked_where(s.lon.mask==1, s.x)
            s.y = np.ma.masked_where(s.lon.mask==1, s.y)
            s.get_ellipses()
        self.time = s.times  # use time of last simulation
        self.hours = [(t-self.time[0]).total_seconds()/3600 for t in self.time]
        self.xmin = min([s.x.min() for s in self.simulations])
        self.xmax = max([s.x.max() for s in self.simulations])
        self.ymin = min([s.y.min() for s in self.simulations])
        self.ymax = max([s.y.max() for s in self.simulations])
        self.centerx = np.mean([np.mean(s.x) for s in self.simulations])
        self.centery = np.mean([np.mean(s.y) for s in self.simulations])

        self.num_timesteps = self.simulations[0].num_timesteps

    def mme_analysis(self, outfile):
        '''Perform multi-model-ensemble analysis of simulations'''

        # Dictionary to be written to JSON file
        # Adding first some overview data
        pg = {'type': 'FeatureCollection',
              'simulations': {},
              'coverage': {
                    'lonmin': self.lonmin,
                    'lonmax': self.lonmax,
                    'latmin': self.latmin,
                    'latmax': self.latmax,
                    'centerlon': self.centerlon,
                    'centerlat': self.centerlat,
                    },
              'features': []}

        # Generating a string reference for each model:
        # ModelName-CurrentName-WindName
        for i,s in enumerate(self.simulations):
            pg['simulations'][i] = '%s-%s-%s' % (
                                    s.model, s.current, s.wind)
        
        # Adding then some data for each timestep
        for i in range(0, self.num_timesteps, 1):
            tf = {
                'time': s.times[i].isoformat('T')+'Z',
                'ellipses': {},
                'clusters': {}
                }

            # Find clusters
            X = [ [s.centerlon[i], s.centerlat[i]]
                    for s in self.simulations]
            X = np.array(X)
            ms = MeanShift(bandwidth=.08)
            ms.fit(X)
            labels = ms.labels_
            cluster_centers = ms.cluster_centers_
            n_clusters_ = len(np.unique(labels))

            for c, s in enumerate(self.simulations):
                tf['ellipses'][c] = {
                    'ellipsis_major_axis': s.major_axis[i],
                    'ellipsis_minor_axis': s.minor_axis[i],
                    'centerlon': s.centerlon[i],
                    'centerlat': s.centerlat[i],
                    'ellipsis_major_axis_azimuth_angle': s.angle[i]}
            for c in range(n_clusters_):
                members = np.where(labels == c)[0].tolist()
                tf['clusters'][c] = {
                    'members': members,
                    'memberlons': [self.simulations[sn].centerlon[i]
                                    for sn in members],
                    'memberlats': [self.simulations[sn].centerlat[i]
                                    for sn in members],
                    'centerlon': cluster_centers[c][0],
                    'centerlat': cluster_centers[c][1]}

                distance_from_centre = [
                    haversine(tf['clusters'][c]['memberlons'][j],
                              tf['clusters'][c]['memberlats'][j],
                              tf['clusters'][c]['centerlon'],
                              tf['clusters'][c]['centerlat'])
                    for j in range(len(members))]
                distance_std = np.std(distance_from_centre)
                tf['clusters'][c]['distance_from_cluster_centre'] = \
                    distance_from_centre
                tf['clusters'][c]['distance_std'] = distance_std
                longest_axis = 0
                for j in members:
                    longest_axis = np.maximum(longest_axis,
                        tf['ellipses'][j]['ellipsis_major_axis'])
                tf['clusters'][c]['longest_ellipsis_axis'] = \
                    longest_axis

            pg['features'].append(tf)

        # Write output of MME analysis to JSON file
        if outfile is not None:
            with open(outfile, 'w') as of:
                json.dump(pg, of, cls=MyEncoder, indent=2)

    def plot(self):
        '''Plotting a simulation collection'''
        for i in range(0, self.num_timesteps, 1):
            fig, ax = plt.subplots() 
            ax.axis('equal')
            ax.set_xlabel('Distance [m]')
            ax.set_ylabel('Distance [m]')
            buffer = .2*(self.xmax-self.xmin)
            ax.set_xlim([self.xmin-buffer, self.xmax+buffer])
            ax.set_ylim([self.ymin-buffer, self.ymax+buffer])
            for s in self.simulations:
                s.plot_timestep(ax, i)
            plt.title('Duration: ' + str(s.times[i] - s.times[0]))
            plt.show()

    def plot_metrics(self):
        '''Plotting some scalar metrics of a simulation collection'''
        fig, (axdist, axazimuth, axarea) = plt.subplots(3)

        for s in self.simulations:
            axdist.plot(self.hours, s.distance/1000.,
                    s.color, label=s.label)
            axazimuth.plot(self.hours, s.azimuth, s.color)
            axarea.plot(self.hours, s.area/1e6, s.color)

        axdist.set_ylabel('Distance [km]')
        axazimuth.set_ylabel('Direction [deg azimuth]')
        axarea.set_ylabel('Area [km2]')
        axdist.legend()

        for ax in (axdist, axazimuth, axarea):
            ax.set_xlabel('Time [hours]')
            ax.set_xlim([0, self.hours[-1]])

        plt.show()

    def __repr__(self):
        r = 'NoosID: %s\n' % self.noosID

        r = r + '%12s%12s%12s\n' % ('Model', 'Current', 'Wind')
        r = r + '-'*36 + '\n'
        for s in self.simulations:
            r = r + '%12s%12s%12s\n' % (s.model, s.current, s.wind)

        return r


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='inputfolder',
                        default='sample_simulations',
                        help='Folder with netCDF files from simulations')
    parser.add_argument('-o', dest='outputfolder',
                        default='mme-output',
                        help='Folder with MME output')

    args = parser.parse_args()

    print('Input folder: ' + args.inputfolder)
    print('Output folder: ' + args.outputfolder)

    if not(os.path.exists(args.inputfolder)):
        sys.exit('Input folder does not exist: ' + args.inputfolder)
    if not(os.path.exists(args.outputfolder)):
        print('Output folder does not exist, creating: ' +
              args.outputfolder)
        try:
            os.mkdir(args.outputfolder)
        except:
            sys.exit('Could not create output folder: ' +
                     args.outputfolder)

    process_folder(args.inputfolder, args.outputfolder)
