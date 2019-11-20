#!/usr/bin/env python

# Multi-Model-Ensemble (MME) analysis tool 
# developed within the Noos-Drift project
# by Knut-Frode Dagestad (MET Norway) 
# Sept 2019

import matplotlib
matplotlib.use("TkAgg")

import sys
import argparse
import os
import glob
from math import radians, cos, sin, asin, sqrt
import json
import pyproj
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
from netCDF4 import Dataset, num2date
from sklearn.cluster import MeanShift


colors = ['g', 'b', 'm', 'y', 'c', 'saddlebrown', 'deeppink',
          'coral', 'darkgrey', 'lime']


def process_folder(inputfolder, outputfolder=None, animation=False):
    '''Perform MME analysis of all simulations stored in a folder

    :param inputfolder:
    :param outputfolder:
    :param animation:
    :return:
    '''

    simulation_files = glob.glob(inputfolder + '/*.nc')
    s = SimulationCollection(simulation_files)

    try:
        s.noosID = os.path.basename(simulation_files[0].split('/')[-1]).split('_')[1]
    except:
        s.noosID = 'requestID'

    if outputfolder is None:
        outputfolder = inputfolder
    mmefile = outputfolder + '/noosdrift_%s.json' % (s.noosID)

    print('Writing MME-analysis to file: %s' % mmefile)
    s.mme_analysis(mmefile)

    if animation is not False:
        animfile = outputfolder + '/noosdrift_%s.mp4' % (s.noosID)
        animate_mme_analysis(filename=mmefile, animfile=animfile, simulationcollection=s)

    # Produce JSON file for each simulation
    print('Writing point JSON files for each simulation')
    for sim in s.simulations:
        print('\t' + sim.filename)
        sim.write_point_geojson(filename =
                                outputfolder + '/noosdrift_%s_%s_%s_%s.json' % (
                                s.noosID, sim.model, sim.current, sim.wind))

def animate_mme_analysis(filename, animfile=None, simulationcollection=None):
    '''
    :param filename:
    :param animfile:
    :param simulationcollection:
    :return:
    '''

    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        import matplotlib.patches as mpatches
        from matplotlib import animation
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

    fig = plt.figure(figsize=(10., 10.))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mercator(
                         central_longitude=centerlon))
    ax.gridlines(draw_labels=True)
    buffer = .1
    ax.set_extent([lonmin-buffer, lonmax+buffer,
                   latmin-buffer, latmax+buffer],
                  crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.GSHHSFeature('high', edgecolor='black',
                   facecolor=cfeature.COLORS['land']))

    def plot_timestep(i):
        '''
        :param i:
        :return:
        '''

        print(i)
        feature = j['features'][i]
        clusters = feature['clusters']
        ellipses = feature['ellipses']
        for sn, si in enumerate(sim_list):
            # Points
            si['particles'].set_offsets(np.c_[
                si['lon'][range(si['lon'].shape[0]), i],
                si['lat'][range(si['lat'].shape[0]), i]])
            # Ellipses
            f = ellipses[str(sn)]
            xy = ccrs.Mercator().transform_points(
                ccrs.Geodetic(), np.array(f['centerlon']), np.array(f['centerlat']))[0]
            si['ellipse'].set_center((xy[0], xy[1]))
            si['ellipse'].height = f['ellipsis_major_axis'] * 8
            si['ellipse'].width = f['ellipsis_minor_axis'] * 8
            si['ellipse'].angle = -f['ellipsis_major_axis_azimuth_angle']
            # Clusters
            #if sn < len(clusters):
            if sn < len(clusters):
                f = feature['clusters'][str(sn)]
                xy = ccrs.Mercator().transform_points(
                    ccrs.Geodetic(), np.array(f['centerlon']), np.array(f['centerlat']))[0]
                radius = 3 * (f['longest_ellipsis_axis'] +
                              np.mean(f['distance_from_cluster_centre']))
                if len(f['distance_from_cluster_centre']) <= 1:
                    xy = (0, 0)
                    radius = 0
            else:
                radius = 0
                xy = (0, 0)
            si['cluster'].set_center((xy[0], xy[1]))
            si['cluster'].radius = radius

        # Super-ellipse
        f = feature['super-ellipse']
        xy = ccrs.Mercator().transform_points(
            ccrs.Geodetic(), np.array(f['centerlon']), np.array(f['centerlat']))[0]
        si['super-ellipse'].set_center((xy[0], xy[1]))
        si['super-ellipse'].height = f['ellipsis_major_axis'] * 8
        si['super-ellipse'].width = f['ellipsis_minor_axis'] * 8
        si['super-ellipse'].angle = -f['ellipsis_major_axis_azimuth_angle']

        return si['particles'], si['ellipse'], si['super-ellipse'], si['cluster']

    if simulationcollection is not None:
        num_sim = len(simulationcollection.simulations)
        sim_list = [{}] * num_sim

        for sn, sim in enumerate(simulationcollection.simulations):
            if sn == 0:
                ax.scatter(sim.lon[:, 0], sim.lat[:, 0], color='k', zorder=150,
                           transform=ccrs.Geodetic(), label='initial position')
            sim_list[sn] = {}
            sd = sim_list[sn]  # pointer
            sd['name'] = sim.label
            sd['lon'] = sim.lon  #.copy()
            sd['lat'] = sim.lat  #.copy()
            # Plot points
            sd['particles'] = ax.scatter(
                [], [], color=sim.color,
                transform=ccrs.Geodetic(),
                label=sim.label, zorder=90)
            # Plot bounding ellipse
            sd['ellipse'] = ax.add_patch(mpatches.Ellipse(
                xy=[0, 0], height=0, width=0, angle=0, fill=False,
                color='k', alpha=1, transform=ccrs.Mercator(), zorder=100))
            # Plot cluster ellipse
            sd['cluster'] = ax.add_patch(mpatches.Circle(
                xy=[0, 0], radius=0, fill=False,
                color='r', alpha=1, transform=ccrs.Mercator(), zorder=120))

        # Plot super ellipse
        sd['super-ellipse'] = ax.add_patch(mpatches.Ellipse(
            xy=[0, 0], height=0, width=0, angle=0, fill=False,
            linewidth=3,
            color='k', alpha=1, transform=ccrs.Mercator(), zorder=200))

    anim = animation.FuncAnimation(plt.gcf(), plot_timestep, blit=False,
                                   frames=len(times), interval=50)

    plt.legend(loc='upper right')

    if animfile is None:
        plt.show()
    else:
        print('Saving animation to: ' + animfile)
        FFwriter = animation.FFMpegWriter(
            fps=5, codec='libx264', bitrate=1800,
            extra_args=['-profile:v', 'baseline', '-pix_fmt', 'yuv420p', '-an'])
        anim.save(animfile, writer=FFwriter)

def plot_mme_analysis(filename, simulationcollection=None):
    '''
    Import and plot the contents of a MME output JSON file

    :param filename:
    :param simulationcollection:
    :return:
    '''

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
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mercator(
                    central_longitude=centerlon))
        ax.gridlines(draw_labels=True)
        buffer = .1
        ax.set_extent([lonmin - buffer, lonmax + buffer,
                       latmin - buffer, latmax + buffer],
                      crs=ccrs.Geodetic())
        # Draw coastlines.
        ax.add_feature(cfeature.NaturalEarthFeature(
            'physical', 'land', '50m', edgecolor='black',
            facecolor=cfeature.COLORS['land']))
        #ax.add_feature(cfeature.GSHHSFeature('high', edgecolor='black',
        #                                     facecolor=cfeature.COLORS['land']))
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
                height=major_axis * 10, width=minor_axis * 10,
                angle=-angle, fill=False,
                color='k', alpha=1,
                transform=ccrs.Mercator(), zorder=30))
            ax.plot(lon, lat, '*k', zorder=50,
                    transform=ccrs.Geodetic())

        for cluster in f['clusters']:
            c = f['clusters'][cluster]
            if len(c['members']) > 1:
                # Plot each cluster, if more than one member
                xy = ccrs.Mercator().transform_points(
                        ccrs.PlateCarree(),
                        np.array(c['centerlon']),
                        np.array(c['centerlat']))[0]
                ax.add_patch(mpatches.Circle(
                    xy=xy, radius = 3 * (c['longest_ellipsis_axis'] + 
                                         np.mean(c['distance_from_cluster_centre'])),
                    fill=False, color='r', lw=2,
                    transform=ccrs.Mercator(), zorder=50))

        # Plot embounding super-ellipse
        se = f['super-ellipse']
        lon = np.array(se['centerlon'])
        lat = np.array(se['centerlat'])
        minor_axis = se['ellipsis_minor_axis']
        major_axis = se['ellipsis_major_axis']
        angle = se['ellipsis_major_axis_azimuth_angle']
        xy = ccrs.Mercator().transform_points(
            ccrs.PlateCarree(), lon, lat)[0]
        ax.add_patch(mpatches.Ellipse(
            xy=[xy[0], xy[1]],
            height=major_axis*10, width=minor_axis*10,
            angle=-angle, fill=False,
            linewidth=3,
            color='k', alpha=1,
            transform=ccrs.Mercator(), zorder=30))
        ax.plot(lon, lat, '*k', zorder=50,
                    transform=ccrs.Geodetic())

        if simulationcollection is not None:
            s = simulationcollection.simulations[0]
            ax.plot(s.lon[:,0], s.lat[:,0], 'y.',
                    transform=ccrs.Geodetic(),
                    label='initial position')
            for s in simulationcollection.simulations:
                ax.plot(s.lon[:,ti], s.lat[:,ti], '.',
                        transform=ccrs.Geodetic(),
                        label=s.label)
            
        plt.legend()
        plt.show()

def haversine(lon1, lat1, lon2, lat2):
    '''
    Calculate the great circle distance between two points

    :param lon1:
    :param lat1:
    :param lon2:
    :param lat2:
    :return:
    '''
    
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371000 # Radius of earth in m
    return c * r

def get_ellipse(lons, lats):
    '''
    Calculate best-fit ellipse for a cloud of lon,lat positions
    
    :param lons:
    :param lats:
    :return:
    '''

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
    '''
    To encode JSON

    '''

    def default(self, obj):
        '''
        
        :param obj:
        :return:
        '''

        if isinstance(obj, np.floating):
            return np.round(float(obj), 4)
        else:
            return super(MyEncoder, self).default(obj)


class Simulation():
    '''
    Containing data for a single simulation

    '''

    def __init__(self, filename):
        '''

        :param filename:
        '''

        print('Importing: ' + filename)

        self.filename = os.path.basename(filename)
        infile = Dataset(filename, 'r')
        self.noosID = self.filename.split('_')[1]

        # Time
        self.times = num2date(infile.variables['time'][:],
                              infile.variables['time'].units)
        self.start_time = self.times[0]
        self.end_time = self.times[-1]
        self.time_step = self.times[1] - self.times[0]
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
        shortname = os.path.basename(self.filename)
        shortname = os.path.splitext(shortname)[0]
        self.model = shortname.split('_')[2].lower()
        self.current = 'Unknown'
        self.current = shortname.split('_')[3].lower()
        #if 'without_nws' in shortname:
        #    self.current = 'tidal'
        #if 'with_ibi' in shortname:
        #    self.current = 'cmems-ibi'
        #if 'with_psy4' in shortname:
        #    self.current = 'cmems-global'
        self.wind = 'Unknown'
        self.wind = shortname.split('_')[4].lower()
        # Override with names from netCDF file, if existing
        if 'current_name' in attributes:
            self.current = attributes['current_name']
        if 'wind_name' in attributes:
            self.wind = attributes['wind_name']

        if self.wind == 'ECMWF':
            self.wind = 'ecmwf'
        if self.current == 'Topaz':
            self.current = 'topaz'

        #self.label = '{0:<12}{1:<8}{2:<8}'.format(
        #            self.model, self.current, self.wind)
        self.label = '%-15s%-15s%-10s' % (self.model, self.current, self.wind)

        geod = pyproj.Geod(ellps='WGS84')
        res = [geod.inv(self.centerlon[0], self.centerlat[0],
            lon, lat) for lon, lat in zip(self.centerlon, self.centerlat)]
        self.azimuth, backaz, self.distance = zip(*res)
        self.azimuth = np.array(self.azimuth)
        self.azimuth[0] = 0
        self.distance = np.array(self.distance)

    def get_ellipses(self):
        '''
        Get best-fit-ellipsis for each timestep of a simulation

        :return:
        '''

        self.major_axis = np.ones(self.num_timesteps)
        self.minor_axis = np.ones(self.num_timesteps)
        self.angle = np.ones(self.num_timesteps)
        for i in range(self.num_timesteps):
            self.major_axis[i], self.minor_axis[i], self.angle[i] = \
                get_ellipse(self.lon[:,i], self.lat[:,i])
        self.area = np.pi*self.major_axis*self.minor_axis/4

    def plot_timestep(self, ax, i):
        '''
        Plot timestep of a simulation

        :param ax:
        :param i:
        :return:
        '''

        x = self.x[:,i]
        y = self.y[:,i]
        ax.plot(x, y, '.', color=self.color,
                alpha=.1, markeredgewidth=0)
        ellipse = patches.Ellipse(
            (x.mean(), y.mean()),
            3 * self.minor_axis[i], 3 * self.major_axis[i],
            angle=-self.angle[i], linewidth=2, fill=False,
            zorder=10, color=self.color)
        ax.add_patch(ellipse)
        ax.plot([self.x[:, 0].mean(), x.mean()],
                [self.y[:, 0].mean(), y.mean()], color=self.color,
                label=self.label)
        ax.legend()

    def json_summary(self):
        '''
        Return some common JSON properties
        :return:        
        '''

        s = {'time_step': '%sH' % (self.time_step.seconds/3600.),
             'start_time': self.times[0].isoformat('T')+'Z',
             'end_time': self.times[-1].isoformat('T')+'Z',
             'number_of_times': len(self.times)}
        return s

    def write_point_geojson(self, filename=None):
        '''
        
        :param filename:
        :return:
        '''

        pg = {'type': 'FeatureCollection',
              'properties': self.json_summary(),
              'features': []}

        # Temporarily add hardoded forcing names
        pg['properties']['model_name'] = self.model
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
            major_axis, minor_axis, angle = get_ellipse(lon, lat)

            # individual points
            coords = [ [lo,la] for lo,la in
                       zip(self.lon[:, i], self.lat[:, i]) ]
            # Removing invalid coordinates (masked values in "lat" & "lon").
            # Replacing masked value with 0 will likely give strange results
            coords[:] = [coord for coord in coords if not (coord[0] > 180 or coord[1] > 90)]

            pg['features'].append({
                'time': self.times[i].isoformat('T') + 'Z',
                'latitude_of_center': np.mean(lat),
                'longitude_of_center': np.mean(lon),
                'ellipsis_major_axis': major_axis,
                'ellipsis_minor_axis': minor_axis,
                'ellipsis_major_axis_azimuth_angle': angle,
                'distance_of_center_from_start': self.distance[i],
                'azimuth_direction_of_center_from_start':
                    self.azimuth[i],
                'geometry': {'type': 'MultiPoint',
                             'coordinates': coords}})

        if filename is not None:
            with open(filename, 'w') as outfile:
                json.dump(pg, outfile, cls=MyEncoder, indent=2)

        major_axis = [f['ellipsis_major_axis']
                      for f in pg['features']]

        return pg, major_axis, minor_axis

    def __repr__(self):
        '''
        
        :return:
        '''

        return 'Simulation: ' + self.filename


class SimulationCollection():
    '''
    Contains a collection of individual trajectory simulations

    '''

    def __init__(self, *simulations):
        '''
        
        :param simulations:
        '''

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
            s.x = np.ma.masked_where(s.lon.mask == 1, s.x)
            s.y = np.ma.masked_where(s.lon.mask == 1, s.y)
            s.get_ellipses()
        self.time = s.times  # use time of last simulation
        self.hours = [(t - self.time[0]).total_seconds()/3600 for t in self.time]
        self.xmin = min([s.x.min() for s in self.simulations])
        self.xmax = max([s.x.max() for s in self.simulations])
        self.ymin = min([s.y.min() for s in self.simulations])
        self.ymax = max([s.y.max() for s in self.simulations])
        self.centerx = np.mean([np.mean(s.x) for s in self.simulations])
        self.centery = np.mean([np.mean(s.y) for s in self.simulations])

        self.num_timesteps = self.simulations[0].num_timesteps

    def mme_analysis(self, outfile):
        '''
        Perform multi-model-ensemble analysis of simulations

        :param outfile:
        :return:
        '''

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
        for i, s in enumerate(self.simulations):
            pg['simulations'][i] = '%s_%s_%s' % (
                                    s.model.lower(), s.current, s.wind)
        
        # Adding then some data for each timestep
        for i in range(0, self.num_timesteps, 1):
            tf = {
                'time': s.times[i].isoformat('T')+'Z',
                'super-ellipse': {},
                'ellipses': {},
                'clusters': {}
                }

            # Find clusters
            X = [ [s.centerlon[i], s.centerlat[i]]
                    for s in self.simulations]
            X = np.array(X)
            hours = (s.times[i]-s.times[0]).total_seconds()/3600.
            distances = [haversine(s.centerlon[0], s.centerlat[0],
                                   s.centerlon[i], s.centerlat[i])
                         for s in self.simulations]
            avg_distance = np.mean(distances)
            bandwidth = .10 - .02*(avg_distance/30000.)
            ms = MeanShift(bandwidth=bandwidth)
            ms.fit(X)
            labels = ms.labels_
            cluster_centers = ms.cluster_centers_
            n_clusters_ = len(np.unique(labels))

            all_lons = np.array([])
            all_lats = np.array([])
            for c, s in enumerate(self.simulations):
                all_lons = np.concatenate((all_lons, s.lon[:,i]))
                all_lats = np.concatenate((all_lats, s.lat[:,i]))
                tf['ellipses'][c] = {
                    'ellipsis_major_axis': s.major_axis[i],
                    'ellipsis_minor_axis': s.minor_axis[i],
                    'centerlon': s.centerlon[i],
                    'centerlat': s.centerlat[i],
                    'ellipsis_major_axis_azimuth_angle': s.angle[i]}
            all_major_axis, all_minor_axis, all_angle = \
                get_ellipse(all_lons, all_lats)
            tf['super-ellipse']['centerlon'] = np.mean(all_lons)
            tf['super-ellipse']['centerlat'] = np.mean(all_lats)
            tf['super-ellipse']['ellipsis_major_axis'] = all_major_axis
            tf['super-ellipse']['ellipsis_minor_axis'] = all_minor_axis
            tf['super-ellipse']['ellipsis_major_axis_azimuth_angle'] = all_angle
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
                tf['clusters'][c]['distance_from_cluster_centre'] = distance_from_centre
                tf['clusters'][c]['distance_std'] = distance_std
                longest_axis = 0
                for j in members:
                    longest_axis = np.maximum(longest_axis,
                        tf['ellipses'][j]['ellipsis_major_axis'])
                tf['clusters'][c]['longest_ellipsis_axis'] = longest_axis

            pg['features'].append(tf)

        # Write output of MME analysis to JSON file
        if outfile is not None:
            with open(outfile, 'w') as of:
                json.dump(pg, of, cls=MyEncoder, indent=2)

    def plot(self):
        '''
        Plotting a simulation collection

        :return:
        '''

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
        '''
        Plotting some scalar metrics of a simulation collection

        :return:
        '''

        fig, (axdist, axazimuth, axarea) = plt.subplots(3)

        for s in self.simulations:
            axdist.plot(self.hours, s.distance / 1000.,
                    s.color, label=s.label)
            axazimuth.plot(self.hours, s.azimuth, s.color)
            axarea.plot(self.hours, s.area / 1e6, s.color)

        axdist.set_ylabel('Distance [km]')
        axazimuth.set_ylabel('Direction [deg azimuth]')
        axarea.set_ylabel('Area [km2]')
        axdist.legend()

        for ax in (axdist, axazimuth, axarea):
            ax.set_xlabel('Time [hours]')
            ax.set_xlim([0, self.hours[-1]])

        plt.show()

    def __repr__(self):
        '''

        :return:
        '''

        r = 'NoosID: %s\n' % self.noosID

        r = r + '%12s%12s%12s\n' % ('Model', 'Current', 'Wind')
        r = r + '-' * 36 + '\n'
        for s in self.simulations:
            r = r + '%12s%12s%12s\n' % (s.model, s.current, s.wind)

        return r


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='inputfolder',
                        default='sample_simulations',
                        help='Folder with netCDF files from simulations')
    parser.add_argument('-o', dest='outputfolder',
                        default=None,
                        help='Folder with MME output')
    parser.add_argument('-f', '--write_animation', action="store_true",
                        help='Write animation (mp4) to output folder')

    args = parser.parse_args()

    print('Input folder: ' + args.inputfolder)
    if args.outputfolder is None:
        args.outputfolder = args.inputfolder + '/mme_output/'
    print('Output folder: ' + args.outputfolder)

    if not (os.path.exists(args.inputfolder)):
        sys.exit('Input folder does not exist: ' + args.inputfolder)
    if not (os.path.exists(args.outputfolder)):
        print('Output folder does not exist, creating: ' +
              args.outputfolder)
        try:
            os.mkdir(args.outputfolder)
        except:
            sys.exit('Could not create output folder: ' +
                     args.outputfolder)

    process_folder(args.inputfolder, args.outputfolder, args.write_animation)
