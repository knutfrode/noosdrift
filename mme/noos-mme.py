#!/usr/bin/env python

import os
import glob
import json
import pyproj
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import patches
from netCDF4 import Dataset, num2date

matplotlib.colors.colorConverter.to_rgba('mediumseagreen', alpha=.1)

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
    angle = np.degrees(np.math.atan2(eigvecs[1][major],
                                     eigvecs[1][1-major]))
    if major == 1:
        angle = -angle
    return major_axis, minor_axis, angle


# To encode JSON
class MyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.floating):
            return np.round(float(obj), 4)
        else:
            return super(MyEncoder, self).default(obj)


class Simulation():

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
        self.azimuth[0] = np.nan
        self.distance = np.array(self.distance)

    def get_ellipses(self):
        self.major_axis = np.ones(self.num_timesteps)
        self.minor_axis = np.ones(self.num_timesteps)
        self.angle = np.ones(self.num_timesteps)
        for i in range(self.num_timesteps):
            self.major_axis[i], self.minor_axis[i], self.angle[i] = \
                get_ellipse(self.lon[:,i], self.lat[:,i])
        self.area = np.pi*self.major_axis*self.minor_axis/4

    def plot_timestep(self, ax, i):
        x = self.x[:,i]
        y = self.y[:,i]
        ax.plot(x, y, '.', color=self.color,
                alpha=.1, markeredgewidth=0)
        #for val, vec in zip(eigvals, eigvecs.T):
        #    val*=2
        #    val = np.sqrt(val)
        #    xv, yv = np.vstack((center + val*vec, center,
        #                        center - val*vec)).T
        #    ax.plot(xv,yv, 'b-', lw=3)
        ellipse = patches.Ellipse(
            (x.mean(), y.mean()),
            3*self.major_axis[i], 3*self.minor_axis[i],
            angle=self.angle[i], linewidth=2, fill=False,
            zorder=10, color=self.color)
        #import pdb; pdb.set_trace()
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

    def get_points_geojson(self, filename=None):
        '''Return gejson dictionary with all particles for each timestep'''
        
        pg = {'type': 'FeatureCollection',
              'properties': self.json_summary(),
              'features': []}

        for i in range(len(self.times)):
            coords = [ [lo,la] for lo,la in
                       zip(self.lon[:,i], self.lat[:,i]) ]
            pg['features'].append({'type': 'Feature',
                'properties': {'time':
                                self.times[i].isoformat('T')+'Z'},
                'geometry': {'type': 'MultiPoint',
                             'coordinates': coords}})

        if filename is not None:
            with open(filename, 'w') as outfile:
                json.dump(pg, outfile, cls=MyEncoder, indent=2)

        return pg


    def get_analysis_geojson(self, filename=None):

        pg = {'type': 'FeatureCollection',
              'properties': self.json_summary(),
              'features': []}

        for i in range(len(self.times)):
            lon = self.lon[:,i]
            lat = self.lat[:,i]
            localproj = pyproj.Proj(
                '+proj=stere +lat_0=%f +lat_ts=%f +lon_0=%s' %
                (lat.mean(), lat.mean(), lon.mean()))
            x, y = localproj(lon, lat, inverse=False)
            major_axis, minor_axis, angle = \
                get_ellipse(lon, lat)

            pg['features'].append({
                'time': self.times[i].isoformat('T')+'Z',
                'ellipsis_major_axis': np.round(major_axis, 2),
                'ellipsis_minor_axis': np.round(minor_axis, 2),
                'ellipsis_major_axis_azimuth_angle': np.round(angle, 2),
                'distance_of_center_from_start': np.round(self.distance[i], 2),
                'azimuth_direction_of_center_from_start':
                    np.round(self.azimuth[i], 2)})

        if filename is not None:
            with open(filename, 'w') as outfile:
                json.dump(pg, outfile, cls=MyEncoder, indent=2)

        return pg

    def __repr__(self):
        return 'Simulation: ' + self.filename


class SimulationCollection():

    colors = ['r', 'g', 'b', 'm', 'y', 'k']

    def __init__(self, *simulations):

        simulations = [Simulation(s) if isinstance(s, str) else s for s in simulations[0]]
            
        self.simulations = simulations

        self.lonmin = min([s.lonmin for s in self.simulations])
        self.lonmax = min([s.lonmax for s in self.simulations])
        self.latmin = min([s.latmin for s in self.simulations])
        self.latmax = min([s.latmin for s in self.simulations])
        self.centerlon = np.mean([s.lon[:,0].mean() for s in self.simulations])
        self.centerlat = np.mean([s.lat[:,0].mean() for s in self.simulations])
        self.proj = pyproj.Proj(
            '+proj=stere +lat_0=%f +lat_ts=%f +lon_0=%s' %
            (self.centerlat, self.centerlat, self.centerlon))
        for i, s in enumerate(self.simulations):
            s.color = self.colors[i]  
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

    def plot(self):
        for i in range(5, self.num_timesteps, 8):
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

    def write_geojson(self, filename):
        '''Obsolete, as this will be plotted from axes and angle'''
        for s in self.simulations:
            print(s)
            print(dir(s))
            print(s.times)
            for i in range(0, 3):
                print(s.times[i])
                lon = s.lon[:,i]
                lat = s.lat[:,i]
                localproj = pyproj.Proj(
                    '+proj=stere +lat_0=%f +lat_ts=%f +lon_0=%s' %
                    (lat.mean(), lat.mean(), lon.mean()))
                x, y = localproj(lon, lat, inverse=False)
                #major_axis, minor_axis, angle, localproj = \
                #    get_ellipse(lon, lat)
                #s.get_ellipses()
                print(s.major_axis)
                ellipse = patches.Ellipse((x.mean(), y.mean()),
                    3*s.major_axis[i], 3*s.minor_axis[i],
                    angle=s.angle[i], linewidth=2, fill=False,
                    zorder=10, color=s.color)
                V = ellipse.get_verts()
                X = V[:,0]
                Y = V[:,1]
                lonellipse,latellipse = localproj(X, Y, inverse=True)
                plt.subplot(2,1,1)
                plt.plot(lonellipse, latellipse, '*')
                plt.subplot(2,1,2)
                plt.plot(X, Y, '*')
                plt.show()
                print(X)
                print(len(X))
                print(latellipse)

                polygonStr = '{"type": "Feature", "properties": { "model": "OpenDrift", "time": "20190619T120000Z", "area": 25012.2}, "geometry": { "type": "Polygon", "coordinates": [ ['

                for p in range(len(X)):
                    polygonStr = polygonStr + '[%f,%f],' % (
                        lonellipse[p], latellipse[p])
                polygonStr = polygonStr + '[%f,%f]]]}}' % (lonellipse[0], latellipse[0])

                print(polygonStr)
            stop

if __name__ == '__main__':


    ## Usage demonstration
    #simulation_files = glob.glob('./sample_simulations/*.nc')
    #s = SimulationCollection(simulation_files)
    #print s
    #s.plot()
    #s.plot_metrics()
    #stop
    
    s1 = Simulation('sample_simulations/opendrift_oil_norway_rlw.nc')
    s1.get_points_geojson(filename='noos_points.json')
    s1.get_analysis_geojson(filename='noos_analysis.json')

    #s2 = Simulation('mothy_oil_france_rhw.nc')
    #s3 = Simulation('mothy_oil_france_rhw_without_nws.nc')
    #s = SimulationCollection(s1, s2, s3)

    ##s1 = Simulation('opendrift_oil_norway_rlw.nc')
    ##s2 = Simulation('mothy_oil_norway_rlw.nc')
    ##s3 = Simulation('mothy_oil_norway_rlw_with_ibi.nc')
    ##s4 = Simulation('mothy_oil_norway_rlw_with_psy4.nc')
    ##s = SimulationCollection(s1, s2, s3, s4)
    #print(s.simulations)
    ##s.plot()
    #s.plot_metrics()

#    for sim in ['oil', 'sar']:
#        for loc in ['france', 'belgium', 'norway']:
#            for cas in ['rlw', 'rhw']:
#                if cas != 'rlw' or sim != 'oil' or loc != 'norway':
#                    continue
#                files = glob.glob('../*%s_%s_%s*.nc' % (sim, loc, cas))
#                files = [files[0], files[2]]  # Subset
#                print(files)
#                s = SimulationCollection(files)
#                s.write_geojson('geojson.txt')
#                stop
#                #s.plot()
#                #s.plot_metrics()
#                #stop
#                #from opendrift.models.openoil3D import OpenOil3D
#                #o = opendrift.open(files[0])
#                #label = [si.label for si in s.simulations]
#                #o = OpenOil3D()
#                #o.io_import_file(files[0])
#                #o.animation(compare=files[1:], legend=label,
#                #            filename='anim_%s_%s_%s.mp4' %
#                #                (sim, loc, cas))
