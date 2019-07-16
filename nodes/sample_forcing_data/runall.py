#!/usr/bin/env python

from datetime import datetime, timedelta
from opendrift.models.leeway import Leeway
from opendrift.models.oceandrift3D import OceanDrift3D
from opendrift.models.openoil3D import OpenOil3D

current = ['norkyst_current_14Jul2019.nc', 'nws_current_14Jul2019.nc']
wind = ['ec_wind_14Jul2019.nc', 'arome_wind_14Jul2019.nc']


simulations = []
simnames = []

for w in wind:
    for c in current:
        simname = '%s_%s' % (c.split('_')[0], w.split('_')[0])
        simnames.append(simname)
        #o = Leeway()
        #o = OceanDrift3D()
        o = OpenOil3D()
        o.add_readers_from_list([c, w])
        o.fallback_values['land_binary_mask'] = 0
        o.seed_elements(lon=4.65, lat=60.5, number=1000, radius=1000,
                        wind_drift_factor=.02,
                        time=datetime(2019,7,14,0,30))
        o.run(duration=timedelta(hours=22))
        print(simname)
        o.plot(filename=simname+'.png')
        simulations.append(o)

simulations[0].plot(compare=simulations[1::],
                    legend=simnames)
