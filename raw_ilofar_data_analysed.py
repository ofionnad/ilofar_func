import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from sunpy.time import TimeRange
import matplotlib.dates as mdates
import datetime as dt

import xarray as xr

import datashader as ds
import datashader.transfer_functions as tf
from colorcet import fire, gray, kr

import holoviews as hv

from holoviews import opts
hv.extension('bokeh', 'matplotlib')


def sb_to_f(sbs, obs_mode):
    nyq_dict = {3:1, 5:2, 7:3}
    nyq_zone = nyq_dict[obs_mode]
    clock_dict = {3:200, 4:160, 5:200, 6:160, 7:200} #MHz
    clock = clock_dict[obs_mode]
    nu = (nyq_zone-1. + sbs/512.) * (clock/2.)
    return nu * u.MHz

"""
fname = 'udpoutput/cygA-1s-stokesI-notimedec_0_2020-10-13T17:47:00_ld'
sbs = np.arange(76, 197)
obs_mode = 3
trange = TimeRange("2020-10-13T17:47:00", 1.*u.second)
"""
fname = 'udpoutput/jupiter-stokesI_0_2020-10-13T17:47:00_19563125244140'
sbs = np.arange(76, 320)
obs_mode = 3
trange = TimeRange("2020-10-13T17:47:00", 10.*u.minute)

xlabel = "Time"
ylabel = "Frequency (MHz)"
title = "Jupiter - Stokes I"


#data = np.fromfile(fname) #too slow for longer data
data = np.memmap(fname, np.float32, mode="c")
data = data.reshape(-1, sbs.shape[0])
data = np.flip(data, axis=1)

freqs = sb_to_f(sbs, obs_mode)

#times = pd.date_range("2020-10-13", periods=df.shape[0], freq='81.92U')
times = pd.timedelta_range(start="0 millisecond", periods=data.shape[0], freq='81.92U')

#dx = xr.DataArray(data, coords=[times, freqs], dims=["Time", "Frequency"])
dx = xr.DataArray(data, coords=[np.arange(len(times)),freqs], dims=["Time", "Frequency"])


### Normalise frequency response by 10th quantile across all times
normalised = dx.quantile(0.1, dim='Time')

dx_norm = dx / normalised


# cut out the highest frequencies, no jupiter data there
dx_short = dx_norm[:,27:130]

