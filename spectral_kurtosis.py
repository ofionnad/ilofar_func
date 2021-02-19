import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u

import pandas as pd
import xarray as xr

from scipy.signal import stft


######################
##### Functions ######
######################

def sb_to_f(sbs, obs_mode):
    nyq_dict = {3:1, 5:2, 7:3}
    nyq_zone = nyq_dict[obs_mode]
    clock_dict = {3:200, 4:160, 5:200, 6:160, 7:200} #MHz
    clock = clock_dict[obs_mode]
    nu = (nyq_zone-1. + sbs/512.) * (clock/2.)
    return nu * u.MHz

def plot_spectrum_plt(data, freqs):
    plt.imshow(data.T, aspect='auto', origin='lower',
           vmin=np.percentile(data.T, 5), 
           vmax=np.percentile(data.T, 95),
           extent=[0, data.shape[1], min(freqs[30:150]).value, max(freqs[30:150]).value])

def plot_stft_map(stft_i, i):
	"""
	i - particular index corresponding to individual subband stft
	stft_i - value returned from scipy.signal.stft()
			 should be a 3d tuple if passing more than 1 subband to stft()
	"""
	f, ax = plt.subplots(2)
	ax[0].pcolormesh(stft_i[1], stft_i[0], np.abs(stft_i[2][:,i,:]), shading='auto', 
					 vmin=np.percentile(np.abs(sft_i[2][:,i,:]), 2), 
					 vmax=np.percentile(np.abs(sft_i[2][:,i,:]), 98))

def spectral_kurtosis(stft_arr):
	"""
	Expects an stft array from scipy.signal.stft() with dimensions (freq, subband, time)
	"""
	frac_top = np.mean((np.abs(stft_arr)**4), axis=2)
	frac_bottom = np.mean((np.abs(atft_arr)**2), axis=2)**2
	return (frac_top / frac_bottom) - 2.

######################

######################
#### Reading data ####
######################

fname = 'udpoutput/jupiter-stokesI_0_2020-10-13T17:47:00_19563125244140'
sbs = np.arange(76, 320)
obs_mode = 3
trange = TimeRange("2020-10-13T17:47:00", 10.*u.min)
"""
fname = 'udpoutput/jupiter-stokesI_0_2020-10-13T17:47:00_19563125244140'
sbs = np.arange(76, 319)
obs_mode = 3
trange = TimeRange("2020-10-13T17:47:00", 10.*u.minute)
"""

xlabel = "Time"
ylabel = "Frequency (MHz)"
title = "Jupiter - Stokes I"

data = np.memmap(fname, np.float32, mode="c")
data = data.reshape(-1, sbs.shape[0])
data = np.flip(data, axis=1)

times = pd.timedelta_range(start="0 millisecond", periods=data.shape[0], freq='81.92U')
freqs = sb_to_f(sbs, obs_mode)

######################

#You can use xarray or just numpy, depending on the datatype and size/shape
dx = xr.DataArray(data, coords=[np.arange(len(times)),freqs], dims=["Time", "Frequency"])

# Normalise data (removes frequency dependence)
norms = dx.quantile(0.1, dim="Time")
dxn = dx/norms

## SK = <STFT>**4 / <STFT**2>**2 - 2.

# Get the Short Time Fourier Transform of each subband (just 3 for now)
dt = dxn[:, [30, 55, 130]]

#from scipy.signal
sft = stft(dt, window='hamming', axis=0, nperseg=564, nfft=2048)

sk = spectral_kurtosis(sft[2])