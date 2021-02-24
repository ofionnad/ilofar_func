import numpy as np 
import matplotlib.pyplot as plt
import astropy.units as u 
import astropy.coordinates as coords
from astropy.time import Time
from astropy.timeseries import TimeSeries
from pygdsm import GSMObserver2016, GlobalSkyModel2016
from datetime import datetime
from healpy import ang2vec, ang2pix, get_nside
import lofarantpos.db

from observation import Observation

def sb_to_f(sbs, obs_mode):
    """
    Convert subband number to a frequency in MHz
    """
    nyq_dict = {3:1, 5:2, 7:3}
    nyq_zone = nyq_dict[obs_mode]
    clock_dict = {3:200, 4:160, 5:200, 6:160, 7:200} #MHz
    clock = clock_dict[obs_mode]
    nu = (nyq_zone-1. + sbs/512.) * (clock/2.)
    return nu * u.MHz

def astropytime_to_datetime(t):
    ds = t.value.split('-')
    hs = ds[2].split('T')
    ts = hs[1].split(':')
    ss = ts[2].split('.')
    return datetime(int(ds[0]), int(ds[1]), int(hs[0]), int(ts[0]), int(ts[1]), int(ss[0]), int(ss[1]))

def kondratiev_temp(freqs):
    """
    Get the LBA antenna temperatures from the sixth order polynomial Kondratiev uses from Winholds+2011
    """
    T_inst_poly = [6.2699888333e-05, -0.019932340239, 2.60625093843, -179.560314268, 6890.14953844, -140196.209123, 1189842.07708]
    temps = np.poly1d(T_inst_poly)
    return temps(freqs.value)

def get_min_d(positions):
    ds = []
    for j in positions:
        distances = j - positions
        d = np.sqrt(distances[:,0]**2 + distances[:,1]**2 + distances[:,2]**2)
        dmin = np.sort(d)[1]
        ds.append(dmin)
    return ds


def get_lofar_aeff_max(freqs, n_elem, station_str='IE613LBA'):
    """
    Calculate the max Aeff 
    """

    c = 300.0 #speed of light when converting from MHz
    l = c / freqs #wavelength in metres

    db = lofarantpos.db.LofarAntennaDatabase()
    positions = db.antenna_pqr(station_str)

    ds = get_min_d(positions)
    print(ds)

    aeff = []
    for j in l.value:
        aef = []
        for s in ds:
            #aef.append((j**2.)/3.)
            aef.append(np.minimum((j**2.)/3., np.pi*(s**2.)/4.))
        aeff.append(n_elem*np.mean(aef))

    return np.array(aeff)

def galactic_noise(freqs, gal_coords):
    """
    Get the background sky noise from PyGDSM from D. Price on github (https://github.com/telegraphic/pygdsm)
    """

    gsm = GlobalSkyModel2016(freq_unit='MHz', data_unit='TCMB')
    #freqs = np.linspace(15,30,79)
    map_cube = gsm.generate(freqs.value)
    vec = ang2vec(np.radians(gal_coords.l.value), np.radians(gal_coords.b.value))
    ipix = ang2pix(nside=get_nside(map_cube[0]), theta=np.radians(gal_coords.l.value), phi=np.radians(gal_coords.b.value))

    T_sky = map_cube[:, ipix]

    return T_sky

def sensitivity(freqs, gal_coords, n_elem=96):
    """
    From kondratiev+2016
    """
    kb = 1.38e-16 #in cgs (so output will be erg/m^2)
    beta = 1 #we are using 16 bit observations so digitisation correction is not needed
    n_pol = 2
    #T_obs = 81.92e-6 #observation resolution
    T_obs = 1e-3 #resample data to 1 ms..
    #T_obs = 1 #1s resolution
    bw = (freqs[1].value-freqs[0].value)*1e6 #bandwidth
    print(bw)
    T_sky = galactic_noise(freqs, gal_coords)
    T_ant = kondratiev_temp(freqs)
    print(T_sky, T_ant)
    T_sys = T_ant + T_sky
    Aeff = get_lofar_aeff_max(freqs, n_elem)
    Aeff = Aeff*1e4 #convert to cm^2 from m^2
    print(Aeff)
    gain = Aeff / (2.*kb)

    return beta*T_sys / (gain*np.sqrt(n_elem*(n_elem-1)*n_pol*T_obs*bw))

"""
gsm = GlobalSkyModel2016(freq_unit='MHz', data_unit='TCMB')
#freqs = np.linspace(15,30,79)
map_cube = gsm.generate(freqs.value)

vec = ang2vec(np.radians(gal_coords.l.value), np.radians(gal_coords.b.value))
ipix = ang2pix(nside=get_nside(map_cube[0]), theta=np.radians(gal_coords.l.value), phi=np.radians(gal_coords.b.value))

print(map_cube.shape)
plt.loglog(map_cube[:,ipix])  # Random pixel
#plt.loglog(map_cube[:,ipix]) # Another random pixel
plt.xlabel("Frequency [MHz]")
plt.ylabel("Temperature [K]")
plt.savefig('uranus_skymodel/mollweide_multi_freq.png')
plt.close()

model_freq = 20 #MHz
ov = GSMObserver2016()
ov.lat, ov.long, ov.elev = ilofar.lat.value, ilofar.lon.value, ilofar.height.value
#ov.date = astropytime_to_datetime(t)
#ov.generate(model_freq)


samples = 20
tdelta = duration/samples
tsl = TimeSeries(time_start=start, time_delta=tdelta, n_samples=samples)
"""
if False:
    for i,j in enumerate(tsl):
        print(j[0])
        ov.date = astropytime_to_datetime(j[0])
        print(ov.date)
        ov.generate(model_freq)
        ov.view(logged=True)
        #plt.show()
        plt.savefig('uranus_skymodel/ortho_{}'.format(i))
        plt.close()

if __name__=="__main__":

    #sbs = np.arange(76,198)[:78]
    sbs = np.arange(76,198)
    freqs = sb_to_f(sbs, 3)
    start = '2020-12-15T20:04:00' #isot format
    duration = 176*60*u.second # seconds
    n_elem = 96

    myobs = Observation('uranus', start, 'isot', duration.value, freqs)

    t = Time('2020-12-15T20:04:00', format='isot')
    c = coords.get_body('uranus', t)
    ilofar = coords.EarthLocation(lat='53.095', lon='-7.9218', height=100*u.m)

    aa = coords.AltAz(location=ilofar, obstime=t)

    altaz_coords = c.transform_to(aa) #at ilofar
    gal_coords = c.transform_to(coords.Galactic())
    print(myobs)
    print('AltAz Coordinates\n'+'alt: {}\taz: {}'.format(altaz_coords.alt, altaz_coords.az))
    print('Galactic Coordinates\n'+'l: {}\tb: {}\n'.format(gal_coords.l, gal_coords.b))

    s = sensitivity(freqs, gal_coords, n_elem=n_elem)
    print(s)
    jy = 1e23*s