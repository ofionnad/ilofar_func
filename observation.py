from typing import ClassVar
import astropy.coordinates as coords 
import astropy.time as Time
import astropy.units as u
import os

class Observation:

    def __init__(self, object_name, time, time_format, duration, frequency) -> None:
        self.name = object_name 
        self.time = Time.Time(time, format=time_format)

        if self.name in coords.solar_system_ephemeris.bodies:
            self.coords = coords.get_body(self.name, self.time)
        elif self.name == 'moon':
            self.coords = coords.get_moon(self.time)
        elif self.name == 'sun':
            self.coords = coords.get_sun(self.time)
        else:
            try:
                self.coords = coords.get_icrs_coordinates(self.name)
            except coords.name_resolve.NameResolveError:
                print('I could not find that object...')
                print('Please manually set coordinates')

        self.duration = duration * u.second
        self.frequency = frequency * u.MHz

    def __str__(self):
        return 'Target: {}\n'.format(self.name) + \
        'Time: {}\n'.format(self.time) + \
        'Coordinates (radians):\n\n' + \
        'RA: {}\n'.format(self.coords.ra.rad) + \
        'Dec: {}\n\n'.format(self.coords.dec.rad) + \
        'Duration: {}\n'.format(self.duration) + \
        'Frequency: {}\n'.format(self.frequency)

    def run_dreambeam(self, command, telescope, configuration, station, beam_model, step, frequency):
        os.system("mkdir {}".format(self.name+self.time.value))
        bashcmd = "pointing_jones {} {} {} {} {} {} {} {} {} {} {} > {}".format(command, 
        telescope, configuration, station, beam_model, self.time.value[:-4], self.duration.value,
        step, self.coords.ra.rad, self.coords.dec.rad, frequency, self.name+self.time.value+'/'+str(frequency/1.e6)+'MHz_jones_matrices.dat')
        print(bashcmd)
        os.system(bashcmd)