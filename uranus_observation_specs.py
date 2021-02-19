import astropy.coordinates as coords 
import astropy.time as Time 

t = Time.Time('2020-12-15T20:04:00', format='isot', scale='utc')

c = coords.get_body('uranus', t)

print("RA = {} rad".format(c.ra.rad))
print("Dec = {} rad".format(c.dec.rad))
print("Observation Time = {} {}".format(t.value, t.scale))