import numpy as np
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric
from astropy.constants import c, au
from numpy import pi
from matplotlib import pyplot as plt
import healpy as hp
import pdb

delT = 0.1
theta = 30 * pi / 180

nside = 1024
skymap = np.zeros(hp.nside2npix(nside))
print hp.pixelfunc.nside2pixarea(nside, degrees=True)
hp.mollview(skymap)

dateInt = Time(['2023-01-01T00:00:00.000', '2028-01-01T00:00:00.000'], format='isot', scale='utc')
numPoints = 2000 

jdInt = dateInt.jd
jds = Time(np.linspace(jdInt[0], jdInt[1], numPoints), format='jd', scale='utc');

sunCoords = get_body_barycentric('sun', jds)
earthCoords = get_body_barycentric('earth', jds) - sunCoords
venusCoords = get_body_barycentric('venus', jds) - sunCoords
marsCoords = get_body_barycentric('mars', jds) - sunCoords

veDists = np.zeros((numPoints, 1))
vmDists = np.zeros((numPoints, 1))
emDists = np.zeros((numPoints, 1))
locErrs = np.zeros((numPoints, 1))
locErrsNoVen = np.zeros((numPoints, 1))

ve = (venusCoords - earthCoords)
vm = (venusCoords - marsCoords)
em = (earthCoords - marsCoords)

for i in range(numPoints):
	veDists[i] = au * np.sqrt(ve[i].x ** 2 + ve[i].y ** 2 + ve[i].z ** 2)
	vmDists[i] = au * np.sqrt(vm[i].x ** 2 + vm[i].y ** 2 + vm[i].z ** 2)
	emDists[i] = au * np.sqrt(em[i].x ** 2 + em[i].y ** 2 + em[i].z ** 2)
	locErrs[i] = 3.3 * delT * np.sqrt(2) * c / max(veDists[i], vmDists[i], emDists[i]); 
	locErrsNoVen[i] = 3.3 * delT * np.sqrt(2) * c / emDists[i]; 

plotDates = jds.plot_date - jds.plot_date[0]
plt.plot(plotDates, 206265 * c * delT / (np.sin(theta) * veDists))
plt.plot(plotDates, 206265 * c * delT / (np.sin(theta) * vmDists))
plt.plot(plotDates, 206265 * c * delT / (np.sin(theta) * emDists))
plt.legend(['Venus-Earth', 'Venus-Mars', 'Earth-Mars'])
plt.xlabel('Days since 1 Jan 2023')
plt.ylabel('Pairwise arcsec resolution')
plt.savefig('resolutions.png')

plt.figure()
plt.plot(plotDates, locErrs)
plt.plot(plotDates, locErrsNoVen)
plt.legend(['With Venus Probe', 'Without Venus Probe'])
plt.xlabel('Days since 1 Jan 2023')
plt.ylabel('Upper Bound on Localization Fraction')
plt.savefig('areas.png')
