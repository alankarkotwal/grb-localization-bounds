import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
nside = 1024
skymap = np.zeros(hp.nside2npix(nside))
print hp.pixelfunc.nside2pixarea(nside, degrees=True)

ra0 = 180.0
dec0 = 30.0
ang0a = 44.0
ang0b = 46.0 
# Transient is within 44 and 46 degrees of ra,dec given above

ra1 = 130.0
dec1 = 45.0
ang1a = 60.0
ang1b = 65.0

def addring(nside, ra, dec, anga, angb):
    """
    take a map, and add 1 to elements between anga and angb from ra, dec
    """
    theta = np.deg2rad(90-dec)
    phi = np.deg2rad(ra)
    temp_map = np.zeros(hp.nside2npix(nside))
    assert angb > anga # else error
    # Everything from 0 to angb = 1
    pixlist = hp.query_disc(nside, hp.ang2vec(theta, phi), np.deg2rad(angb))
    temp_map[pixlist] += 1
    # now delete everything from 0 to anga
    pixlist = hp.query_disc(nside, hp.ang2vec(theta, phi), np.deg2rad(anga))
    temp_map[pixlist] -= 1
    return temp_map

skymap += addring(nside, ra0, dec0, ang0a, ang0b)
hp.mollview(skymap)
skymap += addring(nside, ra1, dec1, ang1a, ang1b)
hp.mollview(skymap)

n_selected_pixels = np.sum(skymap == 2) # make this 3 if you have three conditions, of course
sky_area = 360.**2 / np.pi * n_selected_pixels / len(skymap)
print sky_area
plt.show()

