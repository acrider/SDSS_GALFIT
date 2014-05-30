
#--------------------------------------    

# This is the Sersic function, where n=4 for de Vaucouleurs profile (giant ellipticals) 
# and n=1 for exponential profile (spiral, dwarf elliptical).
#
# https://en.wikipedia.org/wiki/Sersic_profile

def sersic(r, i0, k, n, background):
    return np.exp(np.log(i0) - k*r**(1.0/n)) + background

    
#-------------------------------------- 

# http://stackoverflow.com/questions/18304722/python-find-contour-lines-from-matplotlib-pyplot-contour

def get_contour_verts(cn):
    contours = []
    # for each contour line
    for cc in cn.collections:
        paths = []
        # for each separate section of the contour line
        for pp in cc.get_paths():
            xy = []
            # for each segment of that section
            for vv in pp.iter_segments():
                xy.append(vv[0])
            paths.append(np.vstack(xy))
        contours.append(paths)

    return contours

#-------------------------------------- 
# From http://stackoverflow.com/questions/7878398/how-to-extract-an-arbitrary-line-of-values-from-a-numpy-array

def get_slice(r_img, xmin, xmax, ymin, ymax):

    # extract values on line from r1, c1 to r2, c2
    num_points = 100
    xvalues = np.linspace(xmin, xmax, num_points)
    yvalues = np.linspace(ymin, ymax, num_points)
    
    zvalues = scipy.ndimage.map_coordinates(np.transpose(r_img), np.vstack((xvalues,yvalues)))
    
    # Set old version of x and y to the linear extrapolation...
    old_x = np.sqrt((xvalues-50)**2 + (yvalues-50)**2)
    old_y = zvalues
    
    return old_x, old_y
    
#-------------------------------------- 

# IMPORT LIBRARIES   

import sys
import urllib
import urllib2
import pyfits
import scipy.ndimage
import pylab as py
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import os.path
import csv
import scipy.ndimage as ndimage
from scipy.optimize import curve_fit
import math as mth
import subprocess
from get_SDSS_FITS import *

# Pick a few parameters for a sample galaxy.
objid = '587725578033955000'
FITS_directory = '/Users/acrider/Desktop/FITS/'
ugriz = 'r'
obj = ''

# Run retrieve_SDSS_params...
run, rerun, camcol, field, obj, rowc, colc = retrieve_SDSS_params(objid)

# ...so you can name the FIT file...    
put_file = FITS_directory + ugriz + '-' + run + '-' + camcol + '-' + field + '.fits'

# ...and then get it from the SDSS database.    
get_SDSS_FITS(run, rerun, camcol, field, ugriz, put_file)

# Now get the image data itself out of the FITS file.
r_img, gain, sky, dark_var, sky_err = get_SDSS_img(put_file, rowc, colc)

# Plot the 100x100 image with contours.
figure0 = py.figure(0)
figure0.clf()
f1 = py.imshow(r_img, cmap = cm.Greys_r)
N_contours = 20
cn = plt.contour(ndimage.gaussian_filter(r_img, sigma=1.0, order=0), N_contours) 
plt.show()
    
# Find the theta for contours that have 50,50 inside of them. If none exist, use 0.
thetas = []
nContours = len(cn.collections)
    
for jContours in xrange(0,nContours):
    nCenters = len(cn.collections[jContours].get_paths())
    for iPath in xrange(0, nCenters):
        p = cn.collections[jContours].get_paths()[iPath]
        if p.contains_point([50,50]):
            v = p.vertices
            contour_x = v[:,0]
            contour_y = v[:,1]
            contour_r = np.sqrt((contour_x - 50)**2 + (contour_y - 50)**2) 
            contour_theta = np.arctan2(contour_y - 50, contour_x - 50)
            thetas.append(contour_theta[contour_r.argmax()]) #Tack on this best theta if this path contains 50,50.
if len(thetas) > 0:
    theta = np.median(thetas)
else:
    theta = 0.0

# Plot a line for the slice you want to take...
x0=50
y0=50
radius = 40
x1 = x0 + int(radius * np.cos(theta))
y1 = y0 + int(radius * np.sin(theta))
plt.figure(0)
plt.plot([x0,x1],[y0,y1], color='w')
plt.show()

# Throw out the old x, y and grab an arbitrary slice instead!
x, y = get_slice(r_img, x0, x1, y0, y1)
    
figure1 = py.figure(1)
figure1.clf()
f2 = plt.plot(x,y)

plt.show()

#---------------------------------------------------

p = [max(y)-min(y), 1.0, 2.0, min(y)]
sigma = np.sqrt((y + sky)/gain) + 1 * (dark_var + sky_err)
#print 'MIN(y) = ', min(y)
#print 'MAX(y) = ', max(y)

xfine = np.arange(x0-50,x1-50,0.1)
yfine = sersic(xfine, p[0], p[1], p[2], p[3])
figure2 = py.figure(2)
figure2.clf()
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.errorbar(x,y,yerr=sigma, fmt='o')
plt.plot(xfine,yfine,'k--')
plt.show()
    
try:
    popt, pcov = curve_fit(sersic, x, y, p0=p, sigma=sigma)
except RuntimeError:
    popt = p
    fit_fail = 1
else:
    fit_fail = 0
        
yfine = sersic(xfine, popt[0], popt[1], popt[2], popt[3])
plt.plot(xfine, yfine,'r-')

print 'GAIN = ', gain
print 'I0 = ', popt[0]
print 'k  = ', popt[1]
print 'n  = ', popt[2], '(n=1 for spirals, n=4 for ellipticals)'
half_light_radius = (-np.log(0.5) / popt[1])**popt[1]
print

figure2 = py.figure(2)
figure2.clf()
  
plt.errorbar(x,y,yerr=sigma, fmt='o')
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
if fit_fail == 0:
    plt.plot(xfine, yfine,'r-')
plt.show()
        


