#-------------------------------------- 
def run_galfit():

    p = subprocess.Popen(['/Users/acrider/galfit/galfit galfit.feedme'], cwd='/Users/acrider/Desktop/my-example',
        shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in p.stdout.readlines():
        print line,
    retval = p.wait()
    
    return 

#-------------------------------------- 

def create_galfit_feedme(feedme_file, **kwargs):
    
        # These are all the default settings for the keyword arguments (kwargs) passed to the function.
        settings = {    'input_image'    :'galfits.fits',
                        'output_image'   :'imgblock.fits',
                        'xmin'           : '0',
                        'xmax'           : '99',
                        'ymin'           : '0',
                        'ymax'           : '99',
                        'x'              : '50.0',
                        'y'              : '50.0'}

        # The default settings can be overridden with kwargs.
        for item in settings:     
            if item in kwargs:
                settings[item] = str(kwargs[item])
                print 'Setting '  +  item + ' = ' + settings[item]      
               
        # Open the file and write it...
        f = open(feedme_file, 'w\n')
        
	f.write('\n')
	f.write('===============================================================================\n')
	f.write('# IMAGE and GALFIT CONTROL PARAMETERS\n')
   	f.write('A) '    + settings['input_image']     + ' # Input data image (FITS file)\n')
	f.write('B) '    + settings['output_image']    + ' # Output data image block\n')
	f.write('C) none                # Sigma image name (made from data if blank or "none") \n')
	f.write('D) none 	        # Input PSF image and (optional) diffusion kernel\n')
	f.write('E) 1                   # PSF fine sampling factor relative to data \n')
	f.write('F) none                # Bad pixel mask (FITS image or ASCII coord list)\n')
	f.write('G) none                # File with parameter constraints (ASCII file) \n')
	f.write('H) ' + settings['xmin'] + ' ' + settings['xmax'] + ' ' + settings['ymin'] + ' ' + settings['ymax'] + ' # Image region to fit (xmin xmax ymin ymax)\n')
	f.write('I) 100    100          # Size of the convolution box (x y)\n')
	f.write('J) 26.563              # Magnitude photometric zeropoint \n')
	f.write('K) 0.396127 0.396127   # Plate scale (dx dy)    [arcsec per pixel]\n')
	f.write('O) regular             # Display type (regular, curses, both)\n')
	f.write('P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n')
	f.write('\n')
	f.write('# INITIAL FITTING PARAMETERS\n')
	f.write('#\n')
	f.write('#   For object type, the allowed functions are: \n')
	f.write('#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat, \n')
	f.write('#       ferrer, powsersic, sky, and isophote. \n')
	f.write('#  \n')
	f.write('#   Hidden parameters will only appear when they''re specified:\n')
	f.write('#       C0 (diskyness/boxyness), \n')
	f.write('#       Fn (n=integer, Azimuthal Fourier Modes),\n')
	f.write('#       R0-R10 (PA rotation, for creating spiral structures).\n')
	f.write('# \n')
	f.write('# -----------------------------------------------------------------------------\n')
	f.write('#   par)    par value(s)    fit toggle(s)    # parameter description \n')
	f.write('# -----------------------------------------------------------------------------\n')
	f.write('\n')
	f.write('# Object number: 1\n')
	f.write(' 0) sersic                 #  object type\n')
	f.write(' 1) ' + settings['x'] + ' ' + settings['y'] + ' 1 1  	   #  position x, y\n')
	f.write(' 3) 20.0890     1          #  Integrated magnitude	\n')
	f.write(' 4) 20.000      1          #  R_e (half-light radius)   [pix]\n')
	f.write(' 5) 4.0000		1          #  Sersic index n (de Vaucouleurs n=4) \n')
	f.write(' 6) 0.0000      0          #     ----- \n')
	f.write(' 7) 0.0000      0          #     ----- \n')
	f.write(' 8) 0.0000      0          #     ----- \n')
	f.write(' 9) 1.0      	1          #  axis ratio (b/a)  \n')
	f.write('10) 0    		1          #  position angle (PA) [deg: Up=0, Left=90]\n')
	f.write(' Z) 0                      #  output option (0 = resid., 1 = Dont subtract) \n')
	f.write('\n')
	f.write('# Object number: 2\n')
	f.write(' 0) sky                    #  object type\n')
	f.write(' 1) 1147.64246651 1        #  sky background at center of fitting region [ADUs]\n')
	f.write(' 2) 0.0000      0          #  dsky/dx (sky gradient in x)\n')
	f.write(' 3) 0.0000      0          #  dsky/dy (sky gradient in y)\n')
	f.write(' Z) 0                      #  output option (0 = resid., 1 = Dont subtract) \n')
	f.write('\n')
	f.write('================================================================================\n')

        f.close()

#--------------------------------------    
def get_GALFIT_results(filename):

    # These are things stored in the FITS file. r_img is the most important one.
    r_img    = pyfits.getdata(filename)
    #gain     = pyfits.getval(put_file,'GAIN') #GAIN = 4.59999990463257 / Gain averaged over all amplifiers (e/DN)  
    #sky      = pyfits.getval(put_file,'SKY') # SKY = 132.366674648241 / sky level (DN/pixel) 
    #dark_var = pyfits.getval(put_file,'DARK_VAR') # DARK_VAR = 1. / combined variance from dark cur. and read noise
    #sky_err  = pyfits.getval(put_file,'SKYERR') # SKYERR  = 3.76900984380701E-12 / The error of average sky value in the frame    
   
    return r_img
    
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

import os.path
import urllib
import urllib2
import pyfits
import pylab as py
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.ndimage as ndimage

from get_SDSS_FITS import *

# Pick a few parameters for a sample galaxy.

# The user must set this to the directory and file that needs to be read in.
#inputfile = '/Users/acrider/Desktop/AGN 2014/BPT_SF_data.csv'
inputfile = '/Users/acrider/Desktop/AGN 2014/BPT_AGN_data.csv'

# Read in the CSV file.
data = [];
with open(inputfile, 'rb') as f:
    csvReader = csv.reader(f, delimiter=',',skipinitialspace=True)
    headers = csvReader.next()
    for row in csvReader:
        data.append(row);
    data = np.asarray(data)

# Get a list of all of the unique DR7 galaxy IDs.
oldObjID       = data[:,1]
oldObjID       = np.unique(oldObjID)
nUnique        = len(oldObjID)

# Create an empty array to store the Sersic fit parameters.
sersic_n = np.zeros(nUnique)
sersic_n = sersic_n - 1.0

# Set a few other things...
FITS_directory = '/Users/acrider/Desktop/FITS/'
ugriz = 'r'
obj = ''

# MAIN LOOP
for i in xrange(0,nUnique):
#for i in xrange(0,1):

    objid = oldObjID[i]
    
    # Run retrieve_SDSS_params...
    run, rerun, camcol, field, obj, rowc, colc = retrieve_SDSS_params(objid)

    # ...so you can name the FIT file...    
    put_file = FITS_directory + ugriz + '-' + run + '-' + camcol + '-' + field + '.fits'

    # ...and then get it from the SDSS database.    
    get_SDSS_FITS(run, rerun, camcol, field, ugriz, put_file)

    # Now get the image data itself out of the FITS file.
    r_img, gain, sky, dark_var, sky_err = get_SDSS_img(put_file, rowc, colc)

    # Plot the 100x100 image with contours.
    figure_SDSS = py.figure(0)
    figure_SDSS.clf()
    py.imshow(r_img, cmap = cm.Greys_r)
    N_contours = 20
    plt.contour(ndimage.gaussian_filter(r_img, sigma=1.0, order=0), N_contours) 
    plt.show()

    x = colc # Which way should these be?!?
    y = rowc #

    xmin = int(x - 50)
    xmax = xmin + 99
    ymin = int(y - 50)
    ymax = ymin + 99
    
    output_image = '/Users/acrider/Desktop/my-example/imgblock.fits'
    
    create_galfit_feedme('/Users/acrider/Desktop/my-example/galfit.feedme',
                            input_image=put_file,
                            output_image=output_image,
                            xmin=xmin, xmax=xmax,
                            ymin=ymin, ymax=ymax,
                            x=x,y=y )
    run_galfit()

    img_original  = pyfits.getdata(output_image,1)
    img_model     = pyfits.getdata(output_image,2)
    img_residuals = pyfits.getdata(output_image,3)

    tmp_n      = pyfits.getval(output_image,'1_N', 2)
    tmp_n      = tmp_n.replace('*','') 
    sersic_n[i] = float(tmp_n.split(' +/- ')[0])
    sersic_n_error = float(tmp_n.split(' +/- ')[1])
    print 'Sersic n = ', sersic_n[i], '+/-', sersic_n_error
     
    # Two subplots, the axes array is 1-d
    plt.figure(1)
    plt.clf()
    plt.subplot(1,3,1)
    plt.imshow(img_original)
    plt.subplot(1,3,2)
    plt.imshow(img_model)
    plt.subplot(1,3,3)
    plt.imshow(img_residuals)
    plt.show()
    
    # Create histogram of Sersic n values.    
    plt.figure(2)
    plt.clf()

    bins = np.arange(51)/10.0
    plt.hist(sersic_n, bins=bins)
    plt.show()

