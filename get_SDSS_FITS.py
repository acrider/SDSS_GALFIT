def get_SDSS_FITS(run, rerun, camcol, field, ugriz, put_file):
    
    # This function downloads an SDSS image FITS file given four parameters and the color and places it in put_file.
    
    image = urllib.URLopener()

    get_file = 'http://das.sdss.org/www/cgi-bin/drC?RUN=' + run + '&RERUN=' + rerun + '&CAMCOL=' + camcol + '&FIELD=' + field + '&FILTER=' + ugriz

    if not(os.path.isfile(put_file)):
        print 'Retrieving ' + get_file
        image.retrieve(get_file, put_file)
    else:
        print 'File exists: ' + put_file
        
#--------------------------------------    

def retrieve_SDSS_params(objid):
    
    # This function returns seven parameters (run, rerun, camcol, field, obj, rowc, colc)
    # when given a DR7 objid number (e.g. from Allen et al. 2014).

    temp = 'http://cas.sdss.org/dr7/en/tools/explore/OETOC.asp?id='
    
    response = urllib2.urlopen(temp+objid)
    html     = response.read()
    start    = html.find('onLoad="loadSummary(')
    newId    = html[start+21:start+39] # This should look like 0x08290efc61020051
    newSpec  = html[start+49:start+67] # This should look like 0x0000000000000000
    newURL   = 'http://cas.sdss.org/dr7/en/tools/explore/summary.asp?id=' + newId + '&spec=' + newSpec
    response = urllib2.urlopen(newURL)
    html     = response.read()
    start    = html.find('colc')
    
    num_end = start + 1
    
    # Get the RUN
    class_start = html.find('class=',num_end)
    num_start   = html.find('>',class_start) + 1
    num_end     = html.find('<',num_start)
    run =  html[num_start:num_end]
    
    # Get the RERUN
    class_start = html.find('class=',num_end)
    num_start   = html.find('>',class_start) + 1
    num_end     = html.find('<',num_start)
    rerun =  html[num_start:num_end]
    
    # Get the CAMCOL
    class_start = html.find('class=',num_end)
    num_start   = html.find('>',class_start) + 1
    num_end     = html.find('<',num_start)
    camcol =  html[num_start:num_end]
    
    # Get the FIELD
    class_start = html.find('class=',num_end)
    num_start   = html.find('>',class_start) + 1
    num_end     = html.find('<',num_start)
    field =  html[num_start:num_end]
    
    # Get the OBJ
    class_start = html.find('class=',num_end)
    num_start   = html.find('>',class_start) + 1
    num_end     = html.find('<',num_start)
    obj =  html[num_start:num_end]
    
    # Get the ROWC
    class_start = html.find('class=',num_end)
    num_start   = html.find('>',class_start) + 1
    num_end     = html.find('<',num_start)
    rowc =  html[num_start:num_end]
    rowc = float(rowc)
  
    # Get the COLC
    class_start = html.find('class=',num_end)
    num_start   = html.find('>',class_start) + 1
    num_end     = html.find('<',num_start)
    colc =  html[num_start:num_end]
    colc = float(colc)  
        
    return run, rerun, camcol, field, obj, rowc, colc
    
#--------------------------------------    

def get_SDSS_img(put_file, rowc, colc):

    # These are things stored in the FITS file. r_img is the most important one.
    r_img    = pyfits.getdata(put_file)
    gain     = pyfits.getval(put_file,'GAIN') #GAIN = 4.59999990463257 / Gain averaged over all amplifiers (e/DN)  
    sky      = pyfits.getval(put_file,'SKY') # SKY = 132.366674648241 / sky level (DN/pixel) 
    dark_var = pyfits.getval(put_file,'DARK_VAR') # DARK_VAR = 1. / combined variance from dark cur. and read noise
    sky_err  = pyfits.getval(put_file,'SKYERR') # SKYERR  = 3.76900984380701E-12 / The error of average sky value in the frame    

    # Extract a 100x100 piece centered on the galaxy.
    boxsize = 100   
    xmin = int(rowc-boxsize/2)
    xmax = xmin+boxsize-1
    ymin = int(colc-boxsize/2)
    ymax = ymin+boxsize-1
    r_img = r_img[xmin:xmax,ymin:ymax]
    r_img = r_img - 1000 # This is the offset of 1000 built-in to SDSS images.
   
    return r_img, gain, sky, dark_var, sky_err
    
#--------------------------------------    

# IMPORT LIBRARIES 
  
import os.path
import urllib
import urllib2
import pyfits
import pylab as py
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.ndimage as ndimage

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
