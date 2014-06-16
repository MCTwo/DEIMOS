'''
The purpose of this program is to match the spectroscopic objects of a deimos slit
mask with objects from an image catalog.
'''
from __future__ import division
import numpy
import pyfits
import tools
import ds9 # ds9 and pyds9 should be installed from http://ds9.si.edu/site/Home.html

###########################
### USER INPUTS
###########################
path = '/sandbox/deimos/ciza1A/2013jul14/'
maskname = 'ciza1A'
zspecfile = 'zspec.dawson.ciza1A.final.fits'
tolerance = 2 #matching tolerance (arcsec) any objects within this separation will be considered a match
# image catalog in pandas csv format with first row as the header row
imgcat = '/Users/dawson/OneDrive/Research/Clusters/CIZAJ2242/catalog/SubaruCFHT_Sextractor/I_pandas.csv' #path/name of the image catalog
objkey = 'NUMBER' #ttype name of the unique object id column
imgcoord = ('X_WORLD','Y_WORLD') #ttype name of the ra and dec columns in the image catalog
mag = 'MAG_AUTO'
outputfile = '/sandbox/deimos/ciza1A/matchcat_ciza1A_subaru.txt'

### Fits image file input and mask region file

# Subaru image and visulation parameters
file_fits = '/Users/dawson/OneDrive/Research/Clusters/CIZAJ2242/fits/optical/Subaru/I.fits'

# Define image scale
scale = 'log' #e.g. 'linear', 'log', 'power', etc. 
# Define (min, max) pixel scale values
scale_limits = (0,200)

# Define Subaru and HST object region files
region = '/Users/dawson/OneDrive/Observing/Keck2013a/CIZAJ2242/ciza_m1A.reg'

# this scale is divided by the slitlength (in arcsec) to determine the zoom
# for that object
zoom_scale = 30

###########################
### PROGRAM
###########################
## Setup the ds9 image
# setup ds9
# call ds9
d = ds9.ds9()
# turn off colorbar
d.set('colorbar no')

# load the subaru image in the first frame
cmd = 'file '+file_fits
d.set(cmd)

## correct the WCS issues (deleted PV1_ and PV2_ values)
#d.set('wcs replace Subaru.wcs')

# change the scale to log and apply set limits
cmd = 'scale '+scale
d.set(cmd)
cmd = 'scale limits {0} {1}'.format(scale_limits[0],scale_limits[1])
d.set(cmd)

# load Subaru object regions into both frames
cmd = 'regions load all '+region
d.set(cmd)

# Gather the basic slit info tables from the bintabs.fits file
binfile = maskname+'.bintabs.fits'
hdubin = pyfits.open(path+binfile)
# Target table information, similar to dsim .lst information
tb_targets = hdubin[1].data
# Table of slitmask information
tb_mask = hdubin[2].data
# Table of each slit's physical properties
tb_slits = hdubin[3].data
# Table that maps the id's of the previous two tables
tb_map = hdubin[4].data

## Work with the zspec results files

# Read in the zspec file contents
hduzspec = pyfits.open(path+'../../'+zspecfile)
tb_zspec = hduzspec[1].data

# Create an array with all the slit numbers
slitnumbers = tb_slits.field('SLITNAME')

# Load the image catalog
cat_img = tools.readcatalog(imgcat)
key_img = tools.readheader(imgcat)

#Create the ouput file and write header information
fh = open(outputfile,'w')
fh.write('#This catalog was created by slitcatmatch.py and matches deimos spectrographic\n')
fh.write('#traces with a catalog of images.\n')
fh.write('#ttype0 = target_objid\n')
fh.write('#ttype1 = z\n')
fh.write('#ttype2 = zerr\n')
fh.write('#ttype3 = quality\n')
fh.write('#ttype4 = mask\n')
fh.write('#ttype5 = slit\n')
fh.write('#ttype6 = which_trace\n')
fh.write('#ttype7 = y_trace\n')
fh.write('#ttype8 = ra_trace\n')
fh.write('#ttype9 = dec_trace\n')
fh.write('#ttype10 = objid\n')
fh.write('#ttype11 = ra_obj\n')
fh.write('#ttype12 = dec_obj\n')
fh.write('#ttype13 = mag_obj\n')
fh.write('#ttype14 = matchdelta\n')
fh.write('#ttype15 = comment\n')
fh.close()
def match(slit_i,which_trace,cat,key,coord,objkey,mag,tolerance,outputfile):
    #Filter the tables keeping only the current slit
    tb_s = tb_slits[tb_slits.field('SLITNAME')==slit_i]
    # slitid
    dslitid = tb_s.field('DSLITID')
    tb_m = tb_map[tb_map.field('DSLITID')==dslitid]
    # objectid
    objectid = tb_m.field('OBJECTID')
    tb_t = tb_targets[tb_targets.field('OBJECTID')==objectid]
    # object (e.g.: DLS photometric objid)
    obj = tb_t.field('OBJECT')[0]
    # zspec info
    # converting '00#' to '#' notation 
    tb_zs = tb_zspec[tb_zspec.field('SLITNAME')==str(int(slit_i))]

    # Define the slit astrometric/geometric properties
    slitra = tb_s.field('SLITRA')[0]
    slitdec = tb_s.field('SLITDEC')[0]
    slitlen = tb_s.field('SLITLEN')[0]
    slitwid = tb_s.field('SLITWID')[0]
    slitpa = tb_s.field('SLITLPA')[0] #pa for long axis of the slit +ccw from north
    z = tb_zs.field('Z')[0]
    zerr = tb_zs.field('Z_ERR')[0]
    quality = tb_zs.field('ZQUALITY')[0]
    slitcomment = tb_zs.field('COMMENT')[0]
    
    # Define the pa of the mask
    maskpa = tb_mask.field('PA_PNT')[0]
    # Determine the minimum coterminal angle of the slit respect to the mask
    phi = tools.coterminal(slitpa-maskpa)
    
    # read in the primary trace information
    if which_trace == 'primary':
        hdutrace = pyfits.open(path+'spec1d.{0}.{1}.{2}.fits'.format(maskname,slit_i,obj))
    else:
        hdutrace = pyfits.open(path+'spec1d.{0}.{1}.{2}.fits'.format(maskname,slit_i,which_trace))
    tb_trace = hdutrace[1].data #Note using blue side info, redside typically redundent
    # Get the y-position of the trace in pixels from the bottom of the slit
    # and convert to arcsec using pixel scale of 0.1185 arcsec/pix
    pixscale = 0.1185 #arcsec/pix
    y = (tb_trace.field('OBJPOS')[0]+1)*pixscale
    # Get the FWHM of the trace
    y_fwhm = tb_trace.field('FWHM')[0]*pixscale
    # angular distance of the trace from the center of the slit
    dc = y - slitlen/2. #positive toward top of slit and negative towards bottom
    # Calculate the ra,dec of the trace object based on the pa and slit ra,dec
    if phi > 90 and phi < 270:
        #then essentailly zspec treats the slit as upside down and y becomes with repsect to the top of the slit thus dc should have the opposite sign
        dc = -dc
        ra_trace,dec_trace = tools.angendpt(slitra,slitdec,dc/60.,slitpa)
    else:
        #then zspec treats the slit as righ side up and y is measured with respect to the bottom of the slit
        ra_trace,dec_trace = tools.angendpt(slitra,slitdec,dc/60.,slitpa)
    
    # Try to find a match within the image catalog
    
    # Do a quick trim of the catalog to reduce calculation time
    ramin = slitra-slitlen/(60.**2*2.*numpy.cos(slitdec*numpy.pi/180.))
    ramax = slitra+slitlen/(60.**2*2.*numpy.cos(slitdec*numpy.pi/180.))
    decmin = slitdec-slitlen/(60.**2*2.)
    decmax = slitdec+slitlen/(60.**2*2.)
    mask_cat_ra = numpy.logical_and(cat[coord[0]]>ramin,cat[coord[0]]<ramax)
    mask_cat_dec = numpy.logical_and(cat[coord[1]]>decmin,cat[coord[1]]<decmax)
    mask = numpy.logical_and(mask_cat_ra,mask_cat_dec)
    cat_flt = cat[mask]
    N = numpy.shape(cat_flt)[0]
    
    # pan to the current slit
    cmd = 'pan to {0} {1} wcs fk5 degrees'.format(slitra,slitdec)
    d.set(cmd)
    # place the crosshair over the object of interest
    cmd = 'crosshair {0} {1} wcs fk5 degrees'.format(ra_trace,dec_trace)
    d.set(cmd)
    # zoom in on the object
    zoom = zoom_scale/slitlen
    cmd = 'zoom to {0}'.format(zoom)
    d.set(cmd)
    

    # Calculate the angular separation between all objects and the trace object
    j=0
    delta = numpy.zeros(N)
    for i in range(N):
        ra = cat_flt[i,key[coord[0]]]
        dec = cat_flt[i,key[coord[1]]]
        delta[i] = numpy.abs(tools.angdist(ra,dec,ra_trace,dec_trace)*60**2)
        if delta[i] < tolerance:
            j+=1
    if j==1:
        #there was a single match satisfying the tolerence 
        cat_flt = cat_flt[delta<tolerance,:]
        delta = delta[delta<tolerance]
        match_id = cat_flt[0,key[objkey]]
        match_ra = cat_flt[0,key[coord[0]]]
        match_dec = cat_flt[0,key[coord[1]]]
        match_delta = delta[0]
        match_mag = cat_flt[0,key[mag]]
    elif j == 0:
        if numpy.size(delta) != 0:
            # sort match_delta smallest to largest
            index = numpy.argsort(delta)
            delta=delta[index]
            cat_flt = cat_flt[index,:]
            print 'slitcatmatch: No catalog matches were found for this trace.'
            print 'Slit {0} {1}'.format(slit_i,which_trace)
            print 'The closest objects to the trace are:'
            print 'Object\tRA\t\tdec\tSeparation (arcsec)\tMagnitude'
            for k in range(numpy.size(delta)):
                print '{0}\t{1:0.5f}\t{2:0.4f}\t{3:0.3f}\t{4:0.1f}'.format(k,cat_flt[k,key[coord[0]]],cat_flt[k,key[coord[1]]],delta[k],cat_flt[k,key[mag]])
                #display a region at the object's location, with label
                cmd = '"fk5; circle point {0:0.6f} {1:0.5f}'.format(cat_flt[k,key[coord[0]]],cat_flt[k,key[coord[1]]])+' # color=red text{'+'{0}'.format(k)+'}"'
                d.set('regions', cmd)
                
            print '{0}\tSelect none.'.format(numpy.size(delta))
            selection = raw_input('Enter the number of the correct object match: ')
            if numpy.size(numpy.arange(k+1)==int(selection))==0:
                selection = rawinput("Input invalid. Please enter a valid number.: ")
            if selection == str(numpy.size(delta)):
                # Don't associate the trace with an object
                match_id = match_ra = match_dec = match_delta = match_mag = -99
            elif numpy.size(numpy.arange(k+1)==int(selection))!=0:
                selection=int(selection)
                match_id = cat_flt[selection,key_img[objkey]]
                match_ra = cat_flt[selection,key_img[coord[0]]]
                match_dec = cat_flt[selection,key_img[coord[1]]]
                match_delta = delta[selection]
                match_mag = cat_flt[selection,key_img[mag]]
        else:
            print 'slitcatmatch: No catalog matches were found for this trace.'
            print 'Slit {0} {1}'.format(slit_i,which_trace)
            match_id = match_ra = match_dec = match_delta = match_mag = -99
    elif j > 1:
        cat_flt = cat_flt[delta<tolerance,:]
        delta = delta[delta<tolerance]
        print 'slitcatmatch: More than one matches satisfy the separation tolerence.'
        print 'Slit {0} {1}'.format(slit_i,which_trace)
        print 'Match\tRA\t\tdec\tSeparation (arcsec)\tMagnitude'
        for k in range(j):
            print '{0}\t{1:0.5f}\t{2:0.4f}\t{3:0.3f}\t{4:0.1f}'.format(k,cat_flt[k,key[coord[0]]],cat_flt[k,key[coord[1]]],delta[k],cat_flt[k,key[mag]])
            #display a region at the object's location, with label
            cmd = '"fk5; circle point {0:0.6f} {1:0.5f}'.format(cat_flt[k,key[coord[0]]],cat_flt[k,key[coord[1]]])+' # color=red text{'+'{0}'.format(k)+'}"'
            d.set('regions', cmd)            
        print '{0}\tSelect none.'.format(j)
        selection = raw_input('Enter the number of the correct match: ')
        if numpy.size(numpy.arange(k+1)==int(selection))==0:
            selection = rawinput("Input invalid. Please enter a valid number.: ")
        if selection == str(j):
            # Don't associate the trace with an object
            match_id = match_ra = match_dec = match_delta = match_mag = -99
        elif numpy.size(numpy.arange(k+1)==int(selection))!=0:
            selection=int(selection)
            match_id = cat_flt[selection,key_img[objkey]]
            match_ra = cat_flt[selection,key_img[coord[0]]]
            match_dec = cat_flt[selection,key_img[coord[1]]]
            match_delta = delta[selection]
            match_mag = cat_flt[selection,key_img[mag]]
#    print maskname
#    print slit_i
#    print which_trace
#   print y/pixscale
#    print ra_trace
#    print dec_trace
#    print match_id
#    print match_ra
#    print match_dec
#    print match_delta

    fh = open(outputfile,'a')
    fh.write('{0:0.0f}\t{1:0.6f}\t{2}\t{3:0.0f}\t{4}\t{5:03d}\t{6}\t{7:0.1f}\t{8:0.6f}\t{9:0.5f}\t{10:0.0f}\t{11:0.6f}\t{12:0.5f}\t{13:0.2f}\t{14}\n'.format(obj,z,zerr,quality,maskname,slit_i,which_trace,y/pixscale,ra_trace,dec_trace,match_id,match_ra,match_dec,match_delta,slitcomment))
    fh.close()
    
# Analyze each trace for all slits and associate with photometric object
for slit_i in slitnumbers:
    #check if object is a science target or just an alignment star
    slittyp = tb_slits.field('SLITTYP')[tb_slits.field('SLITNAME')==slit_i]
    if slittyp == 'A':
        continue
    elif slittyp != 'P':
        print 'slitcatmatch: Error unexpected SLITTYP for slit {0}'.format(slit_i)
        print 'SLITTYP = "P" expected but {0}. Will still try to match object.'.format(slittyp)

    #Filter the tables keeping only the current slit
    tb_s = tb_slits[tb_slits.field('SLITNAME')==slit_i]
    # slitid
    dslitid = tb_s.field('DSLITID')
    tb_m = tb_map[tb_map.field('DSLITID')==dslitid]
    # objectid
    objectid = tb_m.field('OBJECTID')
    tb_t = tb_targets[tb_targets.field('OBJECTID')==objectid]
    # object (e.g.: DLS photometric objid)
    obj = tb_t.field('OBJECT')[0]

    # check if there is a 1d trace for this slit
    try:
        hdutrace = pyfits.open(path+'spec1d.{0}.{1}.{2}.fits'.format(maskname,slit_i,obj))
    except IOError:
        print 'slitcatmatch: There is no spec1d trace file for slit number {}'.format(slit_i)
        continue
    else:
        match(slit_i,'primary',cat_img,key_img,imgcoord,objkey,mag,tolerance,outputfile)
    
    # check if there is a serendip trace and if so then attempt to match with catalog
    # works for up to 5 serendips
    for i in range(1,6):
        try:
            hdutrace = pyfits.open(path+'spec1d.{0}.{1}.serendip{2}.fits'.format(maskname,slit_i,i))
        except IOError:
            continue
        else:
            # the file exists, try to match the object
            match(slit_i,'serendip{}'.format(i),cat_img,key_img,imgcoord,objkey,mag,tolerance,outputfile)
