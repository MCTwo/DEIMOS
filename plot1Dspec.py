'''
This function is designed to plot reduced DEIMOS 1D spectra.
'''
from __future__ import division
import numpy
import pylab
import pyfits
import sys

## User Input
#datapath = '/sandbox/deimos/1rxs3A/2013sep05/'
#maskname = '1rxs3A'
#slit = 23 # three digit slit number that you want to plot
#redshift = 0.23507
#which_trace = 0
#pixbin = 10

def plot1D(datapath,maskname,slit,redshift,which_trace=0,pixbin=10):
    ###########################
    ### PROGRAM
    ###########################
    binfile = maskname+'.bintabs.fits'
    # Gather the basic slit info tables from the bintabs.fits file
    hdubin = pyfits.open(datapath+binfile)
    # Target table information, similar to dsim .lst information
    tb_targets = hdubin[1].data
    # Table of slitmask information
    tb_mask = hdubin[2].data
    # Table of each slit's physical properties
    tb_slits = hdubin[3].data
    # Table that maps the id's of the previous two tables
    tb_map = hdubin[4].data
    
    #Filter the tables keeping only the current slit
    if slit < 10:
        slit = '00'+str(slit)
    elif slit < 100:
        slit = '0'+str(slit)
    else:
        slit = str(slit)
    tb_s = tb_slits[tb_slits.field('SLITNAME')==slit]
    # slitid
    dslitid = tb_s.field('DSLITID')
    tb_m = tb_map[tb_map.field('DSLITID')==dslitid]
    # objectid
    objectid = tb_m.field('OBJECTID')
    tb_t = tb_targets[tb_targets.field('OBJECTID')==objectid]
    # object (e.g.: DLS photometric objid)
    obj = tb_t.field('OBJECT')[0]
    
    # check if there is a 1d trace for this slit
    if which_trace == 0:
        try:
            hdutrace = pyfits.open(datapath+'spec1d.{0}.{1}.{2}.fits'.format(maskname,slit,obj))
        except IOError:
            print 'plot1Dspec: Error, there is no spec1d trace file for slit number {}, exiting'.format(slit)
            sys.exit()
        else:
            hdutrace = pyfits.open(datapath+'spec1d.{0}.{1}.{2}.fits'.format(maskname,slit,obj))
            fig_title = 'spec1d.{0}.{1}.{2}.fits'.format(maskname,slit,obj)
        
    if which_trace > 0:
        try:
            hdutrace = pyfits.open(datapath+'spec1d.{0}.{1}.serendip{2}.fits'.format(maskname,slit,which_trace))
        except IOError:
            print 'plot1Dspec: Error, there is no serendip{0} spec1d trace file for slit number {1}, exiting'.format(which_trace,slit)
            sys.exit()
        else:
            # the file exists, read in the data
            hdutrace = pyfits.open(datapath+'spec1d.{0}.{1}.serendip{2}.fits'.format(maskname,slit,which_trace))
            fig_title = 'spec1d.{0}.{1}.serendip{2}.fits'.format(maskname,slit,which_trace)
    
    # read in the blue side data
    blue = hdutrace[3].data # this grabs the Horne extraction
    spec_b = blue[0][0] #trace flux
    lambda_b = blue[0][1] #trace observed wavelendth
    
    # read in the red side data
    red = hdutrace[4].data # this grabs the Horne extraction
    spec_r = red[0][0]
    lambda_r = red[0][1]
    
    # get rid of the zeroed out spectrum values (typically just at the very ends of
    # the blue and red sides
    mask_b = spec_b > 0
    mask_r = spec_r > 0
    spec_b = spec_b[mask_b]
    lambda_b = lambda_b[mask_b]
    spec_r = spec_r[mask_r]
    lambda_r = lambda_r[mask_r]
    
    # trim off the first and last 10 pixels of the spectra since these are usually
    # just noisy pixels
    spec_b = spec_b[10:]
    lambda_b = lambda_b[10:]
    spec_r = spec_r[:-10]
    lambda_r = lambda_r[:-10]
    
    
    
    if pixbin > 1:
        # convolve the spectra with a boxcar of size pixbin
        boxcar = numpy.ones(pixbin)
        spec_b = numpy.convolve(spec_b,boxcar)    
        spec_r = numpy.convolve(spec_r,boxcar)
        # trim the new edges from the convolved spectra
        # determine the added number of pixels during the convolution
        pix_added = pixbin-1
        trim_left = pix_added // 2
        trim_right = pix_added - trim_left
        spec_b = spec_b[trim_left:-trim_right]
        spec_r = spec_r[trim_left:-trim_right]
            
    fig = pylab.figure(figsize=(20,5))
    pylab.plot(lambda_b,spec_b,'k')
    pylab.plot(lambda_r,spec_r,'k')
    xl = pylab.xlim()
    yl = pylab.ylim()
    #Plot common spectral lines
    x_Hb = 4861*(1+redshift)
    x_OIII_1 = 4960*(1+redshift)
    x_OIII_2 = 5008*(1+redshift)
    x_Mgb = 5176*(1+redshift)
    x_FeI = 5269*(1+redshift)
    x_NaD = 5893*(1+redshift)
    x_NII_1 = 6548*(1+redshift)
    x_Ha = 6563*(1+redshift)
    x_NII_2 = 6585*(1+redshift)
    x_SII = 6726*(1+redshift)
    
    
    pylab.plot((x_Hb,x_Hb),yl,'--k')
    pylab.plot((x_OIII_1,x_OIII_1),yl,'--k')
    pylab.plot((x_OIII_2,x_OIII_2),yl,'--k')
    pylab.plot((x_Mgb,x_Mgb),yl,'--k')
    pylab.plot((x_FeI,x_FeI),yl,'--k')
    pylab.plot((x_NaD,x_NaD),yl,'--k')
    pylab.plot((x_NII_1,x_NII_1),yl,'--k')
    pylab.plot((x_Ha,x_Ha),yl,'--k')
    pylab.plot((x_NII_2,x_NII_2),yl,'--k')
    
    pylab.text(x_Hb, 0.05*(yl[0]+yl[1]), 'Hb', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    pylab.text(x_OIII_1, 0.05*(yl[0]+yl[1]), '[OIII]', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    pylab.text(x_OIII_2, 0.05*(yl[0]+yl[1]), '[OIII]', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    pylab.text(x_Mgb, 0.05*(yl[0]+yl[1]), 'Mg I(b)', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    pylab.text(x_FeI, 0.05*(yl[0]+yl[1]), 'Fe I', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    pylab.text(x_NaD, 0.05*(yl[0]+yl[1]), 'Na I (D)', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    pylab.text(x_NII_1, 0.05*(yl[0]+yl[1]), '[NII]', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    pylab.text(x_Ha, 0.05*(yl[0]+yl[1]), 'Ha', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    pylab.text(x_NII_2, 0.05*(yl[0]+yl[1]), '[NII]', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    
    
    
    
    pylab.xlabel('$\lambda_{observed}$')
    pylab.ylabel('Flux')
    pylab.title(fig_title)
    
    ## Flats
    ## check if there is a blue side calibration for this slit, and read it in
    #try:
        #hducalib_b = pyfits.open(datapath+'calibSlit.{0}.{1}B.fits.gz'.format(maskname,slit))
    #except IOError:
        #print 'plot1Dspec: Error, there is no blue side calibration file for slit number {}, exiting'.format(slit)
        #sys.exit()
    #else:
        #hducalib_b = pyfits.open(datapath+'calibSlit.{0}.{1}B.fits.gz'.format(maskname,slit))
    ## check if there is a red side calibration for this slit, and read it in
    #try:
        #hducalib_r = pyfits.open(datapath+'calibSlit.{0}.{1}R.fits.gz'.format(maskname,slit))
    #except IOError:
        #print 'plot1Dspec: Error, there is no red side calibration file for slit number {}, exiting'.format(slit)
        #sys.exit()
    #else:
        #hducalib_r = pyfits.open(datapath+'calibSlit.{0}.{1}R.fits.gz'.format(maskname,slit))
    
    #blue_flat = hducalib_b[1].data
    #red_flat = hducalib_r[1].data
    
    #flat_b = blue_flat[0][0]
    #flat_r = red_flat[0][0]
    
    #Nlam_b = numpy.size(lambda_b)
    #Nlam_r = numpy.size(lambda_r)
    
    #dlam_b = (lambda_b[-1]-lambda_b[0])/Nlam_b
    #dlam_r = (lambda_r[-1]-lambda_r[0])/Nlam_r
    
    #lam1_b = lambda_b[-1]+dlam_b
    
    pylab.ylim(yl)
    pylab.show()
    
