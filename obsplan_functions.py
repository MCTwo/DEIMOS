#-------------------------------------------------------------------------------
#Purpose: This is a list of functions for making DEIMOS slitmasks
#Originally written by Will. A. Dawson <will@dawsonresearch.com>
#Rewritten by Karen Y. Ng <karen.yng@ucdavis.edu>
#License: BSD
#Date: 06/07/2013
#-------------------------------------------------------------------------------
'''
This is meant to be used with 
* SDSS catalog 
* stored in the format of a dataframe

'''
from __future__ import division
import sys
import numpy as np,numpy
import pandas as pd 
import matplotlib.pyplot as plt
import numpy.ma as ma
import pyfits 
import aplpy
import pdb

def filter_catalog_dataframe(cat,field,lowerbound=None, upperbound=None, 
                            plot_diag=False,
                            save_diag=False, verbose=True):
    '''
    Stability : works 
    INPUT:
    cat = dataframe containing your data
    field = string denoting the field name 
    lowerbound = some number of the bound
    lbound_
    OUTPUT:
    filtered_cat = filtered catalog in terms of a dataframe
    '''
    dataNo = cat.shape[0]
    if (upperbound==None and lowerbound==None):
        print 'WARNING: both upper and lower bound are not specified'
        return cat
    if (upperbound!=None and lowerbound!=None):
        mask = np.logical_and(cat[field]>=lowerbound, cat[field]<= upperbound)
        cat = cat[mask]
        return cat
    else:
        if lowerbound!=None:
            mask = cat[field] >= lowerbound 
            cat = cat[mask]
        if upperbound!=None:
            mask = cat[field] <= upperbound
            cat = cat[mask]
    if plot_diag == True:
        plt.xlabel('Data')
        plt.ylabel(field)
        plt.plot(cat[field],'x')
        plt.show()
        if save_diag==True:
            plt.save('diagnose_filter_'+field+'.png')
    print 'Filtering according to '+field
    print 'Started out with %d\n filtered out %d data entries\n with %d remaining'%\
            (dataNo,dataNo-cat.shape[0],cat.shape[0])

    return cat 

def determine_weight(cat, x, sigma_x, mu_cl, weight_field='weight', 
                     plot_diag=False):
    '''
    Stability : works 
    INPUT:
    x = estimated photoz for the galaxy 
    sigma_x = estimated error of photoz for the galaxy 
    mu_cl = estimated z for the cluster  
    OUTPUT:
    weight(x=z_cl) = 1 + 100 * N(z_gal, sigma_gal, x = z_cl)
    '''
    cat[weight_field]= 1+100*np.exp(-(x-mu_cl)**2/(2*sigma_x**2))/np.sqrt(2*np.pi)/ sigma_x
    if plot_diag==True:
        plt.ylabel('Weight')
        plt.xlabel('Data')
        plt.hist(cat[weight_field],bins=100)
        plt.show()
    return cat

def determine_sample_no(cat, sample_no, field, first_sample = False, 
                        lower_criteria=None,
                        upper_criteria=None,
                         plot_diag=False, save_diag=False,
                        verbose=True):
    '''
    Description of how slits are selected based on sample no: 
    Sample -- (int) 1, 2, 3, etc. Sample to which object belongs.
    When auto-selecting, objects in Sample 1 are selected first; 
    remaining space is then filled with Sample 2, then Sample 3, etc. 
    Default=1
    (http://www.ucolick.org/~phillips/deimos_ref/masks.html#format)
    
    stability: works
    INPUT:
    cat = data frame catalog object 
    field = (field) criteria which sample number is selected based on 
    first_sample = determine if this is the first sample number being
    initialized 
    lower_critiera, uppercriteria  = criteria between which the data
    is going to assign this particular sample no 
    sample_no = int to be assigned 

    OUTPUT:
    cat = data frame with an extra field called "sample"

    '''
    #create a new field for storing sample number
    if first_sample ==True:
        cat['sample'] = pd.Series(3*np.ones(cat.shape[0]),index=cat.index)
    
    if(lower_criteria ==None and upper_criteria ==None):
        print 'Warning : Both upper and lower critiera are NA'
        return cat
    if (lower_criteria !=None and upper_criteria !=None):
        mask = np.logical_and(cat[field]>=lower_criteria,
                              cat[field]<= upper_criteria)
        cat['sample'][mask] = sample_no
        return cat
    else:
        if lower_criteria != None:
            mask = cat[field] >= lower_criteria
            cat['sample'][mask] = sample_no
        if upper_criteria != None:
            mask = cat[field] <= upper_criteria
            cat['sample'][mask] = sample_no
    
    if plot_diag ==True:
        plt.xlabel('Sample')
        plt.ylabel(field)
        plt.plot(cat['sample'],cat[field],'x')
        plt.show()
        if save_diag ==True:
            plt.savefig('sample_vs_'+field+'_cut.png')

    return cat 

def plot_gal_histogram_as_fits(cat, prefix, rabin=50,
                               rarange=None, decrange=None, zrange=None,
                               N_boot = None, plot_hist=False,
                               plot_contour=False):
    
    ra_min = cat['ra'].min(axis=1)
    ra_max = cat['ra'].max(axis=1)

    dec_min = cat['dec'].min(axis=1)
    dec_max = cat['dec'].max(axis=1)



    #The rest of this code is adapted from the number density function 
    dec_mean = (dec_max + dec_min) / 2 * numpy.pi / 180.0 #radians
    # determine decbin such that the pixes are approximately square
    decbin = rabin*(dec_max-dec_min)//((ra_max-ra_min)*numpy.cos(dec_mean))
    # create the pixel/bin edge arrays
    ra_binwidth = (ra_max-ra_min)/rabin
    dec_binwidth = (dec_max-dec_min)/decbin
    ra_edge = numpy.arange(ra_min,ra_max+ra_binwidth,ra_binwidth)
    dec_edge = numpy.arange(dec_min,dec_max+dec_binwidth,dec_binwidth)
    # due to possible rounding issues it is necessary to redefince decbin to
    # match dec_edge
    decbin = numpy.size(dec_edge)-1

    # Create the blank map_array
    if N_boot != None:
        h = numpy.zeros((3+N_boot,decbin,rabin))
        N = numpy.shape(cat)[0] #number of rows in filtered catalog
        # Random with replacement bootstrap index array 
        b = numpy.random.randint(0,high=N,size=(N,N_boot))
        print 'numberdensity: will perform bootstrap analysis with {0} random iterations'.format(N_boot)
            
    ##Convert the dataframe to numpy array for the next steps
    dec = np.array(cat['dec'])
    ra = np.array(cat['ra'])
            
    #create the 2D histogram    
    if N_boot == None:
        h, tmp_edge1, tmp_edge2 = numpy.histogram2d(dec,ra,
                                                 bins=(dec_edge,ra_edge))
        plt.xlabel('DEC')
        plt.ylabel('RA')
        plt.imshow(h)
        plt.colorbar()
        if plot_hist==True:
            plt.title('Not perfect hist, DEC & RA might have flipped')
            plt.show()
    else:
        h[0,:,:], tmp_edge1, tmp_edge2 = numpy.histogram2d(cat[:,deccol],cat[:,racol],bins=(dec_edge,ra_edge))
        for i in numpy.arange(N_boot):
            h[3+i,:,:], tmp_edge, tmp_edge = numpy.histogram2d(cat[b[:,i],deccol],cat[b[:,i],racol],bins=(dec_edge,ra_edge))
        #Create signal/noise and standard deviation maps
        h[2,:,:] = numpy.std(h[3:,:,:],axis=0,ddof=1)
        h[1,:,:] = h[0,:,:]/h[2,:,:]
    #Create the fits file
    H = h
    hdu = pyfits.PrimaryHDU(numpy.float32(H))

    # Calculate the wcs CR**** for the fits header    
    xscale = (ra_max - ra_min) * numpy.cos(dec_mean) / rabin
    yscale = (dec_max - dec_min) / decbin
    crval1 = (ra_max - ra_min) / 2 + ra_min
    crval2 = (dec_max - dec_min) / 2 + dec_min
    crpix1 = rabin / 2
    crpix2 = decbin / 2

    # Apply the wcs to the fits header
    hdu.header.update('ctype1', 'RA---TAN')
    hdu.header.update('ctype2', 'DEC--TAN')
    hdu.header.update('crval1', crval1)
    hdu.header.update('crval2', crval2)
    hdu.header.update('crpix1', crpix1)
    hdu.header.update('crpix2', crpix2)
    hdu.header.update('cd1_1',xscale)
    hdu.header.update('cd1_2',0)
    hdu.header.update('cd2_1',0)
    hdu.header.update('cd2_2',yscale)
    filename = prefix+'_masked_nodensity'
    hdu.writeto(filename+'.fits',clobber=True) 

def convert_to_sexadec_coord(cat,ra_field='ra',dec_field='dec',
                            dec_sign_field='sign',out_rah_field='rah',
                             out_ram_field='ram', out_ras_field='ras',
                             out_decd_field='decd', out_decm_field='decm',
                             out_decs_field='decs'):
    '''
    INPUT: 
    cat = dataframe catalog object 
    ra_field = string denoting the ra field name
    dec_field = string denoting the dec field name
    similar with the 
    out_*_field = string denoting the converted * field name 

    OUTPUT:
    cat = dataframe but with several more fields 
    the wildcard * denotes the following:
    * 'rah' degrees of ra
    * 'ram' minutes of ra 
    * 'ras' seconds of ra 
    * 'decd' degrees of dec
    * 'decm' minutes of dec 
    * 'decs' seconds of dec 
    * 'sign' sign of dec    
    '''
    #convert deg RA to sexadec RA
    # I should have used deg2ra instead 
    # I should have used deg2dec
    # but then I also have to change write_galaxies 
    ra = cat[ra_field]/15.0
    cat[out_rah_field] = np.floor(ra)
    res = (ra-cat[out_rah_field])*60.
    cat[out_ram_field] = np.floor(res)
    cat[out_ras_field] = (res-cat[out_ram_field])*60.

    #convert deg dec to sexadec dec
    dec = cat[dec_field]
    #initialize a field for storing the sign of dec
    cat[dec_sign_field] = cat[dec_field]

    for i in cat.index:
        if cat[dec_field][i]<0:
            cat[dec_sign_field][i]= -1.
            dec = abs(dec)
        else:
            cat[dec_sign_field][i] = 1.
    cat[out_decd_field] = np.floor(dec)
    res = (dec-cat[out_decd_field])*60.
    cat[out_decm_field] = np.floor(res)
    cat[out_decs_field] = (res-cat[out_decm_field])*60.
    return cat

def return_objects_in_mask_region(cat,regfile):
    '''
    INPUT:
    cat = dataframe catalog object 
    
    regfile = 
    mask file generated from ds9 see warning statement below about 
    how this mask file should be created 

    OUTPUT:
    cat = dataframe with only the objects in the box region 
    box = info specifying how the mask is oriented 
    
    '''
    box = numpy.fromregex(regfile,
        r"box\(([0-9]*\.?[0-9]+),([0-9]*\.?[0-9]+),([0-9]*\.?[0-9]+)\",([0-9]*\.?[0-9]+)\",([0-9]*\.?[0-9]+)",
        [('xc',numpy.float),('yc',numpy.float),('width',numpy.float),
        ('height',numpy.float),('angle',numpy.float)])

    if len(box) == 0: 
        print 'WARNING!!!!! \
               When saving your reg file,choose fk5 > WCS > degrees\n \
               in the popup menu there is something wrong with your\n \
               box region file format~!\n \
               You have put a regfile with wrong coordinate format'

    d2r = numpy.pi/180.0
    #loop through the regions creating masks for galaxy inclusion
    for i in numpy.arange(numpy.shape(box)[0]):
        #phi is the ccw angle from the +East axis
        xc = box[i][0]
        yc = box[i][1]
        w = box[i][2]
        h = box[i][3]
        phi=box[i][4]*d2r
        #rotate the galaxies into the "primed" (p) region coorditate frame centered
        #at the center of the region
        ra_p = (np.array(cat['ra'])-xc)*numpy.cos(yc*d2r)*numpy.cos(-phi)+(np.array(cat['dec'])-yc)*numpy.sin(-phi)
        dec_p = -(np.array(cat['ra'])-xc)*numpy.cos(yc*d2r)*numpy.sin(-phi)+(cat['dec']-yc)*numpy.cos(-phi)
        #determine the min and max bounds of the region
        # min = (box center [deg])-((box height [sec])/(2*60**2))
        ra_p_min = -w/(2*60**2)
        ra_p_max = w/(2*60**2)
        dec_p_min = -h/(2*60**2)
        dec_p_max = h/(2*60**2)
        #create the mask for the i region
        mask_ramin = ra_p >= ra_p_min
        mask_ramax = ra_p <= ra_p_max
        mask_decmin = dec_p >= dec_p_min
        mask_decmax = dec_p <= dec_p_max
        mask_i = mask_ramin*mask_ramax*mask_decmin*mask_decmax
        
        #combine the i mask with the previous masks
        if i == 0:
            mask = mask_i
        else:
            #if the galaxy was in any mask then it should be in the concatenated mask
            mask_tmp = mask + mask_i
            mask = mask_tmp > 0

    #apply the region filter to the input galaxy catalog
    # see http://www.ucolick.org/~phillips/deimos_ref/masks.html for file format
    Nint = cat.shape[0]
    cat = cat[mask]
    Nfin = cat.shape[0]
    Ncut = Nint - Nfin
    print 'obsplan: {0} rows were removed from the catalog with\n \
           {1} initial rows, leaving {2} rows'.format(Ncut,Nint,Nfin)
    return cat,box 

def write_dsim_header(F,prefix,box):
    F.write('#This catalog was created by plan*.py and is intended to be used \n')
    F.write('#as input to the deimos slitmask software following the format \n')
    F.write('#outlined at http://www.ucolick.org/~phillips/deimos_ref/masks.html\n')
    F.write('#Note that the automatic generation of this file does not include\n')
    F.write('#guide or alignment stars.\n')
    F.write('#ttype1 = objID\n')
    F.write('#ttype2 = ra\n')
    F.write('#ttype3 = dec\n')
    F.write('#ttype4 = equinox\n')
    F.write('#ttype5 = magnitude\n')
    F.write('#ttype6 = passband\n')
    F.write('#ttype7 = priority_code\n')
    F.write('#ttype8 = sample\n')
    F.write('#ttype9 = selectflag\n')
    F.write('#ttype10 = pa_slit\n')
    F.write('#ttype11 = len1\n')
    F.write('#ttype12 = len2\n')
    #Write in the Slitmask information line
    F.write('{0}\t{1:0.6f}\t{2:0.6f}\t2000\tPA={3:0.2f}\n'
            .format(prefix,box[0][0]/15.,box[0][1],box[0][4]-90))

#def write_guide_stars(F, cat):
#    '''
#    Status: WORK IN PROGRESS
#
#    Guide stars are in the guide camera region with Preselect code = -1 
#    INPUT:
#    F = filestream for writing the file out to 
#    cat = dataframe catalog of the stars
#    
#    Feature to be implemented: with a use of pyds9 
#    we should be able to select and pick regions without leaving python
#    and outputting the list of selected stars directly for writing 
#    
#    #a sample line should look like
#    #F.write("176676037	10:54:07.305	54:56:42.927	2000	16.33	R	-1	0	1\n")
#    '''
#    for i in cat.index:
#        F.write("{0} {1}:{2} 2000	{3}	R	-1	0	1\n".format(
#            cat['objID'][i],cat['ra'][i],cat['dec'][i], cat['dered_r'][i]))
#
#    return

#def write_align_stars(filestream, align_star_cat):
#    ''''
#    Work in progress
#    #sample line 
#    #F.write("412932112	10:54:44.062	54:48:23.996	2000	17.91	R	-2	0	1\n")
#    '''
#    for i in cat.index:
#        F.write("{0} {1} {2} 2000	{3}	R	-2	0	1\n".format(cat['objID'][i],
#           cat['ra'][i],cat['dec'][i], cat['dered_r'][i]))
#    return

def write_galaxies_to_dsim(F, cat,  sky):
    '''
    Stability: works
    INPUT:
    F = file stream of the file to write to 
    cat = dataframe object of the galaxy 
    '''
    print '!!!!!WARNING!!! Chopping off first 3 digit of SDSS ObjID'
    print 'if you are not using SDSS, modify\
    obsplan_functions.write_galaxies_to_dsim() to disable this'
    cat['shortID'] = cat['objID']-1230000000000000000
    #current dsim input can only accept obj name limited to 16 characters 
    #SDSS ObjID has 18 characters 
    #instead of writting ObjID out we write the index out
    for i in cat.index:
        if cat['sign'][i]==-1:
            F.write('{0:16d}\t{1:02.0f}:{2:02.0f}:{3:06.3f}\t-{4:02.0f}:{5:02.0f}:{6:06.3f}\t2000\t{7:0.2f}\tR\t{8:0.0f}\t{9:0.0f}\t{10:0.0f}\t{11:0.2f}\t{12:0.1f}\t{13:0.1f}\n'
                        .format(cat['shortID'][i],cat['rah'][i],
                                cat['ram'][i],
                                cat['ras'][i], 
                                cat['decd'][i],
                                cat['decm'][i],
                                cat['decs'][i],
                                cat['dered_r'][i],
                                cat['weight'][i],
                                cat['sample'][i],
                                cat['pscode'][i]*1.,
                                cat['PA'][i],
                                cat['deVRad_r'][i]/2.+sky[0],
                                cat['deVRad_r'][i]/2.+sky[1]))
        else:
            F.write('{0:16d}\t{1:02.0f}:{2:02.0f}:{3:06.3f}\t{4:02.0f}:{5:02.0f}:{6:06.3f}\t2000\t{7:0.2f}\tR\t{8:0.0f}\t{9:0.0f}\t{10:0.0f}\t{11:0.2f}\t{12:0.1f}\t{13:0.1f}\n'
                        .format(cat['shortID'][i],cat['rah'][i],
                                cat['ram'][i],
                                cat['ras'][i], 
                                cat['decd'][i],
                                cat['decm'][i],
                                cat['decs'][i],
                                cat['dered_r'][i],
                                cat['weight'][i],
                                cat['sample'][i],
                                cat['pscode'][i]*1,
                                cat['PA'][i],
                                cat['deVRad_r'][i]/2.+sky[0],
                                cat['deVRad_r'][i]/2.+sky[1]))

    return

def write_mask_for_dsim(gal_cat, guide_star_cat, align_star_cat, mask_reg, 
                        dsim_outfile):
    '''
    INPUT:
    gal_cat = catalog dataframe object containing galaxies in the mask region
    star_cat = catalog dataframe object containing guide stars in mask region
    align_star_cat = catalog dataframe object containing alignment stars 
    mask_reg = string denoting the name of the mask region file
    dsim_outfile = string denoting the name of the output file 

    OUTPUT:
    file as input for dsim (an IRAF module on Theta)
    '''
    F = open(dsim_outputfile, 'w')
    write_dsim_header(F)
    write_guide_stars(F,guide_star_cat)
    write_align_stars(F,align_star_cat)
    write_galaxies(F,gal_cat)
    F.close()

    return

def output_target_as_exclude_csv():
    '''
    write out target of one mask as input of exclusion file
    '''


    return

def PAround(cat,field,PAmin,PAmax,PAvalue,maskPA):
    '''
    Inspects each element of the PAarray and 
    if it falls between PAmin and PAmax
    then it redefines the PA to the PAround value.
    Where the min and max bounds are in the masks coordinate system.
    Purpose : Determine the optimal PA for the slits
    '''



    maskPAmin = PAmin+maskPA
    maskPAmax = PAmax+maskPA
    if maskPAmax > 90:
        test1a = cat[field] >= maskPAmin
        test1b = cat[field] <= maskPAmax
        test1 = test1a*test1b
        test2a = cat[field] >= maskPAmin-180
        test2b = cat[field]<= maskPAmax-180
        test2 = test2a*test2b
        test = test1+test2
    elif maskPAmin <90:
        test1a = cat[field]>= maskPAmin
        test1b = cat[field]<= maskPAmax
        test1 = test1a*test1b
        test2a = cat[field] >= maskPAmin+180
        test2b = cat[field] <= maskPAmax+180
        test2 = test2a*test2b
        test = test1+test2
    else:
        test1a = cat[field] >= maskPAmin
        test1b = cat[field] <= maskPAmax
        test = test1a*test1b
    cat[field][test]=PAvalue+maskPA
    return cat

def pick_PA(cat, PA_field, box, axis_angle='deVPhi_r',plot_diag=False):
    '''
    Stability: it runs  
    Modified from the original obsplan.py PA function
    INPUT angle
    axis_angle = the angle of the major axis?
                 by default use sdss name for this angle
    '''
    #maskPA is the orientation of the mask
    #that complies with dsim 's definition of the angles
    maskPA = box[0][4]+90
    if maskPA > 90:
        maskPA -= 180
    if maskPA < -90:
        maskPA += 180

    
    ### Attempt to align the slits with the major axis of each galaxy
    ##PA = cat[:,key['deVPhi_r']]*1.0
    # Attempt to align the slits with the minor axis of each galaxy
    cat[PA_field] = cat[axis_angle]*1.0+90

    # make slits to have angle between 5 and 30(-150) degrees 
    # or between -5 and -30(150) degrees 

    #The 0-5 spec is based on DEEP2 recommendations for
    #better sky subtraction
    cat = PAround(cat, PA_field,0,5,5,maskPA)
    cat = PAround(cat, PA_field,-5,0,-5,maskPA)
    #cat = PAround(cat, PA_field,175,180,175,maskPA)
    #cat = PAround(cat, PA_field,-180,-175,-175,maskPA)

    #make slits NOT to lie perpendicular to the long axis of the mask
    cat = PAround(cat, PA_field,30,90,30,maskPA)
    cat = PAround(cat, PA_field,90,150,150,maskPA)
    cat = PAround(cat, PA_field,-90,-30,-30,maskPA)
    cat = PAround(cat, PA_field,-150,-90,-150,maskPA)
    cat = PAround(cat, PA_field,210,270,210,maskPA)
    cat = PAround(cat, PA_field,270,330,330,maskPA)
    cat = PAround(cat, PA_field,-210,-270,-210,maskPA)
    cat = PAround(cat, PA_field,-270,-330,-330,maskPA)


    if plot_diag ==True:
        plt.axhline(y=5,color='r')
        plt.axhline(y=-5,color='r')
        plt.axhline(y=30,color='r')
        plt.axhline(y=-30,color='r')
        plt.axhline(y=185, color='r')
        plt.axhline(y=175,color='r')
        plt.axhline(y=210,color='r')
        plt.axhline(y=150,color='r')

        plt.ylim(-360,360)
        plt.ylabel('PA angle')
        plt.xlabel('Data index')
        plt.plot(cat[PA_field]-maskPA,'x')
        plt.show()
    return cat

def exclude_objects(cat, exclude_file):
    '''
    
    '''
    #actually we can selectively import columns using pandas to make read_csv
    #faster 
    #i just didn't bother to figure this out
    col_name = ['objID','sex_ra','sex_dec','equinox','dered_r','R','weight',
               'sample','pscode','stuff1','stuff2','stuff3']
    exclude_cat = pd.read_csv(exclude_file,skiprows=7,delimiter=r"\s*",names=col_name)
    print 'If you see error message for the reading in of exclude_cat'
    print 'check if there are non-data rows in the middle of the file'
    exclude_cat = exclude_cat.dropna()

    print 'obsplan exclude_objects: apply exclusion list to further filter catalog'
    mask_ex = numpy.zeros(cat.shape[0])
    i = 0
    for i in exclude_cat.index:
        oid = int(exclude_cat['objID'][i])
        if i == 0:
            mask_ex = cat['objID']-1230000000000000000  == oid 
            i=1
        else:
            mask_i = cat['objID']-1230000000000000000  == oid
            mask_tmp = mask_ex+mask_i
            mask_ex = mask_tmp > 0
    mask_ex = mask_ex == False
    Nint = cat.shape[0]
    cat = cat[mask_ex]
    Nfin = cat.shape[0]
    Ncut = Nint - Nfin
    print 'obsplan exclude_objects:'
    print '{0} rows were removed from the catalog with {1} initial rows, leaving {2} rows'.format(Ncut,Nint,Nfin)
    return cat

#def output_candidate_star_to_dsim(cat, output_prefix=None, F=None):
#    '''
#    Status: work in progress
#    Write the lines of the stars to dsim input format
#    To add more documentation
#    '''
#    if output_prefix == None and F == None:
#        print 'You have to specify where to output your stars'
#    if output_prefix != None: 
#        output = output_prefix+'.txt'
#        print 'Writing list of candidate stars to file '+output
#        F=open(output,'w')
#
#    #To add warning if both prefix and F is present
#        
#    for i in cat.index:
#        obj = cat['objID']
#        ra = cat['ra']
#        ra = str(tools.deg2ra(ra,":"))
#        dec = cat['dec']
#        dec = str(tools.deg2dec(dec,":"))
#        #rmag = cat[i,key['dered_r']]
#        output= 'F.write("'+'{0:4.0f}\t'.format(obj)+ra+'\t'+ dec+'\t2000\t{0:2.2f}'.format(cat['dered_r'])+'\tR\t-2\t0\t1' +r'\n")'+'\n'
#        F.write(output)
#
#    if prefix !=None: 
#        F.close()

def write_circle_reg(cat,output_prefix):
    '''
    Status: work in progress
    '''
    outputname = output_prefix+'_circles.reg'
    F = open(outputname,'w')
    F.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+'\n')
    F.write('fk5'+'\n')
    for i in cat.index:
        ra = cat['ra'][i]
        dec = cat['dec'][i]
        size = cat['deVRad_r'][i]+1
        obj = cat['objID'][i]
        F.write('circle({0:1.5f},{1:1.5f},{2:1.1f}") # text={{'.format(ra,dec,size)+'{0:0.0f}'.format(obj)+'}\n')
    F.close()

def write_slit_reg(cat,output_prefix,sky):
    outputname = output_prefix+'_slits.reg'
    print "reg file with name of:",outputname,"has been written."
    F = open(outputname,'w')
    F.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+'\n')
    F.write('fk5'+'\n')
    ###
    ### Survey Galaxies
    ###
    for i in cat.index:
        ra = cat['ra'][i]
        dec = cat['dec'][i]
        height = cat['deVRad_r'][i] +sky[0]+sky[1]
        #there is a 90 degree discrepancy between dsim and ds9 definitions
        #of angles
        angle = cat['PA'][i] +90
        if cat['sample'][i] == 1:
            color = 'green'
        else:
            color = 'blue'
        F.write('box({0:1.5f},{1:1.5f},{2:1.1f}",1",{3:0.2}) # color={4}'.format(ra,dec,height,angle,color)+'\n')
       
    F.close()

