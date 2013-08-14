from __future__ import division
import numpy
import tools
import sys

def readMaskRegion(regfile):
    box = numpy.fromregex(regfile,r"box\(([0-9]*\.?[0-9]+),([0-9]*\.?[0-9]+),([0-9]*\.?[0-9]+)\",([0-9]*\.?[0-9]+)\",([0-9]*\.?[0-9]+)",
                          [('xc',numpy.float),('yc',numpy.float),('width',numpy.float),('height',numpy.float),('angle',numpy.float)])
    return box[0]

def createSlitmaskMask(regfile,ra,dec):
    '''
    Input:
    regfile = [string] name of the ds9 region file. Note region should be
        defined using Coordinate/WCS/Degrees and Size/WCS/Arcmin options, with
        the Size 5 by 16.1 arcmin, Angle will then correspond to the slitmask's
        parallactic angle (i.e. +CCW from north towards east) with the guider
        camera in the North-east quadrent at Angle=0.
    ra = [1D array of floats; units:degrees] RA of the galaxies
    dec = [1D array of floats; units:degrees] Dec of the galaxies
    '''
    box = readMaskRegion(regfile)
    d2r = numpy.pi/180.0
    #phi is the ccw angle from the +North axis
    xc = box[0]
    yc = box[1]
    w = box[2]
    h = box[3]
    phi=box[4]*d2r
    #rotate the galaxies into the "primed" (p) region coorditate frame centered
    #at the center of the region
    ra_p = (ra-xc)*numpy.cos(yc*d2r)*numpy.cos(-phi)+(dec-yc)*numpy.sin(-phi)
    dec_p = -(ra-xc)*numpy.cos(yc*d2r)*numpy.sin(-phi)+(dec-yc)*numpy.cos(-phi)
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
    mask_slitmask = mask_ramin*mask_ramax*mask_decmin*mask_decmax
    return mask_slitmask

def photoz_PriorityCode(z_cluster,gal_photoz,photo_z_err,plot_diag=False):
    '''
    INPUT:
    x = estimated photoz for the galaxy 
    sigma_x = estimated error of photoz for the galaxy 
    mu_cl = estimated z for the cluster  
    OUTPUT:
    weight(x=z_cl) = 1 + 100 * N(z_gal, sigma_gal, x = z_cl)
    '''
    wght = 1+100*numpy.exp(-(gal_photoz-z_cluster)**2/(2*photo_z_err**2))/numpy.sqrt(2*numpy.pi)/ photo_z_err
    if plot_diag==True:
        plt.ylabel('Weight')
        plt.xlabel('Data')
        plt.hist(wght,bins=100)
        plt.show()
    return wght

def assignSample(param_array,samplebounds):
    '''
    Currently obsplan is only setup to break samples according to one object
    variable (e.g. magnitude).  The sampel_param_ttype is the ttype name of the
    vairable in the catalog to be used to make the sample division (e.g.
    magnitude). samplebounds defines the min and max of sample_param for each
    sample: e.g. (sample1 lowerbound, sample1 upperbound, sample2 lower bound,
    sample2 upper bound, etc., etc.).
    Input:
    param_array = [1D array of floats] a 1D array containing some parameter
        information for each target object (e.g. magnitude).
    samplebounds = (1D list of floats with length in multiple of 2) The sample
        bounds that define the min and max of values of param for an arbitrary
        number of samples, e.g: (sample1 lowerbound, sample1 upperbound
        sample2 lower bound,  sample2 upper bound, etc., etc.)
    '''
    N_gal = numpy.size(param_array)
    # Make sure that samplebounds in input in pair
    if numpy.size(samplebounds)%2 != 0:
        print 'obsplan.assignSample: error, samplebounds must contain an lower and upper bound for each sample. An odd number of bounds detected. Exiting.'
        sys.exit()
    N_sample = numpy.size(samplebounds)//2
    # if for some reason the sample bounds are not all inclusive, assign the by
    # default the next highest sample to all galaxies.
    sample = numpy.ones(N_gal)*(N_sample+1)
    # loop through each galaxy to see if it falls within a given sample bound
    for i in numpy.arange(N_sample):
        for j in numpy.arange(N_gal):
            # Check if in lowest priority sample first, then overwrite if makes
            # it in higher priority sample.
            if param_array[j] > samplebounds[2*N_sample-2-2*i] and param_array[j] < samplebounds[2*N_sample-1-2*i]:
                sample[j] = N_sample - i
    return sample
                

def assignSelectionFlag(objid,psfile=None,psobjid_ttype=None):
    '''
    Input:
    objid = [1D array int] object id's of the target catalog
    psfile = ['string'] ttype catalog of objects to preselect in dsim
    psobjid_ttype = ['string'] ttype name of the objid column in the psfile,
        the objid's should correspond to some of the objid's in the objid array
    Output:
    sflag = [1D array int] each object has value 1 if it is to be preselected in
        dsim or 0 if it is not to be preselected
    '''
    #sflag is the preselection code, should be 0 for not being preselected
    sflag = numpy.zeros(numpy.size(objid))    

    if psfile != None:
        #read in the preselection catalog and header
        pskey = tools.readheader(psfile)
        print 'obsplan: determining preselections'
        pslist = numpy.loadtxt(psfile,usecols=(pskey[psobjid_ttype],pskey[psobjid_ttype]))
        i = 0
        for oid in pslist[:,0]:
            if i == 0:
                sflag = objid == oid
                i=1
            else:
                mask_i = objid == oid
                sflag += mask_i
                i += 1    
        print 'obsplan: {0} slits preselected'.format(numpy.sum(sflag))
    # Since want a list of zeros and 1's need to convert bool array
    sflag = sflag*numpy.ones(numpy.size(sflag))
    return sflag

def createExclusionMask(objid,exfile,exobjid_ttype):
    '''
    Input:
    objid = [1D array int] object id's of the target catalog
    exfile = ['string'] ttype catalog of objects to exclude
    exobjid_ttype = ['string'] ttype name of the objid column in the exfile,
        the objid's should correspond to some of the objid's in the objid array
    Output:
    mask_ex = [1D boolean array] objects to be excluded with have a
        corresponding False value
    '''
    print 'obsplan: apply exclusion list to further filter catalog'
    exkey = tools.readheader(exfile)
    exlist = numpy.loadtxt(exfile,usecols=(exkey[exobjid_ttype],exkey[exobjid_ttype]))
    mask_ex = numpy.zeros(numpy.size(objid))
    i = 0
    for oid in exlist[:,0]:
        if i == 0:
            mask_ex = objid == oid
            i=1
        else:
            mask_i = objid == oid
            mask_tmp = mask_ex+mask_i
            mask_ex = mask_tmp > 0
    mask_ex = mask_ex == False
    return mask_ex

def objectPA(H,delta,phi=19.82525):
    '''
    This function calculates the paralactic angle (PA) of an astronomical object
    a the instant of its given hour angle (H), the declination of the object 
    (delta), and the geographical latitude of the observer (phi). The
    calculation is based on Equation 14.1 of Jean Meeus' Astronomical
    Algorithms (2nd edition).
    
    Input:
    H = [float; units=hours] object's hour angle (H is + if west of the meridian)
    delta = [float; units=degrees] object's declination
    phi = [float; units=degrees] observers latitude, defaults to Mauna Kea
    
    Output:
    parallactic angle of the object (-180 to 180 degrees) measured + from north
    ccw towards east
    '''
    from math import atan, tan, cos, sin, pi
    d2r = pi/180.
    phi *= d2r    
    H = H*15*d2r
    delta *= d2r
    if H < 0:
        sign = -1
        H = -H
    else:
        sign = 1  
    denom = tan(phi)*cos(delta)-sin(delta)*cos(H)
    q = atan(sin(H)/denom)
    q = q/d2r
    if denom < 0:
        q += 180
    return sign*q

def optimalPA(pa_mask,H,delta,phi=19.82525,relPA_min=5,relPA_max=30):
    '''
    As recommended by: 
    Filippenko, A.V., 1982. The importance of atmospheric differential refraction in spectrophotometry. Publications of the Astronomical Society of the Pacific, 94, pp.715-721. Available at: http://adsabs.harvard.edu/abs/1982PASP...94..715F.
    
    The spectral slit should be as much aligned with the axis along the
    horizon, object and zenith. Thus the position angle of the slit (as defined 
    by the angle + from north toward east) should equal the parallactic angle
    of the object. This is not always possible given the bounds placed on the
    slit orentation with respect to the slitmask, thus this funciton determines
    the best possible slit position angle with respect to the slitmask.
    
    Input:
    pa_mask = [float; units=degrees] parallactic angle of the mask    
    H = [float; units=hours] object's hour angle (H is + if west of the meridian)
    delta = [float; units=degrees] object's declination
    phi = [float; units=degrees] observers latitude, defaults to Mauna Kea
    relPA_min = [float; units=degrees] minimum absolute angle between the slit
       pa and the mask pa
    relPA_max = [float; units=degrees] maximum absolute angle between the slit
       pa and the mask pa
    '''
    from math import pi,sin,cos,asin
    # test that pa_mask is defined between 0 and 360 degrees
    test_pa_mask = numpy.logical_and(pa_mask >= 0, pa_mask <= 360)
    if ~test_pa_mask:
        print 'obsplan.optimalPA: error, mask_pa must be defined between 0 and 360 degrees,check that ds9 mask region is defined appropriately, exiting'
        sys.exit()
        
    # the optimal slit position angle is the parallactic angle of the object
    pa_obj = objectPA(H,delta,phi)
    
    # test to make sure that the pa_obj is in the range -180 to 180 degrees
    test_pa_obj = numpy.logical_and(pa_obj >=-180, pa_obj <= 180)
    if ~test_pa_obj:
        print 'obsplan.optimalPA: error, the pa_obj returned from objectPA is not in the expected range of -180 to 180 degrees, exiting'
        sys.exit()
        
    # due to symmetery we can simplify the problem
    if pa_obj < 0:
        pa_obj += 180
    if pa_mask > 180:
        # we do not want to redefine pa_mask since it is not really symmetric
        # due to the offcenter guider cam and other asymmetries
        pa_mask_prime = pa_mask-180
    else:
        pa_mask_prime = pa_mask
    
    # Determine the best allowable slit PA
    if pa_mask_prime >= pa_obj:
        if pa_mask_prime-pa_obj < relPA_min:
            pa_slit = pa_mask_prime-relPA_min
        elif pa_mask_prime-pa_obj < relPA_max:
            pa_slit = pa_obj
        else:
            pa_slit = pa_mask_prime-relPA_max
    elif pa_mask_prime < pa_obj:
        if pa_obj-pa_mask_prime < relPA_min:
            pa_slit = pa_mask_prime+relPA_min
        elif pa_obj-pa_mask_prime < relPA_max:
            pa_slit = pa_obj
        else:
            pa_slit = pa_mask_prime+relPA_max
    return pa_slit

def slitsize(pa_slit,sky,A_gal,B_gal=None,pa_gal=None):
    '''
    This function determines the slit size based on the size of the target and
    the requested sky on either side.  It returns desim input len1 and len2.
    
    Input:
    pa_slit = [1D list of floats; units:degrees] the position angle of the long
        axis of the slit, +CCW from North towards East
    sky = [1D float list; (len1 sky, len2 sky); units:arcsec] the amount of sky
        either side of the object size.
    A_gal = [1D list of floats; units:arcsec] galaxy radius along its major axis
    B_gal = [1D list of floats; units:arcsec] galaxy radius along its minor axis,
        if None then the galaxy is assumed circular with radius A_gal
    pa_gal = [1D list of floats; units:degrees] the position angle of the
        galaxy's major axis +CCW from North towards East
    '''
    from numpy import cos, sin, pi, sqrt
    deg2rad = pi/180.
    if B_gal == None and pa_gal == None:
        objrad = A_gal
    elif B_gal != None and pa_gal != None:
        objrad = sqrt((A_gal*cos(pa_slit-pa_gal))**2+
                       (B_gal*sin(pa_slit-pa_gal))**2)
    len1 = objrad+sky[0]
    len2 = objrad+sky[1]
    return len1, len2
    

def write_dsim_header(F,regfile,prefix):
    '''
    F = an opened file (e.g. F=open(filename,'w')
    box = 
    '''
    box = readMaskRegion(regfile)
    F.write('#This catalog was created by obsplan.py and is intended to be used \n')
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
            .format(prefix,box[0]/15.,box[1],box[4]))

def write_guide_stars(F,gs_ids,objid,ra,dec,magnitude,equinox='2000',passband='R'):
    '''
    F = dsim file
    gs_ids = (list of integers)
    
    '''
    N = numpy.size(gs_ids)
    # Possibly reformat radii array to enable single float input
    if N == 1:
        gs_ids = numpy.reshape(gs_ids,(1,))
    for i in gs_ids:    
        mask_s = objid == i
        if numpy.sum(mask_s) ==0:
            print 'obsplan.write_guide_stars: no objects in catalog match the guide star id:{0}, please check your input, exiting'.format(i)
            sys.exit()
        ra_i = tools.deg2ra(ra[mask_s],':')
        dec_i = tools.deg2dec(dec[mask_s],':')
        mag_i = magnitude[mask_s][0]
        F.write('{0}  {1}  {2}  {3}  {4:0.2f}  {5}  -1  0  1\n'
                .format(i,ra_i,dec_i,equinox,mag_i,passband))    

def write_align_stars(F,as_ids,objid,ra,dec,magnitude,equinox='2000',passband='R'):
    N = numpy.size(as_ids)
    # Possibly reformat radii array to enable single float input
    if N == 1:
        as_ids = numpy.reshape(as_ids,(1,))
    for i in as_ids:
        mask_s = objid == i
        if numpy.sum(mask_s) ==0:
            print 'obsplan.write_align_stars: no objects in catalog match the alignment star id:{0}, please check your input, exiting'.format(i)
            sys.exit()        
        ra_i = tools.deg2ra(ra[mask_s],':')
        dec_i = tools.deg2dec(dec[mask_s],':')
        mag_i = magnitude[mask_s][0]
        F.write('{0}  {1}  {2}  {3}  {4:0.2f}  {5}  -2  0  1\n'
                .format(i,ra_i,dec_i,equinox,mag_i,passband)) 

def write_galaxies_to_dsim(F,objid,ra,dec,magnitude,priority_code,sample,selectflag,pa_slit,len1,len2,equinox='2000',passband='R'):
    '''
    
    '''
    from math import floor
    for i in numpy.arange(numpy.size(objid)):
        #convert deg RA to sexadec RA
        ra_i = ra[i]/15.0
        rah = floor(ra_i)
        res = (ra_i-rah)*60
        ram = floor(res)
        ras = (res-ram)*60.
        #convert deg dec to sexadec dec
        dec_i = dec[i]
        if dec_i<0:
            sign = -1.
            dec_i = abs(dec_i)
        else:
            sign = 1.
        decd = floor(dec_i)
        res = (dec_i-decd)*60.
        decm = floor(res)
        decs = (res-decm)*60.
        if numpy.size(pa_slit) == 1:
            pa_slit_i = pa_slit
        else:
            pa_slit_i = pa_slit[i]
        if sign==-1:
            F.write('{0:0.0f}\t{1:02.0f}:{2:02.0f}:{3:06.3f}\t-{4:02.0f}:{5:02.0f}:{6:06.3f}\t{7}\t{8:0.2f}\t{9}\t{10:0.0f}\t{11:0.0f}\t{12:0.0f}\t{13:0.2f}\t{14:0.1f}\t{15:0.1f}\n'
                    .format(objid[i],rah,ram,ras,decd,decm,decs,equinox,magnitude[i],passband,priority_code[i],sample[i],selectflag[i],pa_slit_i,len1[i],len2[i]))
        else:
            F.write('{0:0.0f}\t{1:02.0f}:{2:02.0f}:{3:06.3f}\t{4:02.0f}:{5:02.0f}:{6:06.3f}\t{7}\t{8:0.2f}\t{9}\t{10:0.0f}\t{11:0.0f}\t{12:0.0f}\t{13:0.2f}\t{14:0.1f}\t{15:0.1f}\n'
                    .format(objid[i],rah,ram,ras,decd,decm,decs,equinox,magnitude[i],passband,priority_code[i],sample[i],selectflag[i],pa_slit_i,len1[i],len2[i]))

def makeSlitmaskRegion(prefix,ra,dec,pa_slit,length,sample,width=1):
    '''
    create a region file that maps the suggested slit of each galaxy
    '''
    outputname = prefix+'_slits.reg'
    out = open(outputname,'w')
    out.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+'\n')
    out.write('fk5'+'\n')
    for i in numpy.arange(numpy.size(ra)):
        ra_i = ra[i]
        dec_i = dec[i]
        length_i = length[i]
        if numpy.size(pa_slit) == 1:
            pa_slit_i = pa_slit
        else:
            pa_slit_i = pa_slit[i]        
        if sample[i] == 1:
            color = 'green'
        else:
            color = 'blue'
        out.write('box({0:1.5f},{1:1.5f},{2:1.1f}",{3:1.1f}",{4:0.0f}) # color={5}'.format(ra_i,dec_i,width,length_i,pa_slit_i,color)+'\n')
    out.close()