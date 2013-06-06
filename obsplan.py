'''
The purpose of this program is to help with planning the DEIMOS observations, in
particular the galaxy selection part of the mask making.

Note that while this program can handle multiple mask regions at one time,
the output deimos input will not be meaningful.

Version 0 of this code is based on obsplan.py that was used in the planning of
the Musket Ball 2011A DEIMOS run.

INPUT:
catalog = [string] file name of a ttype index, white space delimited, catalog
        of objects in the area surounding the input slit mask region file.
regfile = [string] file name of a ds9 region file for the box region
        representing the deimos slit mask.
prefix = [string] the prefix to add to all output file names
R_bounds = [(float,float) units:(magnitude,magnitude)] the (min,max) R-band
        magnitudes to allow (inclusive).  
dlsqc_bounds = [(float,float)] the (min,max) dlsqc values to allow (inclusive).
        The purpose of this filter is to remove stars based on thier psf shape.
JK_Rz = [(float,float)] ("J-K cut","R-z cut") exclued all objects with
        J-K < "J-K cut" and R-z > "R-z cut".  The purpose of this filter is to
        remove stars based on their colors.
redshift_bounds = [(float,float)] (min,max) photo-z to allow (inclusive)
wghtlist = [((string,float,float,real),(string,float,float,real),etc.)] a list
        of lists that specify the priorit_code weighting for each galaxy. Format
        (property,minprop,maxprop,additive weight for galaxies within prop range)
        where property is the ttype name of the property
sky = [(float,float) units:(arcsec,arcsec)] length in arcsec (above,below) with
        which to add to aR/2 for determining L1 and L2 of the minimum slit length
exfile = [string] file name of a ttype indexed file that contains all the
        objid's that should be excluded from the deimos input file
OUTPUT:
prefix_maskcat.txt = deimos slitmask input text file containing all the suitable
        galaxies for a deimos slit.  Includes PA and priority_code, but does not
        contain alignment stars
prefix_zgaldist.png = figure of the suitable galaxies' RA-dec distribution,
        color coded based on their photo-z
prefix_zhist.png = figure of the suitable galaxies' redshift histogram
prefix_prioritycodegaldist.png = figure of the suitable galaxies' RA-dec
        distribution, color coded based on their priority_code
prefix_circles.reg = ds9 region file that circles the suitable galaxies
prefix_slits.reg = ds9 region file that makes a slit for the suitable galaxies
        according to the prefix_maskcat.txt specifications
'''
import numpy
import pylab
from math import floor
import tools #module written by Will

###########################################################################
## USER INPUT 
###########################################################################
catalog = 'catalog.txt'
regfile = 'Mask78.reg'
maskNumber = 8 #Note that stars are currently clugely incorporated, search for maskNumber below to see how
#Prefix for all output files
prefix = 'Mask8_rev4'
#List of preselected galaxies to exclude from mask (excludes matching ttypeX = obsid
exfile = 'exfile_rev4.txt'
#Preselection list
presel ='preselect_mask8_rev4.txt' #a string (e.g. 'preselect_mask3_rev3.txt') or None

## Define hard catalog cuts/limits
R_bounds = (0,23.5)
dlsqc_bounds = (5,None) #dlsqc values to allow
JK_Rz = (1.2,0.75) #exclued all objects with J-K < 1.2 and R-z > 0.75
redshift_bounds = (None,None)

##Specify the priority_code weighting
#(property,minprop,maxprop,additive weight for galaxies within prop range)
wghtlist = (('R',0,23.5,10),('z_b',0.46,0.60,20))
#wghtlist = (('R',0,23.5,5),('R',0,23.5,5))

# Define the sample cut if the object has property less than this value assign
# it to sample 1 otherwise sample 3
samplecut = ('R',22.6)
#The amount of sky on either side of the galaxy to include in slit (arcsec)
sky = (1,1)


############################################################################
## PROGRAM
############################################################################

#########################################
#### Standard Survey Galaxies
#########################################

#read in the catalog and header
cat = tools.readcatalog(catalog)
key = tools.readheader(catalog)

## hard catalog filters
#filter the catalog
print 'obsplan: filter out all galaxies that have available spec-z'
cat = tools.catfilter(None,0,cat,key['z_spec'],max_inc=False)
print 'obsplan: apply R filter'
cat = tools.catfilter(R_bounds[0],R_bounds[1],cat,key['R'])
print 'obsplan: apply DLSQC filter'
cat = tools.catfilter(dlsqc_bounds[0],dlsqc_bounds[1],cat,key['dlsqc'])
print 'obsplan: apply the J-K R-z color cut'
#keep all galaxies that have a null magnitude
mask_J = cat[:,key['J']] < -10 
mask_K = cat[:,key['K']] < -10
mask_z = cat[:,key['z']] < -10
mask_JK = cat[:,key['J']]-cat[:,key['K']] > JK_Rz[0]
mask_Rz = cat[:,key['R']]-cat[:,key['z']] < JK_Rz[1]
mask_color = mask_J+mask_K+mask_z+mask_JK+mask_Rz > 0
Nint = numpy.shape(cat)[0]
cat = cat[mask_color,:]
Nfin = numpy.shape(cat)[0]
Ncut = Nint-Nfin
print 'obsplan: {0} rows were removed from the catalog with {1} initial rows, leaving {2} rows'.format(Ncut,Nint,Nfin)
print 'obsplan: apply the redshift filter'
cat = tools.catfilter(redshift_bounds[0],redshift_bounds[1],cat,key['z_b'])

## region file catalog filter
print 'obsplan: apply the regions filter'
# find all the box regions
box = numpy.fromregex(regfile,r"box\(([0-9]*\.?[0-9]+),([0-9]*\.?[0-9]+),([0-9]*\.?[0-9]+)\",([0-9]*\.?[0-9]+)\",([0-9]*\.?[0-9]+)",
                      [('xc',numpy.float),('yc',numpy.float),('width',numpy.float),('height',numpy.float),('angle',numpy.float)])

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
    ra_p = (cat[:,key['ra']]-xc)*numpy.cos(yc*d2r)*numpy.cos(-phi)+(cat[:,key['dec']]-yc)*numpy.sin(-phi)
    dec_p = -(cat[:,key['ra']]-xc)*numpy.cos(yc*d2r)*numpy.sin(-phi)+(cat[:,key['dec']]-yc)*numpy.cos(-phi)
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
Nint = numpy.shape(cat)[0]
cat = cat[mask,:]
Nfin = numpy.shape(cat)[0]
Ncut = Nint - Nfin
print 'obsplan: {0} rows were removed from the catalog with {1} initial rows, leaving {2} rows'.format(Ncut,Nint,Nfin)

####################################################
## Remove exclusion file objects
####################################################

#If an exclusion list was input then further filter the catalog
if exfile != None:
    print 'obsplan: apply exclusion list to further filter catalog'
    exkey = tools.readheader(exfile)
    exlist = numpy.loadtxt(exfile,usecols=(exkey['objid'],exkey['objid']))
    mask_ex = numpy.zeros(numpy.shape(cat)[0])
    i = 0
    for oid in exlist[:,0]:
        if i == 0:
            mask_ex = cat[:,key['objid']] == oid
            i=1
        else:
            mask_i = cat[:,key['objid']] == oid
            mask_tmp = mask_ex+mask_i
            mask_ex = mask_tmp > 0
    mask_ex = mask_ex == False
    Nint = numpy.shape(cat)[0]
    cat = cat[mask_ex,:]
    Nfin = numpy.shape(cat)[0]
    Ncut = Nint - Nfin
    print 'obsplan: {0} rows were removed from the catalog with {1} initial rows, leaving {2} rows'.format(Ncut,Nint,Nfin)

####################################################
## Determine slit parameters for each object
####################################################

#Determine the optimal PA for the slits
def PAround(PAarray,PAmin,PAmax,PAvalue,maskPA):
    '''
    Inspects each element of the PAarray and if it falls between PAmin and PAmax
    then it redefines the PA to the PAround value.
    Where the min and max bounds are in the masks coordinate system.
    '''
    maskPAmin = PAmin+maskPA
    maskPAmax = PAmax+maskPA
    if maskPAmax > 90:
        test1a = PAarray >= maskPAmin
        test1b = PAarray <= maskPAmax
        test1 = test1a*test1b
        test2a = PAarray >= maskPAmin-180
        test2b = PAarray <= maskPAmax-180
        test2 = test2a*test2b
        test = test1+test2
    elif maskPAmin <90:
        test1a = PAarray >= maskPAmin
        test1b = PAarray <= maskPAmax
        test1 = test1a*test1b
        test2a = PAarray >= maskPAmin+180
        test2b = PAarray <= maskPAmax+180
        test2 = test2a*test2b
        test = test1+test2
    else:
        test1a = PAarray >= maskPAmin
        test1b = PAarray <= maskPAmax
        test = test1a*test1b
    PAarray[test] = PAvalue+maskPA
    return PAarray
maskPA = box[0][4]-90
if maskPA > 90:
    maskPA -= 180
elif maskPA < -90:
    maskPA += 180

### Attempt to align the slits with the major axis of each galaxy
##PA = cat[:,key['theataR']]*1.0
# Attempt to align the slits with the minor axis of each galaxy
PA = cat[:,key['theataR']]*1.0+90
#The 0-5 spec is based on DEEP2 recommendations for better sky subtraction
PA = PAround(PA,0,5,5,maskPA)
PA = PAround(PA,-5,0,-5,maskPA)
PA = PAround(PA,30,90,30,maskPA)
PA = PAround(PA,-90,-30,-30,maskPA)

################################################################
## Determine the "priority_code" value for each object
################################################################

# calculate the priority_code values to assign to each galaxy
wght = numpy.zeros(numpy.shape(cat)[0])
for i in numpy.arange(numpy.shape(wghtlist)[0]):
    mask_wght_gt = cat[:,key[wghtlist[i][0]]]>=wghtlist[i][1]
    mask_wght_lt = cat[:,key[wghtlist[i][0]]]<=wghtlist[i][2]
    mask_wght = mask_wght_gt*mask_wght_lt
    wght += mask_wght*wghtlist[i][3]

####################################################
## Determine the "sample" assignment for each object
####################################################

# Determine the sample assignments
sample = numpy.ones(numpy.shape(cat)[0])*3
test = cat[:,key[samplecut[0]]] <= samplecut[1]
sample[test] = 1

####################################################
## Determine if object is preselected
####################################################

#preselect objects
pscode = numpy.zeros(numpy.shape(cat)[0])
if presel != None:
    #read in the preselection catalog and header
    pskey = tools.readheader(presel)
    print 'obsplan: determining preselections'
    pslist = numpy.loadtxt(presel,usecols=(pskey['objid'],pskey['objid']))
    i = 0
    for oid in pslist[:,0]:
        if i == 0:
            pscode = cat[:,key['objid']] == oid
            i=1
        else:
            mask_i = cat[:,key['objid']] == oid
            pscode += mask_i
            i += 1
    print 'obsplan: {0} slits preselected'.format(numpy.sum(pscode))



##### Stars ###################################################################
## This was the start at an attempt to generalize the guide and alignment star
## selection for the dsim input.  Should incorporate this in the future.
##if starcatalog != None:
##    #read in the catalog and header
##    starcat = tools.readcatalog(starcatalog)
##    starkey = tools.readheader(starcatalog)
##
##    print 'obsplan: apply the regions filter to the star catalog'
##    # find all the box regions
##    starbox = numpy.fromregex(starregfile,r"box\(([0-9]*\.?[0-9]+),([0-9]*\.?[0-9]+),([0-9]*\.?[0-9]+)\",([0-9]*\.?[0-9]+)\",([0-9]*\.?[0-9]+)",
##                          [('xc',numpy.float),('yc',numpy.float),('width',numpy.float),('height',numpy.float),('angle',numpy.float)])
##
##    d2r = numpy.pi/180.0
##    #loop through the regions creating masks for galaxy inclusion
##    for i in numpy.arange(numpy.shape(starbox)[0]):
##        #phi is the ccw angle from the +East axis
##        xc = starbox[i][0]
##        yc = starbox[i][1]
##        w = starbox[i][2]
##        h = starbox[i][3]
##        phi= starbox[i][4]*d2r
##        #rotate the galaxies into the "primed" (p) region coorditate frame centered
##        #at the center of the region
##        ra_p = (starcat[:,key['ra']]-xc)*numpy.cos(yc*d2r)*numpy.cos(-phi)+(starcat[:,key['dec']]-yc)*numpy.sin(-phi)
##        dec_p = -(starcat[:,key['ra']]-xc)*numpy.cos(yc*d2r)*numpy.sin(-phi)+(starcat[:,key['dec']]-yc)*numpy.cos(-phi)
##        #determine the min and max bounds of the region
##        # min = (box center [deg])-((box height [sec])/(2*60**2))
##        ra_p_min = -w/(2*60**2)
##        ra_p_max = w/(2*60**2)
##        dec_p_min = -h/(2*60**2)
##        dec_p_max = h/(2*60**2)
##        #create the mask for the i region
##        mask_ramin = ra_p >= ra_p_min
##        mask_ramax = ra_p <= ra_p_max
##        mask_decmin = dec_p >= dec_p_min
##        mask_decmax = dec_p <= dec_p_max
##        mask_i = mask_ramin*mask_ramax*mask_decmin*mask_decmax
##        #combine the i mask with the previous masks
##        if i == 0:
##            mask = mask_i
##        else:
##            #if the galaxy was in any mask then it should be in the concatenated mask
##            mask_tmp = mask + mask_i
##            mask = mask_tmp > 0
##
##    #apply the region filter to the input galaxy catalog
##    # see http://www.ucolick.org/~phillips/deimos_ref/masks.html for file format
##    Nint = numpy.shape(starcat)[0]
##    starcat = starcat[mask,:]
##    Nfin = numpy.shape(starcat)[0]
##    Ncut = Nint - Nfin
##    print 'obsplan: {0} rows were removed from the star catalog with {1} initial rows, leaving {2} rows'.format(Ncut,Nint,Nfin)


####################################################
## Create the dsimulator input file
####################################################

outcatname = prefix+'_maskcat.txt'
F = open(outcatname,'w')
F.write('#This catalog was created by obsplan.py and is intended to be used \n')
F.write('#as input to the deimos slitmask software following the format \n')
F.write('#outlined at http://www.ucolick.org/~phillips/deimos_ref/masks.html\n')
F.write('#Note that the automatic generation of this file does not include\n')
F.write('#guide or alignment stars.\n')
F.write('#ttype1 = objid\n')
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
##########
### Write in the Guide Star and Alignment Star lines
##########
# This is a very clugey way of incorporating stars into the dsim input that is
# difinitly not fit for general application.  I should eventually rewrite this
# aspect of the program.
if maskNumber == 1 or maskNumber == 2 or regfile == 'Mask12.reg':
    F.write('62224587489339  09:16:20.172   29:50:38.14   2000   16.41402   R  -1  0 1\n')
    F.write('62224587554914  09:16:36.821   29:56:18.59   2000   16.79682   R  -2  0 1\n')
    F.write('64092903178358  09:16:01.442   29:48:27.47   2000   15.82015   R  -2  0 1\n')
    F.write('64092903178264  09:16:05.707   29:50:16.60   2000   15.64186   R  -2  0 1\n')
    F.write('62224587423875  09:15:49.858   29:48:04.29   2000   16.66378   R  -2  0 1\n')
    F.write('64092903178317  09:16:15.058   29:56:02.61   2000   15.59187   R  -2  0 1\n')
elif maskNumber == 3 or maskNumber == 4 or regfile == 'Mask34.reg':
    F.write('64092903178307   09:16:28.651   29:48:46.64   2000   16.8204    R  -1  0 1\n')
    F.write('64092903178483   09:16:36.151   29:50:45.86   2000   16.42034   R  -1  0 1\n')
    F.write('64092903243910   09:16:36.821   29:56:18.59   2000   16.79812   R  -2  0 1\n')
    F.write('62224587489383   09:16:02.234   29:47:48.42   2000   16.4677    R  -2  0 1\n')
    F.write('64092903178358   09:16:01.442   29:48:27.47   2000   15.82015   R  -2  0 1\n')
    F.write('64092903178264   09:16:05.707   29:50:16.60   2000   15.64186   R  -2  0 1\n')
elif maskNumber == 5 or maskNumber == 6 or regfile == 'Mask56.reg':
    F.write('62224587489311   09:16:05.707   29:50:16.60   2000   15.64775   R  -1  0 1\n')
    F.write('62224587489399   09:16:06.646   29:50:41.63   2000   15.94431   R  -1  0 1\n')
    F.write('64092903178307   09:16:28.651   29:48:46.64   2000   16.8204    R  -2  0 1\n')
    F.write('64092903178424   09:16:03.026   29:56:55.92   2000   16.64095   R  -2  0 1\n')
    F.write('62224587489296   09:15:45.122   29:56:31.56   2000   16.10877   R  -2  0 1\n')
    F.write('64092903178272   09:16:01.913   29:55:02.29   2000   15.21739   R  -2  0 1\n')
    F.write('64092903178483   09:16:36.151   29:50:45.86   2000   16.42034   R  -2  0 1\n')
elif maskNumber == 7 or maskNumber == 8 or regfile == 'Mask78.reg':
    F.write('62224587489382   09:16:02.234   29:47:48.42   2000   16.46766   R  -1  0 1\n')
    F.write('64092903178358   09:16:01.442   29:48:27.47   2000   15.82015   R  -1  0 1\n')
    F.write('64092903178307   09:16:28.651   29:48:46.64   2000   16.8204    R  -2  0 1\n')
    F.write('64092903178260   09:15:56.345   29:52:13.05   2000   16.12409   R  -2  0 1\n')
    F.write('62224587489296   09:15:45.122   29:56:31.56   2000   16.10877   R  -2  0 1\n')
    F.write('64092903178272   09:16:01.913   29:55:02.29   2000   15.21739   R  -2  0 1\n')
    F.write('64667891859546   09:16:22.668   29:47:53.46   2000   16.7627    R  -2  0 1\n')
##########
### Write in the survey galaxies
##########
    
for i in numpy.arange(numpy.shape(cat)[0]):
    #convert deg RA to sexadec RA
    ra = cat[i,key['ra']]/15.0
    rah = floor(ra)
    res = (ra-rah)*60
    ram = floor(res)
    ras = (res-ram)*60.
    #convert deg dec to sexadec dec
    dec = cat[i,key['dec']]
    if dec<0:
        sign = -1.
        dec = abs(dec)
    else:
        sign = 1.
    decd = floor(dec)
    res = (dec-decd)*60.
    decm = floor(res)
    decs = (res-decm)*60.
## If attempting to align slit pa with major axis of each galaxy
    #if sign==-1:
        #F.write('{0:0.0f}\t{1:02.0f}:{2:02.0f}:{3:06.3f}\t-{4:02.0f}:{5:02.0f}:{6:06.3f}\t2000\t{7:0.2f}\tR\t{8:0.0f}\t{9:0.0f}\t0\t{10:0.2f}\t{11:0.1f}\t{12:0.1f}\n'
                #.format(cat[i,key['objid']],rah,ram,ras,decd,decm,decs,cat[i,key['R']],wght[i],sample[i],PA[i],cat[i,key['aR']]/2.+sky[0],cat[i,key['aR']]/2.+sky[1]))
    #else:
        #F.write('{0:0.0f}\t{1:02.0f}:{2:02.0f}:{3:06.3f}\t{4:02.0f}:{5:02.0f}:{6:06.3f}\t2000\t{7:0.2f}\tR\t{8:0.0f}\t{9:0.0f}\t0\t{10:0.2f}\t{11:0.1f}\t{12:0.1f}\n'
                #.format(cat[i,key['objid']],rah,ram,ras,decd,decm,decs,cat[i,key['R']],wght[i],sample[i],PA[i],cat[i,key['aR']]/2.+sky[0],cat[i,key['aR']]/2.+sky[1]))
# If attempting to align slit pa with minor axis of each galaxy
    if sign==-1:
        F.write('{0:0.0f}\t{1:02.0f}:{2:02.0f}:{3:06.3f}\t-{4:02.0f}:{5:02.0f}:{6:06.3f}\t2000\t{7:0.2f}\tR\t{8:0.0f}\t{9:0.0f}\t{10:0.0f}\t{11:0.2f}\t{12:0.1f}\t{13:0.1f}\n'
                .format(cat[i,key['objid']],rah,ram,ras,decd,decm,decs,cat[i,key['R']],wght[i],sample[i],pscode[i]*1,PA[i],cat[i,key['bR']]/2.+sky[0],cat[i,key['bR']]/2.+sky[1]))
    else:
        F.write('{0:0.0f}\t{1:02.0f}:{2:02.0f}:{3:06.3f}\t{4:02.0f}:{5:02.0f}:{6:06.3f}\t2000\t{7:0.2f}\tR\t{8:0.0f}\t{9:0.0f}\t{10:0.0f}\t{11:0.2f}\t{12:0.1f}\t{13:0.1f}\n'
                .format(cat[i,key['objid']],rah,ram,ras,decd,decm,decs,cat[i,key['R']],wght[i],sample[i],pscode[i]*1,PA[i],cat[i,key['bR']]/2.+sky[0],cat[i,key['bR']]/2.+sky[1]))

F.close()

####################################################
## Create various output plots
####################################################

#Plot the RA-dec distribution color coded based on photo-z
fig1 = pylab.figure()
pylab.scatter(cat[:,key['ra']],cat[:,key['dec']],s=20,
              c=cat[:,key['z_b']],cmap='spectral',alpha=0.75,edgecolors='none')
pylab.title('{0} Galaxies in Slit Mask Regions'.format(numpy.shape(cat)[0]))
pylab.xlabel('Right Ascension')
pylab.ylabel('Declination')
#Invert the RA-axis
pylab.xlim((pylab.xlim()[1],pylab.xlim()[0]))
cb = pylab.colorbar()
cb.set_label('photo-z')
figname = prefix+'_zgaldist'
#pylab.savefig(figname)

#Plot the redshift distribution of the selected galaxies
fig2 = pylab.figure()
#histgal = numpy.histogram(cat[:,key['z_b']],bins=40,range=(0,4))
#pylab.step(histgal[1][:-1],histgal[0],where='post',linewidth=1)
pylab.hist(cat[:,key['z_b']],bins=40,range=(0,4),log=True,histtype='bar',facecolor='g')
pylab.title('Redshift Histogram of Selected Galaxies')
pylab.xlabel('Photo-z')
pylab.ylabel('$N_{galaxy}$')
figname = prefix+'_zhist'
#pylab.savefig(figname)

#Plot the priority_code spatial distribution for this mask
fig3 = pylab.figure()
pylab.scatter(cat[:,key['ra']],cat[:,key['dec']],s=20,
              c=wght,cmap='Blues',alpha=0.75,edgecolors='none')
pylab.title('{0} Galaxies in Slit Mask Regions'.format(numpy.shape(cat)[0]))
pylab.xlabel('Right Ascension')
pylab.ylabel('Declination')
#Invert the RA-axis
pylab.xlim((pylab.xlim()[1],pylab.xlim()[0]))
cb = pylab.colorbar()
cb.set_label('priority_code')
figname = prefix+'_prioritycodegaldist'
#pylab.savefig(figname)

####################################################
## Create ds9 region files
####################################################

#create a region file that circles the selected galaxies
outputname = prefix+'_circles.reg'
F = open(outputname,'w')
F.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+'\n')
F.write('fk5'+'\n')
for i in numpy.arange(numpy.shape(cat)[0]):
    ra = cat[i,key['ra']]
    dec = cat[i,key['dec']]
    size = cat[i,key['aR']]+1
    obj = cat[i,key['objid']]
    F.write('circle({0:1.5f},{1:1.5f},{2:1.1f}") # text={{'.format(ra,dec,size)+'{0:0.0f}'.format(obj)+'}\n')
F.close()

#create a region file that maps the suggested slit of each galaxy
outputname = prefix+'_slits.reg'
F = open(outputname,'w')
F.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+'\n')
F.write('fk5'+'\n')
###
### Survey Galaxies
###
for i in numpy.arange(numpy.shape(cat)[0]):
    ra = cat[i,key['ra']]
    dec = cat[i,key['dec']]
    height = cat[i,key['aR']]+sky[0]+sky[1]
    angle = PA[i]+90
    if sample[i] == 1:
        color = 'green'
    else:
        color = 'blue'
    F.write('box({0:1.5f},{1:1.5f},{2:1.1f}",1",{3:0.2}) # color={4}'.format(ra,dec,height,angle,color)+'\n')
    
F.close()

   
#pylab.show()
# I am just playing with the following for now.

def optimalPA(h,delta,phi):
    '''
    This is based on:
    Filippenko, A.V., 1982. The importance of atmospheric differential refraction in spectrophotometry. Publications of the Astronomical Society of the Pacific, 94, pp.715–721. Available at: http://adsabs.harvard.edu/abs/1982PASP...94..715F.
    
    Input:
    phi = [float; units=degrees] observers latitude
    h = [float; units=hours] object's hour angle (h is + if west of the meridian)
    delta = [float; units=degrees] object's declination
    '''
    from math import pi,sin,cos,asin
    d2r = pi/180.
    phi *= d2r
    if h < 0:
        sign = -1
        h = -h
    else:
        sign = 1
    h = h*15*d2r
    delta *= d2r
    eta_rad = sign*asin(sin(h)*cos(phi) / 
                        (1-(sin(phi)*sin(delta) +
                            cos(phi)*cos(delta)*cos(h))**2)**(0.5))
    eta_deg = eta_rad/d2r
    return eta_deg