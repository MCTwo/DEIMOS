
########################
'''
This code was first written by Karen Ng <karenyng@ucdavis.edu> 
on 12/13/2012 by modifying from Will Dawson 's code

Purpose: for selecting guide & alignment stars of DEIMOS slitmasks.
For details of how to make the slitmasks,
see http://www.ucolick.org/~phillips/deimos_ref/masks.html#retrieve
'''


import numpy 
import pylab

from math import floor
from math import exp
from math import pi

import tools
import pdb

import make_masks_functions as fcn
import pandas as pd 
import pdb

###########################################################################
#USER INPUT
###########################################################################
catalog = 'zwcl2341_sdsscat.csv'

#this region file should represent your slitmask location 
regfile = 'mask0_rev0.reg'

#name to be added to all the output files of this program
prefix = 'stars_m0'




###########################################################################
#read in the SDSS star catalog 
###########################################################################
cat = pd.read_csv(catalog,na_values='null')
#R magnitude bounds
R_bounds = (15,19) 


#########################################
#filter the catalog
#########################################
print 'obsplan: filter out galaxies so that only stars remain'
#read sdss description of having 6 as stars , 3 as galaxies
#http://tdc-www.harvard.edu/catalogs/sdss.html
cat = fcn.filter_catalog_dataframe(cat,'type',6,6)
cat = fcn.filter_catalog_dataframe(cat,'dered_r',15,19.5)

################################################################
# find all the box regions
# ideally this box should represent your mask area 
################################################################
print 'obsplan: apply the regions filter'
box = numpy.fromregex(regfile,r"box\(([0-9]*\.?[0-9]+),([0-9]*\.?[0-9]+),([0-9]*\.?[0-9]+)\",([0-9]*\.?[0-9]+)\",([0-9]*\.?[0-9]+)",
                      [('xc',numpy.float),('yc',numpy.float),('width',numpy.float),('height',numpy.float),('angle',numpy.float)])

if len(box[0]) == 0: 
    print 'WARNING!!!!! there is something wrong with your box region file format~!\n You have put a regfile with wrong coordinate format'

d2r = numpy.pi/180.0
#loop through the regions creating masks for galaxy inclusion
for i in numpy.arange(numpy.shape(box)[0]):
    #phi is the ccw angle from the +East axis
    xc = box[i][0]
    yc = box[i][1]
    # make the box wider to mimic the area of the guide camera 
    # somehow the variables for the width and height are flipped  
    w = box[i][2]
    h = box[i][3]*1.5 
    phi=box[i][4]*d2r
    #rotate the galaxies into the "primed" (p) region coorditate frame centered
    #at the center of the region
    ra_p = (cat['ra']-xc)*numpy.cos(yc*d2r)*numpy.cos(-phi)+(cat['dec']-yc)*numpy.sin(-phi)
    dec_p = -(cat['ra']-xc)*numpy.cos(yc*d2r)*numpy.sin(-phi)+(cat['dec']-yc)*numpy.cos(-phi)
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
print 'obsplan: {0} rows were removed from the catalog with {1} initial rows, leaving {2} rows'.format(Ncut,Nint,Nfin)

#-------------
#process the disgusting SDSS 19-digit id so that there 's only 16-digits
#remaining for dsim input
#------------
cat['objID'] = cat['objID']-1230000000000000000
for i in cat.index:
    print '{0:19d}'.format(cat['objID'][i])

#######################################################
#create a file that writes out properties of the candidate stars
#######################################################
#one should copy and paste output of this file to obsplan.py
#for the section for writing out prefix_maskcat.txt
output = 'candidate_'+prefix+'.txt'
print 'Writing list of candidate stars to file "'+output+'"'
F=open(output,'w')
for i in cat.index:
    #Warning if your OBJID is longer than 18 digit 
    #You need to change your output line 
    obj = cat['objID'][i]
    ra = cat['ra'][i]
    ra = str(tools.deg2ra(ra,":"))
    dec = cat['dec'][i]
    dec = str(tools.deg2dec(dec,":"))
    #rmag = cat[i,key['dered_r']]
    #!!!!Warning!!!! if your OBJID is longer than 19 digit 
    #You need to change the digit precision of your output line 
    output= 'F.write("'+'{0:19d}\t'.format(obj)+ra+'\t'+dec+'\t2000\t{0:2.2f}'.format(cat['dered_r'][i])+'\tR\t-2\t0\t1' +r'\n")'+'\n'
    F.write(output)



#######################################################
#create a region file that circles the selected stars
#######################################################
outputname = prefix+'_circles.reg'
print 'Writing region of candidate stars to file "'+outputname+'"'

F = open(outputname,'w')
F.write('global color=yellow dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+'\n')
F.write('fk5'+'\n')
for i in cat.index:
    ra = cat['ra'][i]
    dec = cat['dec'][i]
    size = cat['deVRad_r'][i]+10
    #!!!!Warning!!!! if your OBJID is longer than 18 digit 
    #You need to change the digit precision of your output line 
    obj = cat['objID'][i]
    R = cat['dered_r'][i]
    #type = cat['type'][i]
    #!!!!Warning!!!! if your OBJID is longer than 19 digit 
    #You need to change the digit precision of your output line 
    F.write('circle({0:1.5f},{1:1.5f},{2:1.1f}")#text={{'.format(ra,dec,size)+'{0:19d}--{1:1.2f}'.format(obj,R)+'}\n')
F.close()


