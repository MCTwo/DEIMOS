'''
This file provides an example of how to interact and run obsplan.py

All position angles are defined +CCW from North towards East
'''
from __future__ import division
import tools
import obsplan

###########################################################################
## USER INPUT
###########################################################################

## General inputs

#Prefix for all output files
prefix = 'Mask1_rev0'
# Hour angle of the target for the mask (float; unit:hours)
HA = -40/60.

## Star and Galaxy catalog inputs

catalog = 'catalog.txt'
objid_ttype = 'objID'
ra_ttype = 'ra'
dec_ttype = 'dec'
mag_ttype = 'dered_r'
passband = 'R'
equinox = '2000'

# Go ahead and read in the catalog since this will be needed to create the
# galaxy mask later in the user input section
cat = tools.readcatalog(catalog)
key = tools.readheader(catalog)

## Slitmask ds9 region input

# Path/Name of the ds9 region file defining the bound and orientation of the
# slit mask. Note region should be defined using Coordinate/WCS/Degrees and 
# Size/WCS/Arcmin options, with the Size 5 by 16.1 arcmin, Angle will then
# correspond to the slitmask's parallactic angle (i.e. +CCW from north towards
# east) with the guider camera in the North-east quadrent at Angle=0.
regfile = 'Mask1_rev0.reg'

## Slit size inputs

#The amount of sky on either side of the galaxy to include in slit (arcsec)
sky = (1,1)
# The ttype index of the galaxy size. If one value is entered then the galaxy
# will be assumed circular, of three values are entered the the galaxy will be
# assumed elliptical with major axis radius (a), and minor axis radius (b), and
# position angle (pa_gal) of the major axis measured +CCW from North towards
# east
A_gal_ttype = 'deVRad_r'
B_gal_ttype = None # if None then galaxy assumed circular
pa_ga_ttype = None # if None then galaxy assumed circular

## Guide and alignment star inputs

# Enter lists of the object ids for each type of star, these will be matched to
# the object ids in the catalog

# Guide star id's
gs_ids = (440729287)
# Alignment star id's
as_ids = (646957174, 646957223, 646826017, 646826203, 646957144)

## Exclusion list input

# ttype catalog of galaxies to exclude from mask (excludes matching ttype = objid). exobjid_ttype is ['string'] ttype name of the objid column in the exfile, the objid's should correspond to some of the objid's in the objid array
exfile = None #a string (e.g. 'exclusion.txt') or None
exobjid_ttype = 'objid' 

## Preselected list input

# ttype catalog of galaxies to preselect in dsim.
psfile = None #a string (e.g. 'preselect.txt') or None
psobjid_ttype = None

## Create galaxy selection mask

# Some how the user at this point needs to create a boolean type mask for the
# catalog that will filter out all objects that are not galaxies
mask_galaxy = cat[:,key['type']] == 3

## Create a magnitude mask

# It is likely that the a faint end mask should be used
mask_mag = cat[:,key[mag_ttype]] <= 22.5

## Create a sample definition


###########################################################################
## Automated Portion
###########################################################################

# create basic 1D arrays from catalog
objid = cat[:,key[obid_ttype]]
ra = cat[:,key[ra_ttype]]
dec = cat[:,key[dec_ttype]]
mag = cat[:,key[mag_ttype]]

# Create the slitmask mask
mask_slitmask = obsplan.createSlitmaskMask(regfile,ra,dec)

# Create the exclusion list mask
mask_ex = obsplan.createExclusionMask(objid,exfile,exobjid_ttype)

# Determine the priority_code (i.e. weight) for each galaxy

# Determine the sample for each galaxy (sample 1 objects selected first, then
# sample 2, etc.). This order of selection take priority over the priority_code

# Determine the selection flag for each galaxy. If non-zero then the object is
# preselected
selectflag = assignSelectionFlag(objid,psfile,psobjid_ttype)

# determine object declination and the mask PA from the regfile
box = obsplan.readMaskRegion(regfile)
delta = box[1]
pa_mask = box[4]

# Determine the optimal slit PA
pa_slit = obsplan.optimalPA(pa_mask,HA,delta)

# Determine the slit size for each object
len1, len2 = obsplan.slitsize(pa_slit,sky,A_gal,B_gal,pa_gal)

# Create the output dsim file
outcatname = prefix+'_maskcat.txt'
print 'started to write out to ', outcatname
F = open(outcatname,'w')

# Write the dsim header information to the output file
obsplan.write_dsim_header(F,regfile)

# Write the guide star info to the dsim output file
obsplan.write_guide_stars(F,gs_ids,objid,ra,dec,mag,equinox,passband)

# Write the alignment star info to the dsim output file
obsplan.write_guide_stars(F,as_ids,objid,ra,dec,mag,equinox,passband)

# Filter the galaxy catalog before creating dsim input
mask = mask_galaxy*mask_mag*mask_ex*mask_slitmask

# Write the galaxy info to the desim output file
write_galaxies_to_dsim(F,objid[mask],ra[mask],dec[mask],mag[mask],priority_code[mask],sample[mask],selectflag[mask],pa_slit[mask],len1[mask],len2[mask],equinox='2000',passband='R')