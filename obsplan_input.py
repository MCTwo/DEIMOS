'''
This file provides an example of how to interact and run obsplan.py

All position angles are defined +CCW from North towards East
'''

###########################################################################
## USER INPUT 
###########################################################################

## General inputs

#Prefix for all output files
prefix = 'Mask1_rev0'

## Star and Galaxy catalog inputs

catalog = 'catalog.txt'
ra_ttype = 'ra'
dec_ttype = 'dec'
mag_ttype = 'dered_r'

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

#List of preselected galaxies to exclude from mask (excludes matching ttypeX = obsid
exfile = 'exfile_rev4.txt'

## Preselected list input
#Preselection list
presel ='preselect_mask1_rev0.txt' #a string (e.g. 'preselect.txt') or None

