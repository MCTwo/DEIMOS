'''
This file provides an example of how to interact and run obsplan.py
'''
## Slitmask ds9 region input

# Path/Name of the ds9 region file defining the bound and orientation of the
# slit mask. Note region should be defined using Coordinate/WCS/Degrees and 
# Size/WCS/Arcmin options, with the Size 5 by 16.1 arcmin, Angle will then
# correspond to the slitmask's parallactic angle (i.e. +CCW from north towards
# east) with the guider camera in the North-east quadrent at Angle=0.
regfile = 'Mask1_rev0.reg'
