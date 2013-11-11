DEIMOS
======

Programs useful for Keck DEIMOS planning, observing, and reduction

This branch contains:
* an example ipython notebook that interacts with an existng ds9 instance called 'plan_zwcl2341_m3.ipynb' 
* 'obsplan_functions.py' which contains the function that the interactive ipython notebook calls 
* a tar.gz file with all the required input files for the IPython notebook 

Dependencies: 
* pyds9
* pyfits
* pandas 

Given a space separated catalog of galaxies+stars,  the ipython notebook script (with a running instance of ds9) is capable of: 
* weighing how likely each galaxies would be members based on photoz 
* making a number density map of the galaxies(potential members) in fits format 
* make ds9 contour of the number density map of galaxies 
* make ds9-regions that has appropriate angular sizes to mimick the slitmask
and the guider camera 
* filter out stars and display them as ds9-circle-regions on top of the slitmask for easy selection 
* read back in chosen stars and output them appropriately in dsim format
* write dsim input file 

etc. 

To do
---
* make sure that the functions are flexible for different types of catalogs, i.e. not just SDSS catalogs 
so users do not have to modify the code 
* better documentations of the functions 
* move Google Doc "How to create slitmask" to this Github as a wiki page
