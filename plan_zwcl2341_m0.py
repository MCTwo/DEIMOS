#===============================================================================
#Purpose: This makes use of "obsplan_functions.py" to make the masks 
#Author: 
#Will. A. Dawson <will@dawsonresearch.com>
#Karen Y. Ng <karen.yng@ucdavis.edu>
#License: BSD
#Date: 06/07/2013
#===============================================================================
# Import required libraries
#-------------------------------------------------------------------------------
from __future__ import division
import sys
import numpy as np,numpy
import pandas as pd 
import matplotlib.pyplot as plt
#import ds9
import numpy.ma as ma
import pyfits 
#from IPython.core.display import Image  
import aplpy
import obsplan_functions as fcn

# Initialize variables
#-------------------------------------------------------------------------------
catalog = 'zwcl2341_sdsscat.csv'
maskNumber = 1 #2 
prefix = 'zwcl2341_rev2_m0'
regfile1 = 'mask0_rev1.reg'
#regfile2 = 'mask1_rev0.reg'
R_bounds = (0,23.5)
#the amount of sky on either side of galaxy to include in slit (arcsec)
sky = (1.5, 1.5) 
#the known redshift for cluster ZwCL 2341 is 
z_cluster = 0.27
exclude_file = None

#main program
#-------------------------------------------------------------------------------
#By default "pandas.read_csv" reads in the files with ',' as delimiter
#na_values tells pandas to put that string expression as np.na values
#for example of how to read other table types see the end of the file
cat = pd.read_csv(catalog,na_values='null')
cat=cat.drop_duplicates(cols='objID')


#filter out other astronomical objects and leave only galaxies 
#type = 3 is for galaxies 
cat = fcn.filter_catalog_dataframe(cat,field='type',lowerbound=3,upperbound=3)

#filter out entries with spec_z 
cat = fcn.filter_catalog_dataframe(cat,field='z_spec',
                             lowerbound=-sys.float_info.epsilon,
                             upperbound=sys.float_info.epsilon)
#filter out entires with faint unextincted R band magnitude
cat = fcn.filter_catalog_dataframe(cat,field='dered_r',
                             lowerbound=R_bounds[0],
                             upperbound=R_bounds[1])

cat, box = fcn.return_objects_in_mask_region(cat,regfile1)

cat = fcn.determine_weight(cat,cat['z_phot'],cat['z_phot_Err'],
                                z_cluster)
cat = fcn.determine_sample_no(cat,2,'dered_r', first_sample=True,
                              lower_criteria=22.5,
                              upper_criteria=23.0)
cat = fcn.determine_sample_no(cat,1,'dered_r',lower_criteria=None,
                              upper_criteria=22.5)

#pscode is the preselection code, should be 0 for not being preselected
cat['pscode'] = pd.Series(np.zeros(cat.shape[0]),cat.index) 

# Exclude objects from exclude file list 
# no excluded objects for mask 1 
# use fcn.filter_catalog_dataframe(upper_critiera=cat['objID'], )

cat = fcn.convert_to_sexadec_coord(cat)
PA_field = 'PA'
#initialize field in the catalog for storing the PA values 
cat['PA']=pd.Series(np.zeros(cat.shape[0]),cat.index)
cat = fcn.pick_PA(cat, PA_field, box)

dsim_input = prefix+'_maskcat.txt'
F = open(dsim_input,'w') 
fcn.write_dsim_header(F,prefix,box)
#Guide stars with preselection code = -1 
#Had to write a function to crop first 3 digits of  the 19-digit SDSS ObjID to
#16 digits , Dsim can only accept 16 digit objID 
F.write("6408440594592	23:43:37.825	00:14:25.595	2000	18.15	R	-1	0	1\n")
F.write("6408440660224	23:43:44.062	00:14:22.538	2000	18.45	R	-1	0	1\n")
#Alignment stars with preselection code = -2 
F.write("3784195588160	23:43:48.919	00:11:51.639	2000	15.68	R   -2	0	1\n")
F.write("6408440594676	23:43:31.234	00:23:15.976	2000	19.40	R	-2	0	1\n")
F.write("6408440659974	23:43:42.240	00:23:23.067	2000	15.53	R	-2	0	1\n")
F.write("6408440660266	23:44:02.242	00:12:54.387	2000	16.08	R	-2	0	1\n")
F.write("6408440660187	23:43:39.072	00:23:43.484	2000	18.16	R	-2	0	1\n")

fcn.write_galaxies_to_dsim(F, cat, sky)
F.close()
print 'File for inputting to dsim with the name :'
print dsim_input
print 'has been written.'

fcn.write_slit_reg(cat,prefix,sky)


##write out objID for cross matching with the shortened ID 
#temp = []
#for i in cat.index:
#    temp.append('{0:04d}'.format(i))
#cat['shortID'] = pd.Series(temp,cat.index)
#
#cat.to_csv('m0_shortedID.txt',cols=('objID','shortID'),index=False)
#write 


#-----
#example code for reading in Will 's version of the catalog
#note 'sep=' option in the input denotes the regex expression 
#for recognizing separation between 
#different parameters e.g. 'ra', 'dec' (they are called fields) 
#e.g. cat = pd.read_table(catalog, na_values='null', skiprows = 56,
#                   sep=r"\s*",names=array_of_strings_for_column_names)
#-----
