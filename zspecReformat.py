import pyfits
import numpy

#user input
fitsfile = '/sandbox/deimos/zspec.dawson.1rxs1B.2013-02-18.fits'
outputfile = '/sandbox/deimos/1rxs1B/zspec.dawson.1rxs1B.2013-02-18.txt'

hdulist = pyfits.open(fitsfile)

tbdata = hdulist[1].data

#extract desired columns
objid = tbdata.field('OBJNAME')
slit = tbdata.field('SLITNAME')
mask = tbdata.field('MASKNAME')
z = tbdata.field('Z')
zerr = tbdata.field('Z_ERR')
qual = tbdata.field('ZQUALITY')
comm = tbdata.field('COMMENT')

#print the selected columns to a text file
f = open(outputfile,'w')
f.write('This data is extracted from the following zspec output:\n')
f.write(fitsfile+'\n')
f.write('ttype0=objid\n')
f.write('ttype1=slit\n')
f.write('ttype2=mask\n')
f.write('ttype3=z\n')
f.write('ttype4=zerr\n')
f.write('ttype5=quality\n')
f.write('ttype6=comment\n')
for i in numpy.arange(numpy.size(objid)):
    f.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(objid[i],slit[i],mask[i],z[i],zerr[i],qual[i],comm[i]))
    i+=1
f.close()