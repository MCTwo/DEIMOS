'''
The purpose of this program is to take an SDSS color finding chart, the type found at:
http://skyserver.sdss3.org/dr9/en/tools/chart/chart.asp ,
and provided the center (RA,Dec) and image pixel scale, it creates a color fits image
for use with programs like ds9, including WCS.

The original version of this of program was written by Dave Wittman.
'''
import numpy
import pyfits
import os
import sys
import re

#usage = 'Usage: jpgtofits.py jpegname fitsname ra dec arcsecperpix\n'
#if len(sys.argv) != 6:
#    sys.stderr.write(usage)
#    sys.exit(1)

#(inname,outname,ra,dec,scale) = sys.argv[1:]
#ra = float(ra)
#dec = float(dec)
#scale = float(scale)

inname = '/Users/dawson/SkyDrive/Observing/Keck2014a/A1612/A1612findingcart_ra191.93dec-2.79222scale0.79224.ppm'
outname = '/Users/dawson/SkyDrive/Observing/Keck2014a/A1612/sdss079224.fits'
ra = 191.93
dec = -2.79222
scale = 0.79224

#if os.system("jpegtopnm %s>tmp.ppm" % (inname)):
#    sys.stderr.write("jpegtopnm failed, exiting\n")
#    sys.exit(1)

# scan the header
nchars = 0
infile = open(inname)
line = infile.readline()
nchars+=len(line)
if line[:2] != 'P6':
    sys.stderr.write("ppm header failed, exiting\n")
    sys.exit(1)
line = infile.readline()
nchars+=len(line)
matches = re.match('(\d+) (\d+)',line)
xs = int(matches.group(1))
ys = int(matches.group(2))
line = infile.readline()
nchars+=len(line)
# done with hdr, now that we know the number of bytes to skip

# reshape it the way FITS wants it
data = numpy.fromfile(inname,dtype=numpy.uint8)
reddata = data[nchars::3].reshape([xs,ys])
grndata = data[nchars+1::3].reshape([xs,ys])
bluedata = data[nchars+2::3].reshape([xs,ys])
newdata = numpy.empty([3,xs,ys])
newdata[0,:,:] = reddata[::-1,:] # flip it in y
newdata[1,:,:] = grndata[::-1,:]
newdata[2,:,:] = bluedata[::-1,:]

# write the new file
newhdu=pyfits.PrimaryHDU(newdata)
newhdu.header.update('CTYPE1','RA---TAN')
newhdu.header.update('CTYPE2','DEC--TAN')
newhdu.header.update('CRPIX1',xs/2)
newhdu.header.update('CRPIX2',ys/2)
newhdu.header.update('CRVAL1',ra)
newhdu.header.update('CRVAL2',dec)
newhdu.header.update('CD1_1',-scale/3600)
newhdu.header.update('CD2_2',scale/3600)
newhdu.header.update('CD1_2',0.0)
newhdu.header.update('CD2_1',0.0)
newhdu.writeto(outname)
sys.exit(0)
