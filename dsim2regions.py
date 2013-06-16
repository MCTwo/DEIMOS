#!/usr/bin/python

import pyfits
import sys
import numpy

# following dsim2regions.pro conventions
color = {'A':'blue','P':'cyan','G':'yellow'}
tvpolygonx = numpy.array([-1.0,208.0,208.0,-1.0,-1.0,-1.0,208.0])
tvpolygony = numpy.array([ 298.4,302.7,94.0,94.0,298.4,174.0,174.0])
maskpolygonx = numpy.array([-498.0,-498.0,360.0,420.0,460.0,498.0,498.0,-498.0,-498.0,-460.0,-420.0,-360.0,259.7,259.7,259.7,249.3,249.3,5.2,5.2,-5.2,-5.2,-249.3,-249.3,-259.7,-259.7])
maskpolygony = numpy.array([187.0,479.0,479.0,428.0,385.0,332.0,187.0,187.0,332.0,385.0,428.0,479.0,479.0,187.0,479.0,479.0,187.0,187.0,479.0,479.0,187.0,187.0,479.0,479.0,187.0])
arcsec_per_radian = 206264.8
deg_per_radian = 180.0/numpy.pi

def unflatten(dx,dy,ra_center,dec_center):
    # center of projection
    sindeccent = numpy.sin(dec_center/deg_per_radian)
    cosdeccent = numpy.cos(dec_center/deg_per_radian)

    # convert offset from arcsec to radians
    x = dx/arcsec_per_radian
    y = dy/arcsec_per_radian
    
    # do the spherical trig
    d = numpy.arctan(numpy.sqrt(x**2 + y**2))
    theta = numpy.arctan2(-x,y)
    dec = numpy.arcsin(numpy.cos(d)*sindeccent + numpy.sin(d)*cosdeccent*numpy.cos(theta))*deg_per_radian
    delta_ra = numpy.arcsin(numpy.sin(theta)*numpy.sin(d)/numpy.cos(dec))*deg_per_radian
    ra = ra_center + delta_ra
    return ra,dec

def plot_outline(ractr,decctr,maskpa,x,y,color):
    # adjust for mask ctr offset
    x = -x-35.0

    # rotate outline to mask PA
    u=x*numpy.cos(maskpa/deg_per_radian)-y*numpy.sin(maskpa/deg_per_radian)
    v=x*numpy.sin(maskpa/deg_per_radian)+y*numpy.cos(maskpa/deg_per_radian)
    # projection onto sky
    outlinera,outlinedec = unflatten(u,v,ractr,decctr)
    print "polygon(",
    for i in range(len(outlinera)-1):
        mystring = "%.6fd,%.6fd," % (outlinera[i],outlinedec[i])
        print mystring,
    # last point is special
    print "%.6fd,%.6fd) # color=%s" % (outlinera[i+1],outlinedec[i+1],color)
    return

### start main ####
usage = 'Usage: dsim2regions.py infile\n'
if len(sys.argv)!=2:
    sys.stderr.write(usage)
    sys.exit(1)

try:
    dsim = pyfits.open(sys.argv[1])
except:
    sys.stderr.write('could not open file\n')
    sys.exit(1)

# dsim[4] has the slits
ra = dsim[4].data.field('slitRA')
dec = dsim[4].data.field('slitDec')
pa = dsim[4].data.field('slitLPA')+90
slitwidth = dsim[4].data.field('slitWid')
slitlength = dsim[4].data.field('slitLen')
slittype = dsim[4].data.field('slitTyp')
# get ctr ra/dec and PA for the mask/tv regions
ractr = dsim[3].data.field('RA_PNT')[0]
decctr = dsim[3].data.field('DEC_PNT')[0]
maskpa = dsim[3].data.field('PA_PNT')[0] -90
sys.stderr.write("Mask ra/dec/pa: %s %s %s\n" % (ractr,decctr,maskpa))
dsim.close()

# print header
print "# slits from %s via dsim2regions.py" % (sys.argv[1])
print "global color=red\nj2000"

for i in xrange(len(ra)):
    print 'box(%f,%f,%.1f",%.1f",%d) # color=%s' % (ra[i],dec[i],slitlength[i],slitwidth[i],int(pa[i]),color[slittype[i]])

plot_outline(ractr,decctr,maskpa,tvpolygonx,tvpolygony,'yellow')
plot_outline(ractr,decctr,maskpa,maskpolygonx,maskpolygony,'red')

