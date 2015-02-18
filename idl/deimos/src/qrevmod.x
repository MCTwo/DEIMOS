# QREVMOD: simplified model; calls in the prev. gen maps and grating pars, but
# works in reverse to provide equal mappings on the slitmask.

include	<math.h>
include	<math/gsurfit.h>
include	"deimos.h"
include	"instrument.h"


procedure	t_qrevmod()

char	amap[SZ_FNAME], bmap[SZ_FNAME]		# input mappings
real	xmm, ymm, wave			# X,Y in slitmask, wave
real	xpix, ypix			# pixel values on CCD
bool	grid				# generate a grid of points?
pointer	fda, fdb

double	sys[NPARAM]				# system parameters
double	a3[3,3]				# grating transform
# double	r[3]
# double	alpha, beta, gamma
# real	tanx, tany			# should be double; check mapping evals
real	xics, yics			# pixel values in ICS

pointer	asfx, asfy			# pointers to surface fits (amap)
pointer	bsfx, bsfy			# pointers to surface fits (bmap)

int	i, j
real	wave0, wavinc, pixinc, xmmref
real	dely, w		# temp iteration to produce central value

bool	clgetb()
int	clgeti()
real	clgetr()
pointer	open()

begin
	xpix = clgetr ("xpix")
	ypix = clgetr ("ypix")
	wave = 1.e-4 * clgetr ("wave")			# in microns

	MU(sys) = DEGTORAD (clgetr ("mu"))
	GR_YERR(sys) = DEGTORAD (clgetr ("roll3"))
	GR_ZERR(sys) = DEGTORAD (clgetr ("o3"))

	ORDER(sys) = clgeti ("norder")
	GRLINES(sys) = 1.e-3 * clgeti ("gmm")		# in microns

	call clgstr ("amap", amap, SZ_FNAME)
	fda = open (amap, READ_ONLY, TEXT_FILE)

	call clgstr ("bmap", bmap, SZ_FNAME)
	fdb = open (bmap, READ_ONLY, TEXT_FILE)

# Iitialize the maps
	call gs_ingest (fda, asfx, asfy)
	call gs_ingest (fdb, bsfx, bsfy)

#######################################################################
# Below is the actual calculation; gseval simply evaluates the mappings
#######################################################################

# set up the grating transform
	call gsetup (a3, sys)

# get mapping into ICS pixels
# FIX!  BUT eventually want to use ICS
	xics = xpix
	yics = ypix


	call qrevmod (xics, yics, wave, a3, sys, asfx, asfy, bsfx, bsfy, xmm, ymm)


call printf ("%8.1f %7.1f %8.2f -->ICS: %7.1f %7.1f --> %8.3f %7.3f \n")
	call pargr (xpix)
	call pargr (ypix)
	call pargr (wave*1.e4)
	call pargr (xics)
	call pargr (yics)
	call pargr (xmm)
	call pargr (ymm)

	grid = clgetb ("grid")
	if (grid) {
		yics = 0.
		wave0 = wave
		wavinc = 100.e-4	# in microns, thank you
		wavinc = 60.e-4
		pixinc = 4200. / 8.
		do i = -8, 8 {
		    xics = i * pixinc
		    w = wave0
# This is an attempt to move the reference to the 270 arcsec line
		    do j=1, 10 {
			call qrevmod (xics, yics, w, a3, sys, asfx, asfy, bsfx,
								bsfy, xmm, ymm)
			dely = ymm - 270.
			w = w - dely*1.e-4
#			call eprintf ("%6f\n")
#				call pargr (w)
		    }

			
		    xmmref = xmm
		    do j = -3, 7 {
#			wave = wave0 + j * wavinc
			wave = w + j * wavinc
			call qrevmod (xics, yics, wave, a3, sys, asfx, asfy,
							bsfx, bsfy, xmm, ymm)
#			if (ymm < 0. || ymm > 349.)
#				next
#			if (sqrt (xmm*xmm+(ymm+131.**2)) > 436.6)
#				next
			if (ymm < 180. || ymm > 480.)
				next
			if (sqrt (xmm*xmm+ymm*ymm) > 600.)
				next
			call printf ("%8.3f %8.3f  %6.3f  %7.3f\n")
				call pargr (xmm)
				call pargr (ymm)
				call pargr (xmm-xmmref)
				call pargr (xmmref)
		    }
		}
	}
			
end


procedure	qrevmod (xics, yics, wave, a3, sys, asfx, asfy, bsfx, bsfy, xmm, ymm)

real	xics, yics			# pixel values in ICS
real	wave				# wave in microns
double	sys[NPARAM]			# system parameters
double	a3[3,3]				# grating transform
pointer	asfx, asfy, bsfx, bsfy		# pointers to surface fits
real	xmm, ymm			# X,Y in slitmask

double	r[3]
double	alpha, beta, gamma
real	tanx, tany			# should be double; check mapping evals

real	gseval()
begin
# Get mapping and convert to r[3]
	tanx = gseval (bsfx, xics, yics)
	tany = gseval (bsfy, xics, yics)

	r[3] = -1.d0 / sqrt (1. + tanx*tanx + tany*tany)
	r[1] = r[3] * tanx
	r[2] = r[3] * tany

# xform into grating system
	call gen_xfm (r, a3, YES)

# convert to beta,gamma (signs may not be right)
	beta = -atan2 (-r[2], -r[3])
	gamma = atan2 (r[1], sqrt (r[3]*r[3]+r[2]*r[2]))

# Apply the grating equation
	alpha = asin ((ORDER(sys)*GRLINES(sys)*wave / cos (gamma)) - sin (beta))

call eprintf ("alpha, beta: %6f %6f\n")
call pargd (RADTODEG(alpha))
call pargd (RADTODEG(beta))

# convert alpha, gamma into x,y,z (cf Schroeder p259); note sign reversal of alpha
	r[1] = sin (gamma)
	r[2] = sin (-alpha) * cos (gamma)
	r[3] = cos (-alpha) * cos (gamma)

# xform out of grating system
	call gen_xfm (r, a3, NO)

# convert to tanx, tany
	tanx = (-r[1] / -r[3])
	tany = (-r[2] / -r[3])
call eprintf ("tanx,tany: %5f %5f\n")
call pargr (tanx)
call pargr (tany)

# get mapping into Slitmask coord
	xmm = gseval (asfx, tanx, tany)
	ymm = gseval (asfy, tanx, tany)
end
