# TVALIGN:  Work program for estimating height of mirror above mask form

include	<math.h>
include	"deimos.h"
include	"keck.h"

define	D_1	20018.4D0	# see deimos.h PPLDIST
define	R_CURV	2124.71D0	# ditto; r=83.65in in this case = R_IMSURF


procedure	t_tvalign()

double	xmm, ymm

double	xp, yp, ap		# dummy variables
double	hm, hs

real	clgetr()

begin
	xmm = clgetr ("x")
	ymm = clgetr ("y")


	call mask_vs_sph (xmm, ymm, 0.d0, xp, yp, ap, hs, hm)

	call printf ("x, y, hs, hm, delh: %6.2f %6.2f %6.3f %6.3f %6.3f\n")
		call pargd (xmm)
		call pargd (ymm)
		call pargd (hs)
		call pargd (hm)
		call pargd (hs-hm)
end


# slightly hacked version of MASK_TO_PROJ

procedure mask_vs_sph (xc, yc, ac, xp, yp, ap, hs, hm)

double	xc, yc			# x,y values on mask surface
double	ac			# position angle on curved surface
double	xp, yp			# returned x,y in (focal)-plane system
double	ap			# returned position angle in planar system

double	mu, cosm			# mu, cos (mu)
double	cost, tant			# cos, tan of mask tilt angle
double	tanpa				# tan PA

double	rho			# radius from telescope optical axis
double	hs, hm			# height of image surface, mask above datum
double	xx, yy		# Work variables corresponding to xp, yp

begin
	mu = xc / M_RCURV
	cosm = cos (mu)
	cost = cos (DEGTORAD(M_ANGLE))
	tant = tan (DEGTORAD(M_ANGLE))
	xx =  M_RCURV * sin (mu)
	yy =  (yc - M_RCURV * tant * (1. - cosm)) * cost + ZPT_YM

	tanpa = (tan (DEGTORAD(ac)) - tant * xx / M_RCURV) * cost / cosm
	ap = atan (tanpa)

# What follows is a small correction for the fact that the mask does
# not lie exactly in the spherical image surface (where the distortion
# values and gnomonic projection are calculated) and the rays are arriving
# from the pupil image; thus, the exact locations are moved _slightly_
# wrt the telescope optical axis.  Note also that these corrections are
# only calculated to first order.

# Spherical image surface height:
	rho = sqrt (xx * xx + yy * yy)
	hs = R_IMSURF * (1. - sqrt (1. - (rho / R_IMSURF) ** 2))
# Mask surface height:
	hm = MASK_HT0 + yc * sin (DEGTORAD(M_ANGLE)) + M_RCURV * (1. - cosm)
	yp = yy - (hs - hm) * yy / PPLDIST 
	xp = xx - (hs - hm) * xx / PPLDIST 
end



