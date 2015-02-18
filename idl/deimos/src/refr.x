# REFR:  Routines for refraction, plus soon there will be a parabolic collimator
# added,

include	<math.h>

define	D_1	20018.4D0	# see deimos.h PPLDIST
define	R_CURV	2144.40D0	# LRIS codeV

##
# The incoming ray is described as a reference pt (x,y,z) and an angle r[3].
# xs,etc and xc,etc are in the reference frame of the lens axis
# 

procedure	sph_refr (x, r, x0, y0, z0, rcurv, n1, n2)

double	x[3]			# reference point
double	r[3]			# ray angle
double	x0, y0, z0		# x,y,z location of lens cent curv
double	rcurv			# radius of curvature
double	n1, n2			# indices of refraction

double	a[3,3]			# xform matrix
double	xc, yc			# x,y location of ray at cent of curv (zc=0)
double	xs, ys, zs		# x,y,z at lens surface
double	q
double	dz
double	phi, theta		# angles at lens surface

double	phin, thetan
double	rho
double	cosp, sinp, cost, sint

begin
# translate rays into lens system; don't need to change angles, only (x,y,z):
	dz = x[3] - z0 - rcurv
	xc = x[1] - x0 - dz*r[1]/r[3]
	yc = x[2] - y0 - dz*r[2]/r[3]

# Find the location of the intercept of ray with lens surface
	q = r[1] * xc + r[2] * yc
	zs = r[3] * (-q + sqrt (q*q + rcurv*rcurv - (xc*xc + yc*yc)))
	xs = xc + r[1]/r[3] * zs
	ys = yc + r[2]/r[3] * zs

# get normal to lens surface
	rho = sqrt (xs*xs + ys*ys)
	theta = atan ( rho / zs)
	phi = atan2 (ys, xs)

# calc the transform matrix
	thetan = theta
	phin = phi + HALFPI
	cosp = cos (phin)
	sinp = sin (phin)
	cost = cos (thetan)
	sint = sin (thetan)

	a[1,1] = cosp
	a[1,2] = sinp
	a[1,3] = 0.
	a[2,1] = -cost*sinp
	a[2,2] = cost*cosp
	a[2,3] = sint
	a[3,1] = sint*sinp
	a[3,2] = -sint*cosp
	a[3,3] = cost

# transform into lens surface
	call gen_xfm (r, a, YES)

# refract ray (preserve scaling)
	q = r[1]*r[1] + r[2]*r[2] + r[3]*r[3]
	r[1] = n1 / n2 * r[1] / q
	r[2] = n1 / n2 * r[2] / q
	r[3] = sqrt (q - (r[1]*r[1] + r[2]*r[2]))

# transform out of lens surface	
	call gen_xfm (r, a, NO)

# translate coords to old system
	x[1] = xs + x0
	x[2] = ys + y0
	x[3] = zs + z0 + rcurv
end


procedure	flens (xp, yp, x, r)

double	xp, yp			# projected location in slitmask plane
double	x[3], r[3]

double	x0, y0, z0, rclens, rccoll

double	phi, theta, rp, hm		# slitmask params
double	cosp, sinp, cost, sint
double	n1, n2

begin
# calculate the angle vector
	rp = sqrt (xp*xp + yp*yp)
	hm = R_CURV - sqrt (R_CURV*R_CURV - rp*rp)
	theta = atan (rp / (D_1-hm))
	phi = atan2 (yp, xp)

	cosp = cos (phi)
	sinp = sin (phi)
	cost = cos (theta)
	sint = sin (theta)

	r[1] = cosp * sint
	r[2] = sinp * sint
	r[3] = cost

# calculate the reference point
	x[1] = xp
	x[2] = yp
	x[3] = -hm

	n1 = 1.00029d0		#  air
	n2 = 1.45637d0		# "silica_special" at 656nm
# First surface	
	x0 = 0.
	y0 = 278.511d0
	z0 = 76.2
	rclens = -10242.714d0
call eprintf ("x: %8.3f %8.3f %8.3f;     r: %8.5f %8.5f %8.5f\n")
call pargd (x[1]);   call pargd (x[2]);   call pargd (x[3])
call pargd (r[1]);   call pargd (r[2]);   call pargd (r[3])

	call sph_refr (x, r, x0, y0, z0, rclens, n1, n2)
call eprintf ("x: %8.3f %8.3f %8.3f;     r: %8.5f %8.5f %8.5f\n")
call pargd (x[1]);   call pargd (x[2]);   call pargd (x[3])
call pargd (r[1]);   call pargd (r[2]);   call pargd (r[3])

# Second surface
	z0 = 76.2 + 14.2239d0
	rclens = -4036.025d0
	call sph_refr (x, r, x0, y0, z0, rclens, n2, n1)
call eprintf ("x: %8.3f %8.3f %8.3f;     r: %8.5f %8.5f %8.5f\n")
call pargd (x[1]);   call pargd (x[2]);   call pargd (x[3])
call pargd (r[1]);   call pargd (r[2]);   call pargd (r[3])

# Collimator reflection:
	x0 = 0.
	y0 = 0.
	z0 = 76.2 + 14.2239d0 + 1921.577d0
	rccoll = -4009.9136d0
	call para_refl (x, r, x0, y0, z0, rccoll)
call eprintf ("x: %8.3f %8.3f %8.3f;     r: %8.5f %8.5f %8.5f\n")
call pargd (x[1]);   call pargd (x[2]);   call pargd (x[3])
call pargd (r[1]);   call pargd (r[2]);   call pargd (r[3])

end


# PARA_REFL: reflection from a parabolic surface (IT LIES! CURRENTLY SPH!!)

procedure	para_refl (x, r, x0, y0, z0, rcurv)

double	x[3]			# reference point
double	r[3]			# ray angle
double	x0, y0, z0		# x,y,z location of lens cent curv
double	rcurv			# radius of curvature

double	a[3,3]			# xform matrix
double	xc, yc			# x,y location of ray at cent of curv (zc=0)
double	xs, ys, zs		# x,y,z at lens surface
double	q, k
double	dz
double	phi, theta		# angles at lens surface

double	phin, thetan
double	rho
double	cosp, sinp, cost, sint

begin
# translate rays into lens system; don't need to change angles, only (x,y,z):
	dz = x[3] - z0 - rcurv
	xc = x[1] - x0 - dz*r[1]/r[3]
	yc = x[2] - y0 - dz*r[2]/r[3]

# Find the location of the intercept of ray with lens surface
# Spherical -- close enough at low angles
	q = r[1] * xc + r[2] * yc
	zs = r[3] * (-q + sqrt (q*q + rcurv*rcurv - (xc*xc + yc*yc)))	# Sph

# Parabolic -- not needed for x,y determination at low angles
	k = r[1]*r[1] + r[2]*r[2]	# problematic, as very small at times
	q = q - r[3]*rcurv
	zs = r[3]/k * (-q + sqrt (q*q + k * (2*rcurv*rcurv - (xc*xc + yc*yc))))

	xs = xc + r[1]/r[3] * zs
	ys = yc + r[2]/r[3] * zs

# get normal to lens surface
	rho = sqrt (xs*xs + ys*ys)
#	theta = atan (  rho / zs)		# spherical
	theta = atan (  rho / -rcurv)		# parabola
	phi = atan2 (ys , xs)

# calc the transform matrix
	thetan = theta
	phin = phi + HALFPI
	cosp = cos (phin)
	sinp = sin (phin)
	cost = cos (thetan)
	sint = sin (thetan)

	a[1,1] = cosp
	a[1,2] = sinp
	a[1,3] = 0.
	a[2,1] = -cost*sinp
	a[2,2] = cost*cosp
	a[2,3] = sint
	a[3,1] = sint*sinp
	a[3,2] = -sint*cosp
	a[3,3] = cost

# transform into lens surface
	call gen_xfm (r, a, YES)

# reflect ray (preserve scaling)
	r[3] = -r[3]

# transform out of lens surface	
	call gen_xfm (r, a, NO)

# translate coords to old system
	x[1] = xs + x0
	x[2] = ys + y0
	x[3] = zs + z0 + rcurv
end

