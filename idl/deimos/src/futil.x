# This file contains various utilities:
# All are "reasonable" precision, ie i/(o) in real, calc in double; fits return:
# The routines are:
#	GET_LSQF1: iterate LSq Fit for errors in x & y
#	GET_0LSQF1: get LSq Fit with errors in y only
#	GET_LSQF2: iterate LSq Fit for errors in x, y and z
#	GET_0LSQF2: get LSq Fit with errors in z only
#	GET_LSQF0: get LLSq Fit for constant offset
#	G_ELIM: Gaussian elimination
# real= VSUM1
# real= VSUM2
# real= VSUM3
# dbl =	DSUM1 -- vector sum (single vector)
# dbl =	DSUM2 -- vector sum (two vectors, ie dot product)
# dbl =	DSUM3 -- vector sum (three vectors, ie weighted dot product)

include	"futil.h"

#
# GET_LSQF1: iterate LSq Fit for errors in x & y
# NB: xerr, yerr are errors SQUARED.
#

procedure get_lsqf1 (x, y, xerr, yerr, weight, npts, niter, stats)

real	x[npts], y[npts]			# data vectors
real	xerr[npts], yerr[npts]			# error ** 2 vectors
real	weight[npts]				# additional weight factors
int	npts					# vector lengths
int	niter					# no. of iterations
real	stats[NFITPAR]				# returned fit params

int	i, j
real	slope, yintrcpt, me1
pointer	bufr, bufx, bufw
real	msq, wt, dm, db

begin
	call malloc (bufr, npts, TY_REAL)
	call malloc (bufx, npts, TY_REAL)
	call malloc (bufw, npts, TY_REAL)

# initial fit
	call get_0lsqf1 (x, y, weight, npts, stats)
	slope = SLOPE[stats]
	yintrcpt = YINCPT[stats]
	me1 = CHI[stats]
	call printf ("iteration: %2d  slope=%7.4f off=%6.2f (%7.3f) \n")
		call pargi (0)
		call pargr (slope)
		call pargr (yintrcpt)
		call pargr (me1)

# iterate
	do i = 1, niter {
		call altmr (x, Memr[bufr], npts, slope, yintrcpt)
		call asubr (y, Memr[bufr], Memr[bufr], npts)
		msq = slope * slope
		do j = 1, npts {
			wt = 1. / (yerr[j] + msq * xerr[j])
			Memr[bufw+j-1] = weight[j] * wt
			Memr[bufx+j-1] = x[j] + Memr[bufr+j-1] *
							slope * xerr[j] * wt
		}
		call get_0lsqf1 (Memr[bufx], Memr[bufr], Memr[bufw], npts, stats)
		dm = SLOPE[stats]
		db = YINCPT[stats]
		me1 = CHI[stats]
		slope = slope + dm
		yintrcpt = yintrcpt + db
		call printf ("iteration: %2d  slope=%7.4f off=%6.2f (%7.4f) \n")
			call pargi (i)
			call pargr (slope)
			call pargr (yintrcpt)
			call pargr (me1)
	}

	SLOPE[stats] = slope
	YINCPT[stats] = yintrcpt

	call mfree (bufr, TY_REAL)
	call mfree (bufx, TY_REAL)
	call mfree (bufw, TY_REAL)
end

#
# GET_0LSQF1: get LSq Fit with errors in y only; input, returned values are 
# real, fit calculation is double.

procedure get_0lsqf1 (x, y, w, npts, stats)

real	x[npts], y[npts]			# x, y vectors
real	w[npts]					# weight vector
int	npts					# vector length
real	stats[NFITPAR]				# returned

double	sumyy, sumxx, sumxy, sumx, sumy, sumw
double	a, b, det
real	ressq

double 	dsum1(), dsum2(), dsum3()

begin
	sumyy = dsum3 (y, y, w, npts)
	sumxx = dsum3 (x, x, w, npts)
	sumxy = dsum3 (x, y, w, npts)
	sumy = dsum2 (y, w, npts)
	sumx = dsum2 (x, w, npts)
	sumw = dsum1 (w, npts)

	det = sumw * sumxx - sumx * sumx
	if (det == 0.)
		call eprintf ("get_0lsqf1: zero determinant\n")
	a = (sumw * sumxy - sumx * sumy) / det
	b = (sumxx * sumy - sumx * sumxy) / det

# Work out stats:
	SLOPE[stats] = a
	YINCPT[stats] = b
	ressq = sumyy + a * (a*sumxx + 2. * (b*sumx - sumxy)) +
			b*(b*sumw - 2.*sumy)
	ressq = max (ressq, 0.)				# for roundoff
	CHI[stats] = sqrt (ressq / (npts - 2))		# is npts right?
	ESLOPE[stats] = CHI[stats] * sqrt (real (sumw / det))
	EYINCPT[stats] = CHI[stats] * sqrt (real (sumxx / det))
end

#
# VSUM1 -- vector sum (single vector)
#

real	procedure vsum1 (a, n)

real	a[n]			# input vector
int	n			# vector length

int	i
real	sum

begin
	sum = 0.
	do i = 1, n
		sum = sum + a[i]

	return (sum)
end

#
# VSUM2 -- vector sum (two vectors, ie dot product)
#

real	procedure vsum2 (a, b, n)

real	a[n], b[n]		# input vectors
int	n			# vector length

int	i
real	sum

begin
	sum = 0.
	do i = 1, n
		sum = sum + a[i] * b[i]

	return (sum)
end

#
# VSUM3 -- vector sum (three vectors, ie weighted dot product)
#

real	procedure vsum3 (a, b, c, n)

real	a[n], b[n], c[n]	# input vectors
int	n			# vector length

int	i
real	sum

begin
	sum = 0.
	do i = 1, n
		sum = sum + a[i] * b[i] * c[i]

	return (sum)
end

#
# DSUM1 -- vector sum (single vector)
#

double	procedure dsum1 (a, n)

real	a[n]			# input vector
int	n			# vector length

int	i
double	sum

begin
	sum = 0.
	do i = 1, n
		sum = sum + a[i]

	return (sum)
end

#
# DSUM2 -- vector sum (two vectors, ie dot product)
#

double	procedure dsum2 (a, b, n)

real	a[n], b[n]		# input vectors
int	n			# vector length

int	i
double	sum

begin
	sum = 0.
	do i = 1, n
		sum = sum + a[i] * b[i]

	return (sum)
end

#
# DSUM3 -- vector sum (three vectors, ie weighted dot product)
#

double	procedure dsum3 (a, b, c, n)

real	a[n], b[n], c[n]	# input vectors
int	n			# vector length

int	i
double	sum

begin
	sum = 0.
	do i = 1, n
		sum = sum + a[i] * b[i] * c[i]

	return (sum)
end

# Need to decide on naming convention

#
# GET_LSQF2: iterate LSq Fit to z=ax+by+c for errors in x, y and z.
# NB: xerr, yerr, zerr are errors SQUARED.
#

procedure get_lsqf2 (x, y, z, xerr, yerr, zerr, weight, npts, niter, stats)

real	x[npts], y[npts], z[npts]		# data vectors
real	xerr[npts], yerr[npts], zerr[npts]	# error ** 2 vectors
real	weight[npts]				# additional weight factors
int	npts					# vector lengths
int	niter					# no. of iterations
real	stats[NFITPAR]				# returned fit params

int	i, j
real	a, b, c, me1
pointer	bufr, bufx, bufy, bufw
real	asq, bsq, res, wt, da, db, dc

begin
	call malloc (bufr, npts, TY_REAL)
	call malloc (bufx, npts, TY_REAL)
	call malloc (bufy, npts, TY_REAL)
	call malloc (bufw, npts, TY_REAL)

# initial fit; NB needs expansion
	call get_0lsqf2 (x, y, z, weight, npts, stats)
	a = SLOPE1[stats]
	b = SLOPE2[stats]
	c = OFFSET[stats]
	me1 = CHI[stats]
#	call printf ("iteration: %2d  a=%7.4f b=%7.4f off=%6.2f (%7.3f) \n")
#		call pargi (0)
#		call pargr (a)
#		call pargr (b)
#		call pargr (c)
#		call pargr (me1)

# iterate
	do i = 1, niter {
		asq = a * a
		bsq = b * b
		do j = 1, npts {
			res = z[j] - (a * x[j] + b * y[j] + c)
			wt = 1. / (zerr[j] + asq * xerr[j] + bsq * yerr[j])
			Memr[bufr+j-1] = res
			Memr[bufw+j-1] = weight[j] * wt
			Memr[bufx+j-1] = x[j] + res * a * xerr[j] * wt
			Memr[bufy+j-1] = y[j] + res * b * yerr[j] * wt
		}
		call get_0lsqf2 (Memr[bufx], Memr[bufy], Memr[bufr], Memr[bufw], npts, stats)
		da = SLOPE1[stats]
		db = SLOPE2[stats]
		dc = OFFSET[stats]
		me1 = CHI[stats]
		a = a + da
		b = b + db
		c = c + dc
#		call printf ("iteration: %2d  a=%7.4f b=%7.4f off=%6.2f (%7.3f) \n")
#			call pargi (i)
#			call pargr (a)
#			call pargr (b)
#			call pargr (c)
#			call pargr (me1)
	}

	SLOPE1[stats] = a
	SLOPE2[stats] = b
	OFFSET[stats] = c

	call mfree (bufr, TY_REAL)
	call mfree (bufx, TY_REAL)
	call mfree (bufy, TY_REAL)
	call mfree (bufw, TY_REAL)
end

#
# GET_0LSQF2 -- calculate the zeroth order LLSq Fit for 2 independent variables,
# assumming errors in z only
#

	procedure get_0lsqf2 (x, y, z, w, npt, stats)

real	x[npt], y[npt]				# input coords
real	z[npt]					# ref. coord.
real	w[npt]					# weights
int	npt					# number of points
real	stats[NFITPAR]				# fit info struct

real	ga[4, 3]

double	dsum1(), dsum2(), dsum3()

begin
	ga[1,1] = dsum3 (x, x, w, npt)
	ga[2,1] = dsum3 (x, y, w, npt)
	ga[2,2] = dsum3 (y, y, w, npt)
	ga[3,1] = dsum2 (x, w, npt)
	ga[3,2] = dsum2 (y, w, npt)
	ga[4,1] = dsum3 (x, z, w, npt)
	ga[4,2] = dsum3 (y, z, w, npt)
	ga[4,3] = dsum2 (z, w, npt)
	ga[3,3] = dsum1 (w, npt)

	ga[1,2] = ga[2,1]
	ga[1,3] = ga[3,1]
	ga[2,3] = ga[3,2]

	call g_elim(ga, 3)

	SLOPE1[stats] = ga[4,1]
	SLOPE2[stats] = ga[4,2]
	OFFSET[stats]  = ga[4,3]
#need to define errors, me1
	EOFFSET[stats] = INDEF
	ESLOPE1[stats] = INDEF
	ESLOPE2[stats] = INDEF
	ME1[stats] = INDEF
end


#
# GET_LSQF0 -- calculate the offset (0th order coeff)
# fit equation y - x = b

	procedure get_lsqf0 (x, xp, xerr, xperr, w, npt, stats)

int	npt				# vector length
real	xp[npt]				# reference vector
real	x[npt]				# input vector
real	xerr[npt], xperr[npt]		# error vectors
real	w[npt]				# weight vector
real	stats[NFITPAR]			# fit info struct

double	sumxx, sumx, sumw
pointer	bufr, bufw

double	dsum1(), dsum2(), dsum3()

begin
	call malloc (bufr, npt, TY_REAL)
	call malloc (bufw, npt, TY_REAL)

	call asubr (xp, x, Memr[bufr], npt)
	call aaddr (xperr, xerr, Memr[bufw], npt)
	call adivr (w, Memr[bufw], Memr[bufw], npt)

	sumxx = dsum3 (Memr[bufr], Memr[bufr], Memr[bufw], npt)
	sumx = dsum2 (Memr[bufr], Memr[bufw], npt)
	sumw = dsum1 (Memr[bufw], npt)

	OFFSET[stats] = sumx / sumw
	ME1[stats] = sqrt (max (real ((sumxx - sumx*sumx/sumw) / (npt - 1)), 0.))
	EOFFSET[stats] = ME1[stats] / sqrt (real (sumw))

	call mfree (bufr, TY_REAL)
	call mfree (bufw, TY_REAL)
	return
end

#
# G_ELIM: procedure for gaussian elimination, n var's
#

procedure g_elim (a, n)

real	a[n+1,n]			# matrix to be solved
int	n				# number of variables

int	i, j, k
real	den, hold

begin
	do k = 1, n {
		den = a[k,k]
		if (den == 0.) {		# look for non-zero: switch
			do j = (k+1), n {
				if (a[k,k] != 0.) {
					do i = k, (n+1) {
						hold = a[i,j]
						a[i,j] = a[i,k]
						a[i,k] = hold
					}
				den = a[k,k]
				}
			}
			if (den == 0.)			# if still zero, skip
				next
		}
		do i = k, (n+1)
			a[i,k] = a[i,k] / den
		do j = 1, n
			if (j != k) {
				den = a[k,j]
				do i = k, (n+1)
					a[i,j] = a[i,j] - a[i,k] * den
			}
	}
end

