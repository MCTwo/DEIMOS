
define	TVMINMAG	20.8		# min mag to print for TV pos

# TDB immed: insert the NULL values appropriately
# OR, use the cfitsio routines
#
# TBD:  get various params into proper routines
# The malloc's in selector should be replaced by salloc's
# Note that slit-sep should be _added_ to all required lengths,
# The whole auto-selector issue needs to be revisited, with real lengths added.
# Resolution of slit lengths (extend, resolve conficts)
# Binding, slit editing
# Mapping slits to celestial coordinates  -- DONE!
# DEFINE the LEN issue

# Mapping slits to metal [mostly done] ****

# NB: it is probably best to get the central slit coords, etc, and then
# map onto the metal as we did before.  There are a few reasons for doing it
# this way, part of which is the curved-slit problem (assuming the object is
# more likely to be near the center than the edge of the slits); also, this
# is the direction things will likely go in the end; and 3rd it makes dealing
# with the PA=INDEF objects a bit easier, since they can be cut as vertical
# on the metal.

include	<math.h>
include	"deimos.h"
include	"dsimulator.h"

procedure	t_dsimulator()

char	objfile[SZ_FNAME]			# ID, RA, Dec, equ, prior,
char	output[SZ_FNAME]			# output file name
char	mdf[SZ_FNAME]				# mask design file name
char	plotfile[SZ_FNAME]			# optional plotfile
bool	std_fmt				# Standard text input format?

pointer	fda, fdb, fdc				# in/out file descriptors

pointer	indat			# pointer to instr/tel data structure

int	ntarg
pointer	tdat			# pointer to data structure for targets

int	nslit
pointer	sdat			# pointer to data structure for slits

# TMP? XXX find better location
real	delpa

bool	clgetb()
pointer	open()

begin
# Initialize -- no targets or slits
	nslit = 0
	ntarg = 0
	sdat = 0		# Used as flag for whether slits defined

# Read in file names
	call clgstr ("objfile", objfile, SZ_FNAME)
	call clgstr ("output", output, SZ_FNAME)
	call clgstr ("mdf", mdf, SZ_FNAME)
	call clgstr ("plotfile", plotfile, SZ_FNAME)
	std_fmt = clgetb ("std_format")


# Read in telescope, default data (some may get updated through ingested files)
	call data_init (indat)

	if (std_fmt) {
# Open the input list of targets:
		fda = open (objfile, READ_ONLY, TEXT_FILE)
		fdb = open (output, NEW_FILE, TEXT_FILE)
		fdc = open (mdf, NEW_FILE, TEXT_FILE)

# Read in the target data:
		call targ_init (fda, tdat, ntarg, indat)

# refract coordinates (ASSUME PA does not need refraction); to first order
# this is true since PA's are relative to rotator PA, similarly affected.

		call refr_coords (tdat, ntarg, indat)


# Calc the location of the telescope axis
		call fld2telax (indat)

# Calculate position wrt telescope axis:
		call tel_coords (tdat, ntarg, indat)
	} else {
# Open the input list of targets:
		fdc = open (mdf, NEW_FILE, TEXT_FILE)
		fdb = open (output, APPEND, TEXT_FILE)

#		PROJ_LEN(indat) = NO	# should not be needed -- a check

# Read in the target data:
		call marc (objfile, tdat, ntarg, sdat, nslit, indat)

# refract coordinates (ASSUME PA does not need refraction); to first order
# this is true since PA's are relative to rotator PA, similarly affected.

		call refr_coords (tdat, ntarg, indat)
		call refr_coords (sdat, nslit, indat)

# Calc the location of the telescope axis
		call fld2telax (indat)

# Calculate position wrt telescope axis:
		call tel_coords (tdat, ntarg, indat)
		call tel_coords (sdat, nslit, indat)
		call marc2 (sdat, nslit, tdat, ntarg)
call eprintf ("ntarg, nslit = %d %d\n")
call pargi (ntarg)
call pargi (nslit)
	}


# Ready to display and be interactive
	call gselect_targ (tdat, ntarg, sdat, nslit, indat)

# Save output?

# Calculate the SLITS in celestial coordinates.
# It is here that slit lengths, edges must be resolved.  Does binding go here?
# Or is this a step that should actually be included above?
#

	if (std_fmt && nslit == 0) {
		call eprintf ("Running gen_slit -- assumed correct action\n")
		call gen_slits (tdat, ntarg, sdat, nslit, indat, YES)
#####
call eprintf ("Slit lengthening INCLUDED\007!\n")
#####
	}


call eprintf (" calling sky_coord \n")
# Convert adopted slits into coordinates on the sky:
	call sky_coords (sdat, nslit, indat)
call eprintf (" calling unrefr_coord \n")
	call unrefr_coords (sdat, nslit, indat)		# XXX; should go elsewh


## For symmetry, we now recalculate the telescope coords and then mask coords.
## In the end, this operation may be performed at Keck immediately prior to
## mask fabrication.
#  Note that _refraction_ correction should probably go here, as a differential
#  correction, as it is tiny and specific to actual time and date of observation

call eprintf (" calling tel_coord \n")

# Calculate position wrt telescope axis:
	call tel_coords (sdat, nslit, indat)

 
 call eprintf (" calling mask_coords \n")
# Generate the mask coordinates for the slits.

 	call mask_coords (sdat, nslit)


#
call eprintf (" calling write_design ...\n")
#	call write_fitsio (mdf, indat, tdat, ntarg, sdat, nslit, objfile)
  	call write_design (indat, tdat, ntarg, sdat, nslit, fdc, objfile)
# call eprintf (" ASCII Surfcam Format ...\n")
#	call write_surfcam (sdat, nslit, fdc, "DSIM1145")

call eprintf ("Close output\n")
	call close (fdc)


call eprintf ("list_dat\n")
	call list_dat (fdb, indat, tdat, ntarg, sdat, nslit, delpa, std_fmt, plotfile)

	call close (fdb)


####### COPIED FROM MAPMASK -- correction to PA
# One final kludge: precession will cause a change in PA -- until precession
# is worked into the whole package, calculate the updated PA for now.
# Precession approximation from Lang, using m,n for 1975.0

	delpa = atan ((EPOCH(indat) - STD_EQX(indat)) * 9.7157e-5 *
				sin (RA0_FLD(indat)) / cos (DEC0_FLD(indat)))
call eprintf ("DELPA: %6.3f\n")
call pargr (RADTODEG(delpa))





# XXX
call eprintf (" [recall: still must check proper x,y and tan mapping]\n")
call eprintf ("MDF file written: %-s\n")
	call pargstr (mdf)
end

	

#
# FLD2TELAX:  from field center and rotator PA, calc coords of telescope axis
#

procedure	fld2telax (indat)

pointer	indat

double	r, pa_fld

double	cosd, sind, cosa, sina, cost, sint, cosr, sinr

begin
# convert field center offset (arcsec) to radians
	r = DEGTORAD(sqrt (FLDCEN_X*FLDCEN_X + FLDCEN_Y*FLDCEN_Y) / 3600.)

# get PA of field center
	pa_fld = atan2 (FLDCEN_Y, FLDCEN_X)

	cosr = cos (r)
	sinr = sin (r)
	cosd = cos (DEC_FLD(indat))
	sind = sin (DEC_FLD(indat))

	cost = cos (PA_ROT(indat) - pa_fld)
	sint = sin (PA_ROT(indat) - pa_fld)

	sina = sinr * sint / cosd		# ASSUME not at dec=90
	cosa = sqrt (1. - sina*sina)

	RA_TEL(indat) = RA_FLD(indat) - asin (sina)
	DEC_TEL(indat) = asin ((sind*cosd*cosa - cosr*sinr*cost) /
				(cosr*cosd*cosa - sinr*sind*cost))
# call eprintf ("DEBUG: RA,DEC: %8f %8f\n")
# call pargd (RADTODEG(RA_TEL(indat)))
# call pargd (RADTODEG(DEC_TEL(indat)))

end

############################################################################
#

#
# SELECTOR: Does an auto selection of slits
# Should include an option for weighting to keep things toward the center.
# Note that y's sent to sel_rank are relative to starting y
# have it run from bottom to top
# -- need to work in segments to accomodate currently selected objects
# -- there was something else ...

procedure	selector (indat, tdat, ntarg, nlist, minsep, psum)

pointer	indat
pointer	tdat
int	ntarg
int	nlist		# list to work on
real	minsep		# XXX min. separation -- probably should be in DEFDAT
int	psum		# sum of selected priorities (returned)

int	nopt, npre		# number of options, prev. selected objects
int	i, ndx
int	ix			# starting index for search (saves time)
int	nselect				# Number of selected slits
real	xlow, xupp, xskip
pointer	bufx1, bufx2			# TMP? buffers for pre-sel. objs
pointer	bufn, bufx, bufp, bufsel	# TMP, for now

begin
	nopt = 0
	npre = 0
	do i = 0, ntarg-1 {
		if (SEL(tdat,i) == YES) {
			npre = npre + 1
		} else if (SAMPL(tdat,i) == nlist && STAT(tdat,i) == YES) {
			nopt = nopt + 1
		}
	}
	call malloc (bufx1, npre, TY_INT)
	call malloc (bufx2, npre, TY_REAL)
	call malloc (bufn, nopt, TY_INT)
	call malloc (bufx, nopt, TY_REAL)
	call malloc (bufp, nopt, TY_INT)
	call malloc (bufsel, nopt, TY_INT)

# Grep on previously selected objects and suitable options; fill vectors
	nopt = 0
	npre = 0
	do i = 0, ntarg-1 {
		if (PCODE(tdat,i) == CODE_GS)		# GS's don't take space
			next
		if (SEL(tdat,i) == YES) {
#			Memr[bufx1+npre] = XARCS(tdat,i) - LEN1(tdat,i)
#			Memr[bufx2+npre] = XARCS(tdat,i) + LEN2(tdat,i)
			Memr[bufx1+npre] = X1(tdat,i)
			Memr[bufx2+npre] = X2(tdat,i)
			npre = npre + 1
		} else if (SAMPL(tdat,i) == nlist && STAT(tdat,i) == YES && PCODE(tdat,i) > 0) {
			Memi[bufn+nopt] = i	# INDEX(tdat,i) XXX
			Memr[bufx+nopt] = XARCS(tdat,i)
			Memi[bufp+nopt] = PCODE(tdat,i)
			nopt = nopt + 1
		}
	}

# Sort the two lists
	call sel_sort (Memr[bufx1], Memr[bufx2], npre,
				Memi[bufn], Memr[bufx], Memi[bufp], nopt)

# The number of "gaps" to search is npre+1
	ndx = 0
	xlow = XLOW_LIM
	xskip = 0.
	nselect = 0			# triggers init in sel_rank
	if (nopt > 0) {
	    do i = 0, npre {
		if (i < npre) {
			xupp = Memr[bufx1+i]
			xskip = Memr[bufx2+i] - Memr[bufx1+i]
		} else {
			xupp = XUPP_LIM
		}

## old ...
#		if (xupp <= xlow)
#			next
#		call sel_rank (Memr[bufx], Memi[bufp], Memi[bufn],
#		Memi[bufsel], nopt, ix, xlow, xupp, minsep, nselect)
#		xlow = xupp + xskip

		if (xupp > xlow) {
			call sel_rank (tdat, indat, Memi[bufn],
			Memi[bufsel], nopt, ix, xlow, xupp, minsep, nselect)
		}

		xlow = xupp + xskip
	    }
	}


#...select the mask slits
	if (nselect > 0) {
		do i = 0, nselect-1 {
			SEL(tdat,Memi[bufsel+i]) = YES
		}
	}

	psum = 0
	do i = 0, ntarg-1 {
		if (SEL(tdat,i) == YES)
			psum = psum + max (PCODE(tdat,i), 0)	# NO GS, AS
	}

	call mfree (bufsel, TY_INT)
	call mfree (bufp, TY_INT)
	call mfree (bufx, TY_REAL)
	call mfree (bufn, TY_INT)
	call mfree (bufx2, TY_REAL)
	call mfree (bufx1, TY_INT)
end


# SEL_RANK: Select slits with priority ranking.
# The scheme is to find the next possible slit, and then to look up to one
# min-slit width away for higher-priority objects. The higher priorities are
# down-weighted depending on their distance.
#

procedure	sel_rank (tdat, indat, tndex, sel, npt, isel, xlow, xupp, minsep, nsel)

pointer	tdat
pointer	indat
int	tndex[npt]			# Index of selectable slits
int	sel[npt]			# selected objects
int	npt				# Number of objects
int	isel 				# starting index
real	xlow, xupp			# xrange to fill
real	minsep				# minimum separation (arcsec)
int	nsel				# Number of selected objects

int	i, j, ndx

real	x, xj, xnext, xlook, xstop, xlast
real	len
real	prisel, prinorm

begin
# make sure we initialize to the right index ...
	if (nsel <= 0) {
		isel = 0
		nsel = 0
	}

# Can we fit a minimum slit in here?
	if (xupp - xlow < minsep)		# probably too restrictive
		return


# Start at half a slit length; stop inside half slit length
	ndx = tndex[npt]
	x = XARCS(tdat,ndx)
	xstop = min (x, xupp-0.5*minsep)	# last target or upper limit
	xnext = xlow + 0.5 * minsep
	xlast = xlow

# Loop through to end
	for (i = isel + 1; i <= npt; i = i + 1) {
		ndx = tndex[i]
		x = XARCS(tdat,ndx)
		if (x < xnext)
			next
		if (X1(tdat,ndx) < xlast)
			next

		if (x > xstop) {
			isel = i - 1
			break
		}

		isel = i
		len = X2(tdat,ndx) - X1(tdat,ndx)
		prisel = PCODE(tdat,ndx) / (x - xlast) / len
# Now look for higher priority to win out, over range (xlast,xlook)
		xlook = min (x+minsep, xstop)
		if (isel < npt) {
			do j = isel+1, npt {
				ndx = tndex[j]
				if (X1(tdat,ndx) >
					X2(tdat,isel)+SLIT_GAP(indat)) {
					next		# There is no conflict
				}		# XXX but prisel gets higher?
				if (X2(tdat,ndx) > x_upp)
					next		# XXX Can't use
				xj = XARCS(tdat,ndx)
				if (xj >= xlook)
					break
				len = X2(tdat,ndx) - X1(tdat,ndx)
				prinorm = PCODE(tdat,ndx) / (xj - xlast) / len
				if (prinorm > prisel) {
					x = xj
					isel = j
					prisel = prinorm
				}
			}
		}

		nsel = nsel + 1
		ndx = tndex[isel]
		sel[nsel] = ndx
		xlast = X2(tdat,ndx)
		xnext = xlast + 0.5 * minsep
		i = isel			# Reset search start point

	}
end

#
# SEL_SORT: sort the selected/selection lists
#

procedure	sel_sort (px1, px2, npre, index, x, pri, nopt)

real	px1[npre], px2[npre]		# x-pre limits to sort
int	npre				# Number of prev. selections

int	index[nopt]			# Index of optional objects
real	x[nopt]				# x-opt list to sort
int	pri[nopt]			# Priority
int	nopt				# Number of optional objects

int	i, j
int	ihold, phold
real	xhold

begin

# Sort the preselected list in x (low-to-high)

	do i = 1, npre-1 {
	    do j = 1, npre-i {
		if (px1[j] > px1[j+1]) {
			xhold = px1[j+1]
			px1[j+1] = px1[j]
			px1[j] = xhold
			xhold = px2[j+1]
			px2[j+1] = px2[j]
			px2[j] = xhold
		}
	    }
	}

# Sort the list of optional objects in x (low-to-high)
	do i = 1, nopt-1 {
	    do j = 1, nopt-i {
		if (x[j] > x[j+1]) {
			xhold = x[j+1]
			x[j+1] = x[j]
			x[j] = xhold
			ihold = index[j+1]
			index[j+1] = index[j]
			index[j] = ihold
			phold = pri[j+1]
			pri[j+1] = pri[j]
			pri[j] = phold
		}
	    }
	}
end


#
# MASK_COORDS:  Convert (x,y) on sky to xmm,ymm on slitmask
# NB: Requires instrument params, INCLUDING THE APPROP. TEL FOC LENGTH
# Assumes that the XARCS,YARCS of the slit ends are the tan projections
#

procedure	mask_coords (sdat, nslit)

pointer	sdat
int	nslit

double	xfp, yfp			# x,y points in FP (tan projection)
double	xsm, ysm			# x,y points on the mask
double	pa

int	i
double	xoff, yoff			# offset, telaxis to origin of slitmask
double	asec_rad
real	sina, cosa

begin
	asec_rad = 206264.8D0

# offset from telescope axis to slitmask origin, IN SLITMASK COORDS
	yoff = ZPT_YM * (1. - cos (DEGTORAD(M_ANGLE)))
	yoff = 0.	# XXX check!  Am not sure where the above comes from
	xoff = 0.

	do i = 0, nslit-1 {
# XXX For now, carry through the RELPA thing; in end, must be specified!
			if (RELPA(sdat,i) != INDEF) {
				cosa = cos (RELPA(sdat,i))
				sina = sin (RELPA(sdat,i))
			} else {
				cosa = 1.
				sina = 0.
			}
#			cosa = cos (RELPA(sdat,i))	# XXX 
#			sina = sin (RELPA(sdat,i))	# XXX 


# This is a recalculation ... prob not needed
			X1(sdat,i) = XARCS(sdat,i) - LEN1(sdat,i) * cosa * FLIP
			Y1(sdat,i) = YARCS(sdat,i) - LEN1(sdat,i) * sina
			X2(sdat,i) = XARCS(sdat,i) + LEN2(sdat,i) * cosa * FLIP
			Y2(sdat,i) = YARCS(sdat,i) + LEN2(sdat,i) * sina

# XXX cuidado!  I am not sure that the tan-projection of the rel PA is the
# same as the rel PA -- MUST CHECK! (This code comes from gen_slits)


# The focal plane coordinates are now simply a tan projection of (x,y) arcsec
# Need to verify that these are truly symmetric:
#		xfp = FL_TEL * tan (DEGTORAD(X1(sdat,i)/3600.))
#		yfp = FL_TEL * tan (DEGTORAD(Y1(sdat,i)/3600.)) / cos (DEGTORAD(X1(sdat,i)/3600.))


# X1,Y1 are now tan projections already!

		xfp = FL_TEL *  X1(sdat,i) / asec_rad
		yfp = FL_TEL * (Y1(sdat,i) - 0.5*SLWID(sdat,i)) / asec_rad
		pa = 0.
		call gnom_to_dproj (xfp, yfp, xfp, yfp)		# (allowed)
		call proj_to_mask (xfp, yfp, pa, xsm, ysm, pa)

		XMM1(sdat,i) = xsm + xoff
		YMM1(sdat,i) = ysm + yoff

		xfp = FL_TEL *  X2(sdat,i) / asec_rad
		yfp = FL_TEL * (Y2(sdat,i) - 0.5*SLWID(sdat,i)) / asec_rad
		pa = 0.
		call gnom_to_dproj (xfp, yfp, xfp, yfp)		# (allowed)
		call proj_to_mask (xfp, yfp, pa, xsm, ysm, pa)

		XMM2(sdat,i) = xsm + xoff
		YMM2(sdat,i) = ysm + yoff

		xfp = FL_TEL *  X2(sdat,i) / asec_rad
		yfp = FL_TEL * (Y2(sdat,i) + 0.5*SLWID(sdat,i)) / asec_rad
		pa = 0.
		call gnom_to_dproj (xfp, yfp, xfp, yfp)		# (allowed)
		call proj_to_mask (xfp, yfp, pa, xsm, ysm, pa)

		XMM3(sdat,i) = xsm + xoff
		YMM3(sdat,i) = ysm + yoff

		xfp = FL_TEL *  X1(sdat,i) / asec_rad
		yfp = FL_TEL * (Y1(sdat,i) + 0.5*SLWID(sdat,i)) / asec_rad
		pa = 0.
		call gnom_to_dproj (xfp, yfp, xfp, yfp)		# (allowed)
		call proj_to_mask (xfp, yfp, pa, xsm, ysm, pa)

		XMM4(sdat,i) = xsm + xoff
		YMM4(sdat,i) = ysm + yoff

	}

	call metal_check (sdat, nslit)

## Perhaps we want to force YMM4-YMM1 == YMM3-YMM2; for non-tilted slits,
## this should produce a cleaner edge; otherwise, jumps can occur.

end

#
# METAL_CHECK:  Checks to make sure metal limits are not violated.  Fairly
# temporary and filled with hardcodes

procedure	metal_check (sdat, nslit)

pointer	sdat
int	nslit

int	i
real	xmin, xmax, ymin, ymax

begin
	do i = 0, nslit-1 {

		xmin = min (XMM1(sdat,i), XMM4(sdat,i))
		xmax = max (XMM2(sdat,i), XMM3(sdat,i))

		ymin = min (YMM1(sdat,i), YMM2(sdat,i))
		ymax = max (YMM4(sdat,i), YMM4(sdat,i))

		if (xmin < -373.) {
			call eprintf ("slit=%d; xmin=%6f \n")
				call pargi (i)
				call pargr (xmin)
			call fatal (0, "xmin < -373.")
		}

		if (xmax > 373.) {
			call eprintf ("slit=%d; xmax=%6f \n")
				call pargi (i)
				call pargr (xmax)
			call fatal (0, "xmax > 373.")
		}

# XXX should check for guide stars here
#		if (ymin < 2.) {
#			call eprintf ("slit=%d; ymin=%6f \n")
#				call pargi (i)
#				call pargr (ymin)
#			call fatal (0, "ymin < 2.")
#		}

		if (ymax > 225.17) {
			call eprintf ("slit=%d; ymax=%6f \n")
				call pargi (i)
				call pargr (ymax)
			call fatal (0, "ymax > 225.17")
		}

# y+x < 490.
		if (XMM3(sdat,i)+YMM3(sdat,i) > 496.) {
			call eprintf ("offending slit = %d \n")
				call pargi (i)
			call fatal (0, "beyond cut edge!")
		}

	}
end



#
# SKY_COORDS: Convert xarcs,yarcs in tel coords onto sky
#   Note that this routine, called infrequently, does not need to be efficient.
#

procedure	sky_coords (sdat, nslit, indat)

pointer	sdat
int	nslit
pointer	indat

int	i
double	r			# radius of object from tel-axis
double	phi			# PA on sky from tel-axis
double	sind, sina		# sin of dec, delta-RA
double	ra0, dec0, pa0		# RA, Dec and PA on axis

double	x, y

begin
	ra0  = RA_TEL(indat)
	dec0 = DEC_TEL(indat)
	pa0  = PA_ROT(indat)

	do i = 0, nslit-1 {
		x = 0.5 * (X1(sdat,i) + X2(sdat,i))
		y = 0.5 * (Y1(sdat,i) + Y2(sdat,i))

		r = sqrt (x*x + y*y)
		r = atan (r/206264.8D0)

		phi = pa0 - atan2 (y, x)	# WORK

		sind = sin (dec0) * cos (r) + cos (dec0) * sin (r) * cos (phi)

		sina = sin (r) * sin (phi) / sqrt (1. - sind*sind)

		DEC(sdat,i) = asin (sind)
		RA(sdat,i) = ra0 + asin (sina)
# PA(sdat,i) = already assigned 

# calc the centers and lengths of the slits

		XARCS(sdat,i) = 0.5 * (X1(sdat,i) + X2(sdat,i))
		YARCS(sdat,i) = 0.5 * (Y1(sdat,i) + Y2(sdat,i))

# XXX NB: by convention, slit length will be defined as TOTAL length
		x = X2(sdat,i) - X1(sdat,i)
		y = Y2(sdat,i) - Y1(sdat,i)
		LEN1(sdat,i) = 0.5 * sqrt (x*x + y*y)
		LEN2(sdat,i) = LEN1(sdat,i)

	}
end


#
# SLIT_SORT: sort the slits in x; keep track of the indices only
#

procedure	slit_sort (x, index, n)

real	x[n]				# x-pre limits to sort
int	index[n]			# Index of slit
int	n				# Number of prev. selections


int	i, j
int	ihold
real	xhold

begin

# Sort the slit list in x (low-to-high)

	do i = 1, n-1 {
	    do j = 1, n-i {
		if (x[j] > x[j+1]) {
			xhold = x[j+1]
			x[j+1] = x[j]
			x[j] = xhold
			ihold = index[j+1]
			index[j+1] = index[j]
			index[j] = ihold
		}
	    }
	}
end

#
# LEN_SLITS: adjust slit lengths to fit -- perhaps should be integral part of
# gen_slits
# XXX TBD: report conflicts?
#

procedure	len_slits (tdat, ntarg, sdat, nslit, indat)

pointer	tdat			# pointer to target struct
int	ntarg

pointer	sdat			# pointer to slits struct
int	nslit

pointer	indat			# pointer to instrument struct

int	i
pointer	sp
pointer	bufx			# pointer to x-pos
pointer	bufi			# pointer to index vector

int	pc1, pc2
int	ndx1, ndx2
real	xlow, xupp, xcen
real	del1, del2, tana

real	dxlow, dxupp		# corrections to avoid overlap
real	yas, dxavg
pointer	fdx			# pointer to XPROJ mapping
pointer	asfx, asfy		# pointers to surf fits (only asfx used)

real	gseval()
pointer	open()
begin
	call smark (sp)
	call salloc (bufx, nslit, TY_REAL)
	call salloc (bufi, nslit, TY_INT)

	call amovr (XARCS(sdat,0), Memr[bufx], nslit)
	call amovi (INDEX(sdat,0), Memi[bufi], nslit)

	call slit_sort (Memr[bufx], Memi[bufi], nslit)

## XXX must add extension to mask edge in here.

# Open X-Projection mapping
	fdx = open (XPROJ_MAP, READ_ONLY, TEXT_FILE)
	call gs_ingest (fdx, asfx, asfy)

	do i = 0, nslit-2 {
		ndx1 = Memi[bufi+i]
		ndx2 = Memi[bufi+i+1]
		pc1 = PCODE(sdat,ndx1)
		pc2 = PCODE(sdat,ndx2)

# If both are alignment boxes, just go on ...
		if (pc1 == CODE_AS && pc2 == CODE_AS)		# no problem
			next

# We will need to recalculate something ...

		xlow = X2(sdat,ndx1) + SLIT_GAP(indat)
		xupp = X1(sdat,ndx2) - SLIT_GAP(indat)
		xcen = 0.5 * (xlow + xupp)

		yas = Y2(sdat,ndx1)
		dxlow = gseval (asfx, xcen, yas)
		yas = Y1(sdat,ndx2)
		dxupp = gseval (asfx, xcen, yas)
		dxavg = 0.5 * (dxupp + dxlow)
		dxlow = dxlow - dxavg
		dxupp = dxupp - dxavg
call eprintf ("%6.3f %6.3f\n")
call pargr (dxlow)
call pargr (dxupp)

		if (pc1 == CODE_AS) {
			del1 = 0.
			del2 = X1(sdat,ndx2) - xlow - (dxupp - dxlow)
		} else if (pc2 == CODE_AS) {
			del1 =  xupp - X2(sdat,ndx1) + (dxlow - dxupp)
			del2 = 0.
		} else {
			del1 = xcen - 0.5*SLIT_GAP(indat) - X2(sdat,ndx1) + dxlow
			del2 = X1(sdat,ndx2) - (xcen + 0.5*SLIT_GAP(indat)) - dxupp
		}



		X2(sdat,ndx1) = X2(sdat,ndx1) + del1
		if (del1 != 0. && RELPA(sdat,ndx1) != INDEF) {
			tana = tan (RELPA(sdat,ndx1))
			Y2(sdat,ndx1) = Y2(sdat,ndx1) + del1 * FLIP * tana
		}

		X1(sdat,ndx2) = X1(sdat,ndx2) - del2
		if (del2 != 0. && RELPA(sdat,ndx2) != INDEF) {
			tana = tan (RELPA(sdat,ndx2))
			Y1(sdat,ndx2) = Y1(sdat,ndx2) - del2 * FLIP * tana
		}

# The centers and lengths have now changed -- but defer to sky_coords()

# XXX report status here
#		if (del1 < 0) {			# conflict
#		} else if (xlow < xupp) {		# lengthen
#		}
		
	}

# XXX Check against objects, etc, and reflect status
		
	call close (fdx)
	call sfree (sp)
end



#
# LIST_DAT: Text listing of relevant params
#
define	SZ_WRD	60	# length of word for parsing

procedure	list_dat (fd, indat, tdat, ntarg, sdat, nslit, delpa, std_fmt, plotfile)

pointer	fd			# file descriptor of output file
pointer	indat			# pointer to instr/telesc data struct

pointer	tdat			# pointer to targets data struct
int	ntarg

pointer	sdat			# pointer to slits data struct
int	nslit

real	delpa			# updated PA for epoch

int	std_fmt			# Standard format to input list?
char	plotfile[ARB]		# name of (opt) plot file
pointer	fdp			# optional plotfile

char	outline[SZ_LINE]	# remainder of reconstituted output line
char	pars1[SZ_WRD]		# parsed word for output
char	pars2[SZ_WRD]		# parsed word for output
char	pars3[SZ_WRD]		# parsed word for output
char	pars4[SZ_WRD]		# parsed word for output
char	pars5[SZ_WRD]		# parsed word for output

char	ident[SZ_ID]
int	stat
real	xtv, ytv		# Coordinates in the TV system
real	xas, yas		# x,y in arcsec
double	ra_tv, dec_tv		# RA, Dec of TV field
pointer	fdtv1, fdtv2		# pointers to input files for TV surf map
pointer	tvsf1x, tvsf1y, tvsf2x, tvsf2y	# pointers to surface fits for TV

int	ngstar			# number of guide stars
pointer	buftvn, buftvm		# pointers to TV star vectors

int	i, ndx

bool	strne()
int	sscan(), nscan()
real	gseval()
pointer	open()

begin
# Print out basic info:

	call fprintf (fd, "# Mask name, center:\n%-16s %12.2h %12.1h %7.1f PA=%6.3f ##\n")
		call pargstr (DESNAME(indat))
		call pargd (RADTODEG(RA0_FLD(indat))/15.d0)
		call pargd (RADTODEG(DEC0_FLD(indat)))
		call pargd (STD_EQX(indat))
		call pargd (RADTODEG(PA_ROT(indat)))

# Print out TV field (XXX -- testing it here):
	xas = 102.2
	yas = 196.4
	call tel2radec (indat, xas, yas, ra_tv, dec_tv)
	call fprintf (fd, "\n#  Guider center:  %12.2h %12.1h\n")
		call pargd (RADTODEG(ra_tv)/15.)
		call pargd (RADTODEG(dec_tv))

	if (std_fmt == YES) {
		call fprintf (fd, "\n# Selected Objects:\n\n")
		do i = 0, ntarg-1 {
			if (SEL(tdat,i) == NO)
				next

# we now want to reconstitute the output line, changing only the SEL CODE
			stat = sscan (DATLINE(tdat,i))
			call gargwrd (ident, SZ_ID)
			call gargwrd (pars1, SZ_WRD)	# RA
			call gargwrd (pars2, SZ_WRD)	# Dec
			call gargwrd (pars3, SZ_WRD)	# Equinox
			call gargwrd (pars4, SZ_WRD)	# mag
			call gargwrd (pars5, SZ_WRD)	# pband
			call gargi (stat)	# skip Pcode
			call gargi (stat)	# skip Sample (present if more)
			call gargi (stat)	# skip SelCode (present if more)
			call gargstr (outline, SZ_LINE)		# copy rest
			if (nscan() < 10)
				call strcpy ("", outline, SZ_LINE)

			call fprintf (fd, "%-16s %12s %12s %s %5s %s %4d %1d %1d%s\n")
				call pargstr (ident)
				call pargstr (pars1)
				call pargstr (pars2)
				call pargstr (pars3)
				call pargstr (pars4)
				call pargstr (pars5)
				call pargi (PCODE(tdat,i))
				call pargi (SAMPL(tdat,i))
				call pargi (SEL(tdat,i))
				call pargstr (outline)
		}
	}

# OK, now look for guide stars:

	ngstar = 0
	do i = 0, ntarg-1 {
		if (SEL(tdat,i) == NO)
			next
		if (PCODE(tdat,i) != CODE_GS)
			next
		if (MAG(tdat,i) > TVMINMAG)		# XXX kludge
			next
		ngstar = ngstar + 1
	}

	call malloc (buftvn, ngstar, TY_INT)		# index
	call malloc (buftvm, ngstar, TY_REAL)		# TV mag

# first open the TV astrometry mappings
	fdtv1 = open (TV_MIRR_MAP, READ_ONLY, TEXT_FILE)
	fdtv2 = open (TV_MASK_MAP, READ_ONLY, TEXT_FILE)

	call gs_ingest (fdtv1, tvsf1x, tvsf1y)
	call gs_ingest (fdtv2, tvsf2x, tvsf2y)

	call fprintf (fd, "\n# Selected Guide Stars:\n\n")

	ndx = 0
	do i = 0, ntarg-1 {
		if (SEL(tdat,i) == NO)
			next
		if (PCODE(tdat,i) != CODE_GS)
			next
		if (MAG(tdat,i) > TVMINMAG)		# XXX kludge
			next

		stat = sscan (DATLINE(tdat,i))
		call gargwrd (ident, SZ_ID)

		xas = XARCS(tdat,i)
		yas = YARCS(tdat,i)
		if (yas < 176.) {
			xtv = gseval (tvsf1x, xas, yas)
			ytv = gseval (tvsf1y, xas, yas)
		} else {
			xtv = gseval (tvsf2x, xas, yas)
			ytv = gseval (tvsf2y, xas, yas)
		}

#		call fprintf (fd, "# %-11s %7.3f %7.3f %5.2f_%1s  TV(x,y)= %6.1f %6.1f ##\n")
		call fprintf (fd, "# %-11s %11.2h %11.1h %5.2f_%1s  TV(x,y)= %6.1f %6.1f ##\n")
			call pargstr (ident)
			call pargd (RADTODEG(RA0(tdat,i))/15.d0)
			call pargd (RADTODEG(DEC0(tdat,i)))
#			call pargr (XARCS(tdat,i))
#			call pargr (YARCS(tdat,i))
			call pargr (MAG(tdat,i))
			call pargc (PBAND(tdat,i))
			call pargr (xtv)
			call pargr (ytv)
# Does TV star really fall on field?
		if (xtv > 3. && xtv < 1021. && ytv > 3. && ytv < 1019.) {
			Memi[buftvn+ndx] = i
			Memr[buftvm+ndx] = MAG(tdat,i)
			ndx = ndx + 1
		}
	}
	ngstar = ndx

#### For Marc Davis -- print out grid of TV points for fitting:
#	call printf ("Mirr: Xarcs Yarcs, Xpix Ypix\n")
#	do i = -20, 220, 20 {
#	    xas = i
#	    do ndx = 85, 185, 20  {	# MIRR
#		yas = ndx
#		xtv = gseval (tvsf1x, xas, yas)
#		ytv = gseval (tvsf1y, xas, yas)
#		call printf ("     %7.2f %7.2f  %7.2f %7.2f\n")
#			call pargr (xas)
#			call pargr (yas)
#			call pargr (xtv)
#			call pargr (ytv)
#	    }
#	}
#
#	call printf ("Mask: Xarcs Yarcs, Xpix Ypix\n")
#	do i = -20, 220, 20 {
#	    xas = i
#	    do ndx = 170, 310, 20 {	# MASK
#		yas = ndx
#		xtv = gseval (tvsf2x, xas, yas)
#		ytv = gseval (tvsf2y, xas, yas)
#		call printf ("     %7.2f %7.2f  %7.2f %7.2f\n")
#			call pargr (xas)
#			call pargr (yas)
#			call pargr (xtv)
#			call pargr (ytv)
#	    }
#	}
#
##### End Marc Davis Section ###################

	call slit_sort (Memr[buftvm], Memi[buftvn], ngstar)

	if (std_fmt == YES) {
		call fprintf (fd, "\n# Non-Selected Objects:\n\n")
		do i = 0, ntarg-1 {
			if (SEL(tdat,i) == YES)
				next

# we now want to reconstitute the output line, changing only the SEL CODE
# (copied from above)
			stat = sscan (DATLINE(tdat,i))
			call gargwrd (ident, SZ_ID)
			call gargwrd (pars1, SZ_WRD)	# RA
			call gargwrd (pars2, SZ_WRD)	# Dec
			call gargwrd (pars3, SZ_WRD)	# Equinox
			call gargwrd (pars4, SZ_WRD)	# mag
			call gargwrd (pars5, SZ_WRD)	# pband
			call gargi (stat)	# skip Pcode
			call gargi (stat)	# skip Sample
			call gargi (stat)	# skip SelCode
			call gargstr (outline, SZ_LINE)		# copy rest
			if (nscan() < 10)
				call strcpy ("", outline, SZ_LINE)

			call fprintf (fd, "%-16s %12s %12s %s %5s %s %4d %1d %1d%s\n")
				call pargstr (ident)
				call pargstr (pars1)
				call pargstr (pars2)
				call pargstr (pars3)
				call pargstr (pars4)
				call pargstr (pars5)
				call pargi (PCODE(tdat,i))
				call pargi (SAMPL(tdat,i))
				call pargi (SEL(tdat,i))
				call pargstr (outline)

		}
	}

	if (strne (plotfile, "")) {
            fdp = open (plotfile, NEW_FILE, TEXT_FILE)
	    call write_mongo (fdp, tdat, ntarg, sdat, nslit, indat, fdtv1, fdtv2, Memi[buftvn], ngstar, delpa, plotfile)
	    call close (fdp)
	}

	call close (fdtv2)
	call close (fdtv1)

	call mfree (buftvm, TY_REAL)
	call mfree (buftvn, TY_INT)

end

#
# TEL2RADEC:  from tel axis and rotator PA, calc celest coords of any point
# Trial -- put in for Guider field center.  Seems to work OK to within a few
# pixels, but this whole section needs review.  In particular, the unref.
# Tel axis may not be quite right.
# At any rate, Using the DSS 2nd Gen red, apply the following on an eg 5x5' imag
# imlintran dss_im tv_im phi phi 0.203 0.203 ncol=1024 nline=1024
# where phi = 91.4 - (mask_pa)
#

procedure	tel2radec (indat, xas, yas, ra_pt, dec_pt)

pointer	indat
real	xas, yas
double	ra_pt, dec_pt

double	ra0_tel, dec0_tel
double	r, pa_pt, pa_fld

double	cosd, sind, cosa, sina, cost, sint, cosr, sinr

begin
# Find field center in unrefracted coords
# convert field center offset (arcsec) to radians
	r = DEGTORAD(sqrt (FLDCEN_X*FLDCEN_X + FLDCEN_Y*FLDCEN_Y) / 3600.)

# get PA of field center
	pa_fld = atan2 (FLDCEN_Y, FLDCEN_X)

	cosr = cos (r)
	sinr = sin (r)
	cosd = cos (DEC0_FLD(indat))
	sind = sin (DEC0_FLD(indat))

	cost = cos (PA_ROT(indat) - pa_fld)
	sint = sin (PA_ROT(indat) - pa_fld)

	sina = sinr * sint / cosd		# ASSUME not at dec=90
	cosa = sqrt (1. - sina*sina)

	ra0_tel = RA0_FLD(indat) - asin (sina)
	dec0_tel = asin ((sind*cosd*cosa - cosr*sinr*cost) /
				(cosr*cosd*cosa - sinr*sind*cost))
# convert field center offset (arcsec) to radians
	r = sqrt (xas**2 + yas**2) / 206204.8D0

# get PA of field center
	pa_pt = atan2 (yas, xas) + PI

	cosr = cos (r)
	sinr = sin (r)
	cosd = cos (dec0_tel)
	sind = sin (dec0_tel)

	cost = cos (PA_ROT(indat) - pa_pt)
	sint = sin (PA_ROT(indat) - pa_pt)

	sina = sinr * sint / cosd		# ASSUME not at dec=90
	cosa = sqrt (1. - sina*sina)

	ra_pt = ra0_tel - asin (sina)
	dec_pt = asin ((sind*cosd*cosa - cosr*sinr*cost) /
				(cosr*cosd*cosa - sinr*sind*cost))

end
