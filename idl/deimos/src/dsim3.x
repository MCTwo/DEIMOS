# Various subroutines needed in modifications

include	<math.h>
include	"dsimulator.h"
include	"keck.h"
include	"deimos.h"

#
# procedure	refr_coords (tdat, ntarg, indat)
# procedure	prec_coord (ra, dec, equinox, indat, ra, dec)
#


#
# REFR_COORDS: Calculate refracted coords (uses SLALIB).  Also get 
#	corresponding parallactic angle, Atm dispersion vectors and airmass
#	There is some ambiguity: does HA_FLD refer to refracted or unrefracted
#	coords?  (may need HA0_FLD) XXX
#

procedure	refr_coords (tdat, ntarg, indat)

pointer	tdat
int	ntarg
pointer	indat

int	i
double	lst			# LST corresponding to HA
double	lat			# observatory latitude (radians)
double	ha, az, el, zd		# HA, Az, El, ZDist (radians)

double	htm, pmb, tdk, rel_h20, w	# parameters for refraction calc
double	r1, r3				# refraction coeffs

double	w1, w2				# min and max wavelengths
double	a, b				# addn refractions coeffs
double	zd0, zd1, zd2			# addn zenith distii

double	slarms(), slpa()

begin
	lst = RA0_FLD(indat) + HA_FLD(indat)

	lat = DEGTORAD (OBS_LAT)	# radians
	htm = OBS_ALT			# meters
	tdk = TEMP(indat) + 273.15	# temp in K
	pmb = PRES(indat)		# millibars
	rel_h20 = OBS_RH		# relative humidity
	w = WAVER(indat)

# Get the refraction coeffs (2 hardcodes suggested in SLALIB):

	call slrfco (htm, tdk, pmb, rel_h20, w, lat, 0.0065D0, 1D-10, r1, r3)

# Save the refraction coeffs for later use:
	REF1(indat) = r1
	REF3(indat) = r3

# call eprintf ("Refraction Coeffs: %f %f\n")
# call pargd (r1)
# call pargd (r3)

# Apply to field center
	ha = lst - RA0_FLD(indat)
	call slde2h (ha, DEC0_FLD(indat), lat, az, el)
	zd0 = HALFPI - el
	call slrefz (zd0, r1, r3, zd)
	el = HALFPI - zd
	call sldh2e (az, el, lat, ha, DEC_FLD(indat))
	RA_FLD(indat) = lst - ha

# Loop and apply to targets:
	do i = 0, ntarg-1 {
		ha = lst - RA0(tdat,i)
		call slde2h (ha, DEC0(tdat,i), lat, az, el)

		zd = HALFPI - el
		call slrefz (zd, r1, r3, zd)
		el = HALFPI - zd
		
		call sldh2e (az, el, lat, ha, DEC(tdat,i))
		RA(tdat,i) = lst - ha
	}

# Now work out atmospheric dispersion:
	call slrefz (zd0, r1, r3, zd)

	w1 = WAVEMN(indat)
	w2 = WAVEMX(indat)
	call slatmd (tdk, pmb, rel_h20, w, r1, r3, w1, a, b)
	call slrefz (zd0, a, b, zd1)
	call slatmd (tdk, pmb, rel_h20, w, r1, r3, w2, a, b)
	call slrefz (zd0, a, b, zd2)

# ... amount
	AD1(indat) = (zd1 - zd) * 206205.
	AD2(indat) = (zd2 - zd) * 206205.

# ... and paralactic angle
	PAR_ANG(indat) = slpa (HA_FLD(indat), DEC_FLD(indat), lat)

# Finally, airmass:
	AMASS(indat) = slarms (zd)
call eprintf ("DEBUG: parang=%6.2f, AM=%5.3f\n")
call pargd (RADTODEG(PAR_ANG(indat)))
call pargd (AMASS(indat))

end

#
# UNREFR_COORDS: Calculate unrefracted coords (uses SLALIB).  Also get 
#	corresponding parallactic angle, Atm dispersion vectors and airmass
#	There is some ambiguity: does HA_FLD refer to refracted or unrefracted
#	coords?  (may need HA0_FLD) XXX
#	Ambiguity of lst below vs refr_coord
#

procedure	unrefr_coords (sdat, nslit, indat)

pointer	sdat
int	nslit
pointer	indat

int	i
double	lst			# LST corresponding to HA
double	lat			# observatory latitude (radians)
double	ha, az, el, zd		# HA, Az, El, ZDist (radians)
double	tanz

begin
	lst = RA_FLD(indat) + HA_FLD(indat)	# XXX Verify correct/see above
	lat = DEGTORAD (OBS_LAT)	# radians

# Apply to field center
	ha = lst - RA_FLD(indat)	## XXX Clean up (see above)
	call slde2h (ha, DEC_FLD(indat), lat, az, el)
	zd = HALFPI - el
	tanz = tan (zd)
	zd = zd + REF1(indat) * tanz + REF3(indat) * tanz**3
	el = HALFPI - zd
	call sldh2e (az, el, lat, ha, DEC0_FLD(indat))
	RA0_FLD(indat) = lst - ha
call eprintf ("Final Center:  %13.3h %12.2h\n")
call pargd (RADTODEG(RA0_FLD(indat))/15.d0)
call pargd (RADTODEG(DEC0_FLD(indat)))

call eprintf ("IMPT!! No unrefract to TEL -- FIX!! \n")

# Loop and apply to targets:
	do i = 0, nslit-1 {
		ha = lst - RA(sdat,i)
		call slde2h (ha, DEC(sdat,i), lat, az, el)

		zd = HALFPI - el
		tanz = tan (zd)
		zd = zd + REF1(indat) * tanz + REF3(indat) * tanz**3
		el = HALFPI - zd
		
		call sldh2e (az, el, lat, ha, DEC0(sdat,i))
		RA0(sdat,i) = lst - ha
	}


end


######################### CFITSIO routines ##########################
## NB: local "fitsio.h" comes from all the defines in cfitsio.h (?) with
## appropriate characters replaced.  XXX NOTE IN MAINTENANCE DOC

###  %	include "deimos$cfitsio/f77.inc"
include "fitsio.h"

procedure	marc (fname, tdat, ntarg, sdat, nslit, indat)

char	fname[ARB]		# FITS file name w/o extn
pointer	tdat
int	ntarg
pointer	sdat
int	nslit
pointer	indat

char	xname[SZ_FNAME]		# FITS name with extention
char	xline[80]		# XXX hardcode
%	character*80 f77nam
%	character*80 f77lin
%	character*80 f77nul

int	lu
int	stat
int	nc, nr

int	i, j

# double	dval			# general double value
# int	ival			# general int value
real	rval			# general real value

bool	isnul			# Null value(s) present?
double	dnul
int	inul
real	rnul

bool	streq()
int	sscan()
begin
	inul = INDEFI
	rnul = INDEFR
	dnul = INDEFD

# certain setup parameters needed
	ADJ_LEN(indat) = NO
	PROJ_LEN(indat) = NO

# Read in first table [1]: (DESNAME,DESAUTH,RA_PNT,DEC_PNT,PA_PNT,NSLIT)
call eprintf ("DEBUG: table[1]:\n")
	call strcpy (fname, xname, SZ_FNAME)
	call strcat ("[1]", xname, SZ_FNAME)
	call f77pak (xname, f77nam, SZ_FNAME)

	call ftgiou (lu, stat)
	call ftnopn (lu, f77nam, READONLY, stat)

	call ftgncl (lu, nc, stat)
	call ftgnrw (lu, nr, stat)

# check format:
	if (nc != 6 || nr != 1) {
		call eprintf ("%s, nc=%d nr=%d stat=%d\n")
		call pargstr (xname)
		call pargi (nc)
		call pargi (nr)
		call pargi (stat)
		call fatal (0, "Unexpected format in input FITS table")
	}

	call f77pak ("???", f77nul, 80)
	call ftgcvs (lu, 1, 1, 1, 1, f77nul, f77lin, isnul, stat)
	call f77upk (f77lin, DESNAME(indat), 80)
	call ftgcvs (lu, 2, 1, 1, 1, f77nul, f77lin, isnul, stat)
	call f77upk (f77lin, PROJNAME(indat), 80)

	call strcpy ("marc davis <marc@astron.Berkeley.EDU>", DESAUTH(indat), SZ_LINE)
	call strcpy (DESNAME(indat), BLUNAME(indat), SZ_LINE)
	call strcpy ("Marc Davis (UCB) via Dsim", DESCREAT(indat), SZ_LINE)

	call ftgcvd (lu, 3, 1, 1, 1, dnul, RA0_FLD(indat), isnul, stat)
	call ftgcvd (lu, 4, 1, 1, 1, dnul, DEC0_FLD(indat), isnul, stat)
	call ftgcvd (lu, 5, 1, 1, 1, dnul, PA_ROT(indat), isnul, stat)
	call ftgcvj (lu, 6, 1, 1, 1, inul, nslit, isnul, stat)
	RA0_FLD(indat) = DEGTORAD(RA0_FLD(indat))
	DEC0_FLD(indat) = DEGTORAD(DEC0_FLD(indat))
	PA_ROT(indat) = DEGTORAD(PA_ROT(indat))

	call ftclos (lu, stat)

# Read in second table [2]: (OBJID,RA,DEC,Mag,PBand,Majax,PA,minax,objclass)
call eprintf ("DEBUG: table[2]:\n")
	call strcpy (fname, xname, SZ_FNAME)
	call strcat ("[2]", xname, SZ_FNAME)
	call f77pak (xname, f77nam, SZ_FNAME)

	call ftgiou (lu, stat)
	call ftnopn (lu, f77nam, READONLY, stat)

	call ftgncl (lu, nc, stat)
	call ftgnrw (lu, nr, stat)
# check format:
	if (nc != 9 || nr < nslit ) {
		call eprintf ("%s, nc=%d nr=%d stat=%d\n")
		call pargstr (xname)
		call pargi (nc)
		call pargi (nr)
		call pargi (stat)
		call fatal (0, "Unexpected format in input FITS table")
	}
	ntarg = nr
	call targ_alloc (tdat, ntarg)

# Load the vectors:
	call f77pak ("???", f77nul, 80)
#                           (OBJID,RA,DEC,Mag,PBand,Majax,PA,minax,objclass)
	do i = 0, ntarg-1 {
		j = i + 1
# Get the ID
		INDEX(tdat,i) = i
		call ftgcvs (lu, 1, j, 1, 1, f77nul, f77lin, isnul, stat)
		call f77upk (f77lin, DATLINE(tdat,i), 80)

		call ftgcvd (lu, 2, j, 1, 1, dnul, RA0(tdat,i), isnul, stat)
		call ftgcvd (lu, 3, j, 1, 1, dnul, DEC0(tdat,i), isnul, stat)
		call ftgcve (lu, 4, j, 1, 1, rnul, MAG(tdat,i), isnul, stat)
		call ftgcvs (lu, 5, j, 1, 1, f77nul, f77lin, isnul, stat)
		call f77upk (f77lin, xline, 1)
		PBAND(tdat,i) = xline[1]
		if (xline[2] != EOS)
			call eprintf ("Warning; pband truncated\n")
# XXX Majax
		call ftgcve (lu, 7, j, 1, 1, rnul, PA(tdat,i), isnul, stat)
		call ftgcvs (lu, 9, j, 1, 1, f77nul, f77lin, isnul, stat)
		call f77upk (f77lin, xline, 1)
		if (streq (xline, "P"))
			PCODE(tdat,i) = 1
		else if (streq (xline, "G"))
			PCODE(tdat,i) = CODE_GS
		else if (streq (xline, "A"))
			PCODE(tdat,i) = CODE_AS
		else
		    call eprintf ("Unknown ObjClass, line %d\n"); call pargi (i)

# XXX  MinAx objclass
		RA0(tdat,i)  = DEGTORAD(RA0(tdat,i))
		DEC0(tdat,i) = DEGTORAD(DEC0(tdat,i))
		PA(tdat,i)   = DEGTORAD(PA(tdat,i))
	}
# In addition, there are several other vectors; for trial XXX:
	call amovki (YES, SEL(tdat,0),  ntarg)
	call amovki (1,   SAMPL(tdat,0), ntarg)

	call ftclos (lu, stat)

# Read in first table [3]: (DSLITID,DESID,SLITNAME,RA,DEC,TYP,LEN,PA,WID, (wpa))
call eprintf ("DEBUG: table[3]:\n")
	call strcpy (fname, xname, SZ_FNAME)
	call strcat ("[3]", xname, SZ_FNAME)
	call f77pak (xname, f77nam, SZ_FNAME)

	call ftgiou (lu, stat)
	call ftnopn (lu, f77nam, READONLY, stat)

	call ftgncl (lu, nc, stat)
	call ftgnrw (lu, nr, stat)
# check format:
	if (nc != 10 || nr < 1 || nr > ntarg) {
		call eprintf ("%s, nc=%d nr=%d stat=%d\n")
		call pargstr (xname)
		call pargi (nc)
		call pargi (nr)
		call pargi (stat)
		call fatal (0, "Unexpected format in input FITS table")
	}
	nslit = nr
	call slit_alloc (sdat, nslit)

# Load the vectors:
#                           (DSLITID,DESID,SLITNAME,RA,DEC,TYP,LEN,PA,WID)
	call strcpy ("", xline, 80)
	call f77pak ("-9", f77nul, 80)
	do i = 0, nslit-1 {
		j = i + 1
		call ftgcvs (lu, 3, j, 1, 1, f77nul, f77lin, isnul, stat)
		call f77upk (f77lin, xline, 80)
		if (sscan (xline) != EOS) {
			call gargi (INDEX(sdat,i))
			if (INDEX(sdat,i) != i)
				call eprintf ("Sequencing error in Slitname!\n")
		} else {
			INDEX(sdat,i) = nslit		# XXX BOGUS!!
		}

# XXX IDs and lots else
		call ftgcvd (lu, 4, j, 1, 1, dnul, RA0(sdat,i), isnul, stat)
		call ftgcvd (lu, 5, j, 1, 1, dnul, DEC0(sdat,i), isnul, stat)
		call ftgcvs (lu, 6, j, 1, 1, f77nul, f77lin, isnul, stat)
		call f77upk (f77lin, xline, 1)
		if (streq (xline, "P"))
			PCODE(sdat,i) = 1
		else if (streq (xline, "G"))
			PCODE(sdat,i) = CODE_GS
		else if (streq (xline, "A"))
			PCODE(sdat,i) = CODE_AS
		else
		    call eprintf ("Unknown SlitTyp, line %d\n"); call pargi (i)

		call ftgcve (lu, 8, j, 1, 1, rnul, PA(sdat,i), isnul, stat)
		RA0(sdat,i)  = DEGTORAD(RA0(sdat,i))
		DEC0(sdat,i) = DEGTORAD(DEC0(sdat,i))
		PA(sdat,i)   = DEGTORAD(PA(sdat,i))
		call ftgcve (lu, 9, j, 1, 1, rnul, SLWID(sdat,i), isnul, stat)
		call ftgcve (lu, 7, j, 1, 1, rnul, rval, isnul, stat)
		LEN1(sdat,i) = 0.5 * rval
		LEN2(sdat,i) = 0.5 * rval
	}
	call aaddkr (RELPA(sdat,0), -1.*PA_ROT(indat), RELPA(sdat,0), nslit)

	call ftclos (lu, stat)


# Read in fourth table [4]: ()
call eprintf ("DEBUG: table[4]:\n")
	call strcpy (fname, xname, SZ_FNAME)
	call strcat ("[4]", xname, SZ_FNAME)
	call f77pak (xname, f77nam, SZ_FNAME)

	call ftgiou (lu, stat)
	call ftnopn (lu, f77nam, READONLY, stat)

	call ftgncl (lu, nc, stat)
	call ftgnrw (lu, nr, stat)
# check format:
	if (nc != 6 || nr > ntarg || nr < nslit) {
		call eprintf ("%s, nc=%d nr=%d stat=%d\n")
		call pargstr (xname)
		call pargi (nc)
		call pargi (nr)
		call pargi (stat)
		call fatal (0, "Unexpected format in input FITS table")
	}

# Load the vectors:
	call amovki (-1, SLNDX(tdat,0), ntarg)	# fill with nulls

	call strcpy ("", xline, 80)
	call f77pak ("-9", f77nul, 80)
	do i = 0, nr-1 {			# XXX resolve length
		j = i + 1
		call ftgcvs (lu, 6, j, 1, 1, f77nul, f77lin, isnul, stat)
		call f77upk (f77lin, xline, 80)
		if (sscan (xline) != EOS)
			call gargi (SLNDX(tdat,i))
		else
			SLNDX(tdat,i) = -9
	}

	call ftclos (lu, stat)

end

# MARC2: make some calcs that are not made elsewhere

procedure	marc2 (sdat, nslit, tdat, ntarg)

pointer	sdat
int	nslit

pointer	tdat
int	ntarg

int	i, j
real	sina, cosa

begin
# XXX Code lifted from mask_coord, but covers stuff in gen_slit
	do i = 0, nslit-1 {
# XXX For now, carry through the RELPA thing; in end, must be specified!
## XXX RESOLVE FLIP issue here, mark_obj, gen_slits !!!
			if (RELPA(sdat,i) != INDEF) {
				cosa = cos (RELPA(sdat,i))
				sina = sin (RELPA(sdat,i))
			} else {
				cosa = 1.
				sina = 0.
			}

			X1(sdat,i) = XARCS(sdat,i) - LEN1(sdat,i) * cosa
			Y1(sdat,i) = YARCS(sdat,i) - LEN1(sdat,i) * sina * FLIP
			X2(sdat,i) = XARCS(sdat,i) + LEN2(sdat,i) * cosa
			Y2(sdat,i) = YARCS(sdat,i) + LEN2(sdat,i) * sina * FLIP
	}

# Assign some values to the tdat structure
	do i = 0, ntarg-1 {
		j = SLNDX(tdat,i)
		if (j >= 0) {
			SLWID(tdat,i) = SLWID(sdat,j)
#			X1(tdat,i) = X1(sdat,j)
#			X2(tdat,i) = X2(sdat,j)
		}
	}
end


##### RULES FOR MARC:
# 1. Slit index shall start at 0 and slits shall be in order
# 2. columns shall not change order; tables shall not change order

##### Where I am going:
# 1. need to separate GEN_SLIT into allocation and filling parts
# 2. For marc:  allocate slits and targets; fill appropriately.  Note that
#	slits will have to work backwards, ie, we have RA,dec,etc, must
#	convert to arcsec on sky.


#### PROBLEM NOTES:
# with the FITSIO routines, I suspect they are read sequentially.  Thus
# accessing another forward part may cause problems.  If weird format, etc
# occurs, make sure all calls have defined arguments.

# In my FITS-table routines, overfilling a field will cause "No permission
# on String File" error.  This is because IRAF will expand a field to fit
# numbers that are too large.
