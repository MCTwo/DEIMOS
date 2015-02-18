# TBD:  get various params into proper routines

include	<math.h>
# include	<gset.h>
# include	<gim.h>
include	<time.h>
include	"dsimulator.h"
include	"deimos.h"

#
# procedure	data_init (indat)
# procedure	targ_init (fd, tdat, ntarg, indat)
# procedure	targ_alloc (tdat, ntarg)
# procedure	gen_slits (tdat, ntarg, sdat, nslit, indat)
# procedure	slit_alloc (sdat, nslit)
# procedure	tel_coords (tdat, ntarg, indat)
# int procedure	chk_stat (x, y, full_check)


#
# DATA_INIT: initialize data structure for telescope/background data
# 

procedure	data_init (indat)

pointer	indat

int	date_decode()
int	name_check()

double	clgetd ()
real	clgetr()
bool	clgetb(), streq()
begin

# Allocate the data structure and vectors
	call malloc (indat, NINDAT, TY_STRUCT)
	call malloc (PTTELDAT(indat), NTELPAR, TY_DOUBLE)
	call malloc (PTDEFDAT(indat), NDEFPAR, TY_REAL)
	call malloc (PTMSKDAT(indat), NMSKPAR*SZ_LINE, TY_CHAR)

# Read in params
	RA0_FLD(indat)  = DEGTORAD(15. * clgetd ("ra0"))
	DEC0_FLD(indat) = DEGTORAD(clgetd ("dec0"))
	HA_FLD(indat) = DEGTORAD(15. * clgetr ("ha0"))
	PA_ROT(indat) = DEGTORAD(clgetd ("PA0"))
	TEMP(indat) = clgetr ("temp")
	PRES(indat) = clgetr ("pressure")
	WAVER(indat) = clgetr ("lambda_cen") * 1.e-4		# in microns
	WAVEMN(indat) = clgetr ("blue") * 1.e-4
	WAVEMX(indat) = clgetr ("red") * 1.e-4
	SLIT_GAP(indat) = clgetr ("sep_slit")		# sep. bet. slits, asec
	DEF_HLEN(indat) = 0.5 * clgetr ("min_slit")	# min. length in arcsec
	DEF_BOXR(indat) = 0.5 * clgetr ("box_sz")	# box 1/2-length in arcs
	DEF_SLWID(indat) = clgetr ("slit_width")	# slit width in arcsec
	STD_EQX(indat) = clgetr ("equinox")

	if (clgetb ("proj_len"))
		PROJ_LEN(indat) = YES
	else
		PROJ_LEN(indat) = NO

	if (clgetb ("no_overlap"))
		ADJ_LEN(indat) = YES
	else
		ADJ_LEN(indat) = NO

	call clgstr ("maskid",   DESNAME(indat), SZ_LINE)
	call clgstr ("maskid",   BLUNAME(indat), SZ_LINE)
	call clgstr ("guiname",  GUINAME(indat), SZ_LINE)
	call clgstr ("observer", BLUOBSR(indat), SZ_LINE)
	call clgstr ("author",   DESAUTH(indat), SZ_LINE)
	call clgstr ("project",  PROJNAME(indat), SZ_LINE)
	call clgstr ("dateobs",  USEDATE(indat), SZ_LINE)
	call strcpy ("DEIMOS", INSTRUME(indat), SZ_LINE)
	call strcpy ("Dsimulator: Ver 0.0b", DESCREAT(indat), SZ_LINE)

	if (date_decode (USEDATE(indat), EPOCH(indat)) != OK)
		call fatal (0, "USE_DATE format improper!")

	if (streq (GUINAME(indat), ""))
		call fatal (0, "MUST Specify a GUI name!")

	if (name_check (BLUOBSR(indat)) != OK)
	    call fatal (0, "'observer' format not 'name <user@domain>'")
	if (name_check (DESAUTH(indat)) != OK)
	    call fatal (0, "'author' format not 'name <user@domain>'")
	

end


#
# TARG_INIT: initialize data structure for targets; fill
# 

procedure	targ_init (fd, tdat, ntarg, indat)

pointer	fd
pointer	tdat
int	ntarg
pointer	indat

char	tchar
char	idstr[SZ_ID]
char	workstr[SZ_ID], pastr[3]	# for identifying "field" data
char	passband[SZ_PBAND]		# passband
int	prior, nlist, selcode
real	eqnx_std		# TMP, needs input
real	equinox, pangle, l1, l2
real	magn
double	alpha, delta

int	ndx, n

int	line_count()
bool	streq()
int	fscan(), sscan(), nscan(), strlen()

begin
# temp
	eqnx_std = STD_EQX(indat)

# Count the entries
	ntarg = line_count (fd)

# Allocate the vectors
	call targ_alloc (tdat, ntarg)

# Read in data
	ndx = 0
	while (fscan (fd) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#' || nscan() == 0) {
			next
		}
		call reset_scan ()

		if (ndx == 0) {		# check for mask center
			call gargwrd (idstr, SZ_ID)
			call gargd (alpha)
			call gargd (delta)
			call gargr (equinox)
			call gargwrd (workstr, SZ_ID)
			call strcpy (workstr, pastr, 3)
			if (nscan() >= 5 && streq (pastr, "PA=")) {
				call eprintf ("Field data from input file...\n")
				RA0_FLD(indat) = DEGTORAD(15. * alpha)
				DEC0_FLD(indat) = DEGTORAD(delta)
				STD_EQX(indat) = equinox
				eqnx_std = STD_EQX(indat)	# TMP!
				if (sscan (workstr[4]) != EOS) {
				    call gargr (pangle)
				    if (nscan() != 0) {
					PA_ROT(indat) = DEGTORAD(pangle)
					next
				    } else {
					call fatal (0, "Bad format on PA=")
				    }
				}
			}
			call reset_scan ()
		}

		call gargwrd (idstr, SZ_ID)
		call gargd (alpha)
		call gargd (delta)
		call gargr (equinox)
		call gargr (magn)
		call gargwrd (passband, SZ_PBAND)
		call gargi (prior)

		if (nscan() < 7) {
			call eprintf ("Bad format on input line -- skipped\n`%s'\n")
			next
		}

		call gargi (nlist)
		call gargi (selcode)

		call gargr (pangle)
		call gargr (l1)
		call gargr (l2)

# Supply default arguments as needed:
		n = nscan()

# Are lengths present (either both or neither)
		if (n < 12) {
			l1 = DEF_HLEN(indat)
			l2 = DEF_HLEN(indat)
		}
		if (n < 10)
			pangle = INDEF
		if (n < 9)
			selcode = NO
		if (n < 8)
			nlist = PRIMARY

# Check for some special cases:
# XXX should check that nlist = reasonable value

		if (prior == CODE_AS) {
			l1 = DEF_BOXR(indat)
			l2 = DEF_BOXR(indat)
		}

		call reset_scan ()
		call gargstr (DATLINE(tdat,ndx), SZ_LINE-1)

# ASSUME no proper motion updates needed, for now, but add here in future
# XXX  will need epoch, epoch_std

# put on standard equinox, if needed
	if (equinox != eqnx_std) {
		call eprintf ("EQUINOX Conflict: %7.2f vs standard = %7.2f\n")
			call pargr (equinox)
			call pargr (eqnx_std)
		call fatal (0, "No provision for precession yet!")
		equinox = eqnx_std	# XXX
	}


# Assign values, always in radians for angles; arcsec for short lengths (???)
		INDEX(tdat,ndx) = ndx
		RA0(tdat,ndx) = DEGTORAD(alpha*15.)
		DEC0(tdat,ndx) = DEGTORAD(delta)

		MAG(tdat,ndx) = magn
		if (strlen (passband) > 1)
			call eprintf ("WARNING: passband truncated to 1 char\n")
		PBAND(tdat,ndx) = passband[1]

		PCODE(tdat,ndx) = prior
		SAMPL(tdat,ndx) = nlist
		SEL(tdat,ndx) = selcode

		if (pangle != INDEF)
			PA(tdat,ndx) = DEGTORAD(pangle)
		else
			PA(tdat,ndx) = INDEF

		LEN1(tdat,ndx) = l1
		LEN2(tdat,ndx) = l2
			
# XXX Assign slit-width
		if (PCODE(tdat,ndx) == CODE_AS) {
			PA(tdat,ndx) = INDEF
			SLWID(tdat,ndx) = 2.*DEF_BOXR(indat)
		} else {
			SLWID(tdat,ndx) = DEF_SLWID(indat)
		}

		ndx = ndx + 1
	}
	ntarg = ndx

end

#
# TARG_ALLOC: allocate arrays for targets (broken out for MD ingestor)
#
procedure	targ_alloc (tdat, ntarg)

pointer	tdat
int	ntarg

begin
# Allocate the vectors
	call malloc (tdat, TDATLEN, TY_STRUCT)
	call malloc (PTINDEX(tdat), ntarg, TY_INT)
	call malloc (PTRA0(tdat), ntarg, TY_DOUBLE)
	call malloc (PTDEC0(tdat), ntarg, TY_DOUBLE)
	call malloc (PTRA(tdat), ntarg, TY_DOUBLE)
	call malloc (PTDEC(tdat), ntarg, TY_DOUBLE)
	call malloc (PTPA(tdat), ntarg, TY_REAL)
	call malloc (PTLEN1(tdat), ntarg, TY_REAL)
	call malloc (PTLEN2(tdat), ntarg, TY_REAL)
	call malloc (PTWID(tdat), ntarg, TY_REAL)
	call malloc (PTPCODE(tdat), ntarg, TY_INT)
	call malloc (PTSAMPL(tdat), ntarg, TY_INT)

	call malloc (PTSTAT(tdat), ntarg, TY_INT)
	call malloc (PTSEL(tdat), ntarg, TY_INT)
	call malloc (PTXARCS(tdat), ntarg, TY_REAL)
	call malloc (PTYARCS(tdat), ntarg, TY_REAL)
	call malloc (PTRELPA(tdat), ntarg, TY_REAL)
	call malloc (PTX1(tdat), ntarg, TY_REAL)
	call malloc (PTY1(tdat), ntarg, TY_REAL)
	call malloc (PTX2(tdat), ntarg, TY_REAL)
	call malloc (PTY2(tdat), ntarg, TY_REAL)
	call malloc (PTMAG(tdat), ntarg, TY_REAL)
	call malloc (PTPBAND(tdat), ntarg, TY_CHAR)	## XXX should be general

# Allocate slit-index -- cross-reference -- NOTE the zeroing of this array
	call calloc (PTSLNDX(tdat), ntarg, TY_INT)

# Allocate the string array; fill w/ End-of-string
	call malloc (PTLINE(tdat), ntarg*SZ_LINE, TY_CHAR)
	call amovkc (EOS, DATLINE(tdat,0), ntarg*SZ_LINE)

end



#
# TEL_COORDS: Convert (refracted) alpha,dec into offsets from telescope center.
#

procedure	tel_coords (tdat, ntarg, indat)

pointer	tdat
int	ntarg
pointer	indat

int	i
double	r			# radius of object from tel-axis
double	p			# PA on sky from tel-axis
double	dec_obj, del_ra
double	cosr, sinp, cosp
double	ra0, dec0, pa0		# RA, Dec and PA on axis

real	rangle, xgeom, ygeom
real	x, y

int	chk_stat()
begin
	ra0  = RA_TEL(indat)
	dec0 = DEC_TEL(indat)
	pa0  = PA_ROT(indat)

	do i = 0, ntarg-1 {
		dec_obj = DEC(tdat,i)
		del_ra = RA(tdat,i) - ra0
		cosr = sin (dec_obj) * sin (dec0) +
			cos (dec_obj) * cos (dec0) * cos (del_ra)
		r = acos (cosr)

		sinp = cos (dec_obj) * sin (del_ra) / sqrt (1. - cosr*cosr)
		cosp = sqrt (max ((1. - sinp*sinp), 0.))
		if (dec_obj < dec0)
			cosp = -cosp
		p = atan2 (sinp, cosp)

# For now, convert radii to arcsec XXX
# XXX NB: I am not sure this is correct!  We should still be on SPHERICAL surf.
# 	but these are EUCLIDEAN relations.  Options: work in spherical coord
#	OR work in tan projection.
# More:  The difference at 10 arcmin between the tan and angle is < 0.002 arcsec
# If we convert "r" to tan(r) we should have the tan projection
		r = tan(r) * 206264.8
#		r = RADTODEG(r) * 3600.
		XARCS(tdat,i) = r * cos (pa0 - p)
		YARCS(tdat,i) = r * sin (pa0 - p)

		if (PA(tdat,i) == INDEF) {
			RELPA(tdat,i) = INDEF
			rangle = 0.
		} else {
			RELPA(tdat,i) = PA(tdat,i) - pa0
			rangle = RELPA(tdat,i)
		}


# For simplicity, we calculate the endpoints in X here; note use of FLIP
		xgeom = FLIP * cos (rangle)
		ygeom = sin (rangle)
		if (PROJ_LEN(indat) == YES) {
			xgeom = xgeom / abs (cos (rangle))
			ygeom = ygeom / abs (cos (rangle))
		}
# We always want X1 < X2, so:
		if (xgeom > 0) {
			X1(tdat,i) = XARCS(tdat,i) - LEN1(tdat,i) * xgeom
			Y1(tdat,i) = YARCS(tdat,i) - LEN1(tdat,i) * ygeom
			X2(tdat,i) = XARCS(tdat,i) + LEN2(tdat,i) * xgeom
			Y2(tdat,i) = YARCS(tdat,i) + LEN2(tdat,i) * ygeom
		} else {
			X2(tdat,i) = XARCS(tdat,i) - LEN1(tdat,i) * xgeom
			Y2(tdat,i) = YARCS(tdat,i) - LEN1(tdat,i) * ygeom
			X1(tdat,i) = XARCS(tdat,i) + LEN2(tdat,i) * xgeom
			Y1(tdat,i) = YARCS(tdat,i) + LEN2(tdat,i) * ygeom
		}

# Calc STAT
		x = XARCS(tdat,i)	# unclear if XYARCS will be real or dbl
		y = YARCS(tdat,i)
		STAT(tdat,i) = chk_stat (x, y, YES)
	}
end

#
## CHK_STAT: is object within the slitmask? REAL INPUTS; INSTRUMENT SPECIFIC
## This could also be replaced by a bunch of limiting curves
## Note that some slit info should be passed along, also -- width, tilt, ends
#

int procedure	chk_stat (x, y, full_check)

real	x, y
int	full_check

real	r

begin
	r = sqrt (x*x + y*y)

# Is object within 10 arcmin radius?
	if (r > 600.)
		return (NO)

# inner edge of mask
	if (y < YMSKMIN)
		return (NO)

# outer edge of mask
	if (y > YMSKMAX)
		return (NO)

# outer edge of mask
	if (x > XUPP_LIM || x < XLOW_LIM)
		return (NO)

# cut corner
	if (x > -0.98273 * y + 833.0)
		return (NO)

	if (full_check == NO)
		return (YES)		# OK to put slit there

# within radius of camera obscuration/vignetting?
	if (x*x+(y-YCAMCEN)**2 < RADVIGN**2)
		return (NO)
	

# near gaps in mosaic?
# XXX needs definition;  things like this should be contained in defines and
# limits (how close) in parameters
#	if (abs (x+250.) < 4. || abs (x) < 4. || abs(x-250.) < 4.)

	if (abs (x-GAP1CEN) < GAP1HWD)
		return (NO)
	if (abs (x-GAP2CEN) < GAP2HWD)
		return (NO)
	if (abs (x-GAP3CEN) < GAP3HWD)
		return (NO)

# appears OK...
	return (YES)
end

	


#
# GEN_SLITS: initialize data structure for slits; fill (generate slits)
# 

procedure	gen_slits (tdat, ntarg, sdat, nslit, indat)

pointer	tdat
int	ntarg
pointer	sdat
int	nslit
pointer	indat

int	ndx, i
real	x, y

int	chk_stat()
begin

# Count the selected targets
	nslit = 0
	do i = 0, ntarg-1 {
		if (SEL(tdat,i) == YES)
			nslit = nslit + 1
	}

# Allocate the vectors:
	call slit_alloc (sdat, nslit)		# nslit includes GS, etc

# Set up slits for selected targets:
	ndx = 0
	do i = 0, ntarg-1 {
		if (PCODE(tdat,i) == CODE_GS)	# Ignore guide stars
			next
		if (SEL(tdat,i) == YES)	{	# or != 0
			x = XARCS(tdat,i)	# unclear TY of XYARCS
			y = YARCS(tdat,i)
			if (chk_stat (x, y, NO) == NO)
				next		# Not on metal

			INDEX(sdat,ndx) = ndx 
			if (PA(tdat,i) == INDEF) {
				PA(sdat,ndx) = PA_ROT(indat)
			} else {
				PA(sdat,ndx) = PA(tdat,i)
			}
			RELPA(sdat,ndx) = RELPA(tdat,i)
			PCODE(sdat,ndx) = PCODE(tdat,i)

			X1(sdat,ndx) = X1(tdat,i)
			Y1(sdat,ndx) = Y1(tdat,i)
			X2(sdat,ndx) = X2(tdat,i)
			Y2(sdat,ndx) = Y2(tdat,i)

# XXX NB: until the final sky_coords are calc'd, want X/YARCS to repr. objects
			XARCS(sdat,ndx) = XARCS(tdat,i)
			YARCS(sdat,ndx) = YARCS(tdat,i)
# XXX cuidado!  I am not sure that the tan-projection of the rel PA is the
# same as the rel PA -- MUST CHECK!

			SLWID(sdat,ndx) = SLWID(tdat,i)

# This is where we also assign slit index to object
			SLNDX(tdat,i) = ndx

			ndx = ndx + 1
		}
	}
	nslit = ndx

	if (ADJ_LEN(indat) == YES)
		call len_slits (tdat, ntarg, sdat, nslit, indat)
			
end

#
# SLIT_ALLOC: allocate arrays for targets (broken out for MD ingestor)
# SLIT_FREE: free memory
#

procedure	slit_alloc (sdat, nslit)

pointer	sdat
int	nslit

begin
# Allocate the vectors: note that many quantities match the tdat structure
	call malloc (sdat, SDATLEN, TY_STRUCT)
	call malloc (PTINDEX(sdat), nslit, TY_INT)
	call malloc (PTRA0(sdat), nslit, TY_DOUBLE)
	call malloc (PTDEC0(sdat), nslit, TY_DOUBLE)
	call malloc (PTRA(sdat), nslit, TY_DOUBLE)
	call malloc (PTDEC(sdat), nslit, TY_DOUBLE)
	call malloc (PTPA(sdat), nslit, TY_REAL)
	call malloc (PTLEN1(sdat), nslit, TY_REAL)
	call malloc (PTLEN2(sdat), nslit, TY_REAL)
	call malloc (PTWID(sdat), nslit, TY_REAL)
	call malloc (PTPCODE(sdat), nslit, TY_INT)
	call malloc (PTXARCS(sdat), nslit, TY_REAL)
	call malloc (PTYARCS(sdat), nslit, TY_REAL)
	call malloc (PTRELPA(sdat), nslit, TY_REAL)
	call malloc (PTX1(sdat), nslit, TY_REAL)
	call malloc (PTY1(sdat), nslit, TY_REAL)
	call malloc (PTX2(sdat), nslit, TY_REAL)
	call malloc (PTY2(sdat), nslit, TY_REAL)
	call malloc (PTSTAT(sdat), nslit, TY_INT)

	call malloc (PTSCOOR(sdat), nslit*NSCOOR, TY_DOUBLE)
	return

entry	slit_free (sdat, nslit)

	call mfree (PTINDEX(sdat), TY_INT)
	call mfree (PTRA0(sdat), TY_DOUBLE)
	call mfree (PTDEC0(sdat), TY_DOUBLE)
	call mfree (PTRA(sdat), TY_DOUBLE)
	call mfree (PTDEC(sdat), TY_DOUBLE)
	call mfree (PTPA(sdat), TY_REAL)
	call mfree (PTLEN1(sdat), TY_REAL)
	call mfree (PTLEN2(sdat), TY_REAL)
	call mfree (PTWID(sdat), TY_REAL)
	call mfree (PTPCODE(sdat), TY_INT)
	call mfree (PTXARCS(sdat), TY_REAL)
	call mfree (PTYARCS(sdat), TY_REAL)
	call mfree (PTRELPA(sdat), TY_REAL)
	call mfree (PTX1(sdat), TY_REAL)
	call mfree (PTY1(sdat), TY_REAL)
	call mfree (PTX2(sdat), TY_REAL)
	call mfree (PTY2(sdat), TY_REAL)
	call mfree (PTSTAT(sdat), TY_INT)
	call mfree (PTSCOOR(sdat), TY_DOUBLE)

	call mfree (sdat, TY_STRUCT)
	sdat = 0
	nslit = 0
end


#
# BPSTD: Standard Boilerplate
# BPADD: Additional Boilerplate
#

procedure	bpstd (fd, kdat, tabname, npk, nfk, nam1, nam2, nam3)

pointer	fd		# file descriptor
pointer	kdat		# 
char	tabname[ARB]	# tablename
int	npk		# no. of primary keys
int	nfk		# no. of foreign keys
char	nam1[ARB], nam2[ARB], nam3[ARB]

# char	pnam1[ARB]	# first primary name
# char	pnam2[ARB]	# second primary name

char	date[SZ_TIME]					# time/date of file

int	nkey

char	snul
char	creatask[60]
char	person[30]

begin
# XXX
	call strcpy ("", snul, 1)
	call strcpy ("Phillips <phillips@ucolick.org>", person, 35)
	call strcpy ("DSIMULATOR -- June2000", creatask, 35)

	call ft_date (date, SZ_TIME)
# XXX

	call pkwstr (fd, kdat, "EXTNAME", tabname, snul)
	call pkwi   (fd, kdat, "EXTVER", 0, snul)
	call pkwstr (fd, kdat, "DATE", date, snul)
	call pkwstr (fd, kdat, "AUTHOR", person, snul)
	call pkwstr (fd, kdat, "CREATOR", creatask, snul)

	call pkwi   (fd, kdat, "NPRIKEY", npk, snul)

	if (nfk > 0) {
		call pkwi   (fd, kdat, "NFORKEY", nfk, snul)
		call pkwi   (fd, kdat, "NFORHDU", nfk, snul)
	}

	call pkwstr (fd, kdat, "PKTYP1", nam1, snul)
	if (npk == 2)
		call pkwstr (fd, kdat, "PKTYP2", nam2, snul)

	return

entry	bpadd (fd, kdat, nkey, nam1, nam2, nam3)

	if (nkey == 1) {
		call pkwstr (fd, kdat, "FKTYP1", nam1, snul)
		call pkwstr (fd, kdat, "FKFYP1", nam2, snul)
		call pkwi   (fd, kdat, "FKHDU1", nkey, snul)
		call pkwstr (fd, kdat, "HDLOC1", snul, snul)
		call pkwstr (fd, kdat, "HDXTN1", "TABLE", snul)
		call pkwstr (fd, kdat, "HDNAM1", nam3, snul)
	} else if (nkey == 2) {
		call pkwstr (fd, kdat, "FKTYP2", nam1, snul)
		call pkwstr (fd, kdat, "FKFYP2", nam2, snul)
		call pkwi   (fd, kdat, "FKHDU2", nkey, snul)
		call pkwstr (fd, kdat, "HDLOC2", snul, snul)
		call pkwstr (fd, kdat, "HDXTN2", "TABLE", snul)
		call pkwstr (fd, kdat, "HDNAM2", nam3, snul)
	} else if (nkey == 3) {
		call pkwstr (fd, kdat, "FKTYP3", nam1, snul)
		call pkwstr (fd, kdat, "FKFYP3", nam2, snul)
		call pkwi   (fd, kdat, "FKHDU3", nkey, snul)
		call pkwstr (fd, kdat, "HDLOC3", snul, snul)
		call pkwstr (fd, kdat, "HDXTN3", "TABLE", snul)
		call pkwstr (fd, kdat, "HDNAM3", nam3, snul)
	}

	return
end


#
# DATE_DECODE: decode use_date; NB only approx, yet stored in double!!
#
define	SZ_DATE_STR	10

int	procedure date_decode (date_str, epoch)

char	date_str[ARB]		# date-string
double	epoch

char	tchar, rchar
char	wkstr[SZ_DATE_STR]
int	i
int	yyyy, mm, dd

int	stridx(), sscan(), nscan()
begin

	tchar = '-'
	rchar = ' '

	i = 1			# Should actually search for non-blank
	call strcpy (date_str[i], wkstr, SZ_DATE_STR)

	i = stridx (tchar, wkstr)
	wkstr[i] = rchar
	if (i <= 1)
		return (ERR)

	i = stridx (tchar, wkstr)
	wkstr[i] = rchar
	if (i <= 1)
		return (ERR)

	i = sscan (wkstr)
		call gargi (yyyy)
		call gargi (mm)
		call gargi (dd)
	if (nscan() < 3) {
		epoch = INDEF
		return (ERR)
	} else {
		epoch = yyyy + mm/12. + dd/365.		#  Approx; XXX
		return (OK)
	}
end

#
# NAME_CHECK: check name format to ascertain the name has the proper
# format required by the database: First Last <user@domain>
#

int	procedure	name_check (name)

char	name[ARB]		# name to check

int	n
int	patlen, alloclen
pointer	patbuf
pointer	sp

int	patmake()
int	patmatch()

begin
	alloclen = SZ_LINE

	call smark (sp)
	call salloc (patbuf, alloclen+1, TY_CHAR)

#	patlen = patmake ("[A-Z,a-z]?*<[A-Z,a-z]?*@[A-Z,a-z]?*.[A-Z,a-z]?*>", Memc[patbuf], alloclen)
	patlen = patmake ("[A-Z,a-z]?*<[^ ]?*@[^ ]?*.[A-Z,a-z]?*>", Memc[patbuf], alloclen)
	if (patlen > alloclen)
		call fatal (0, "NAME_CHECK: patlen too large")

	n = patmatch (name, Memc[patbuf]) 

	call sfree (sp)

	if (n > 0)
		return (OK)
	else
		return (ERR)
end


