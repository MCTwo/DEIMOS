include	<time.h>
include	"ftab.h"

# FITSTEST:  look for ways to write fits tables
# Path: create struct to handle format, etc?  compose a format string
# do an sprintf to write spp string
# pass string to write routine, which uses chrpak and write OR simply
#      use fprintf.  Keep count of chars; will need to pad on exit.
# 
## NB: No write permission (String_File) results from the number of characters
#  being too short for the requested write in sprintf()

procedure t_fitstest ()

pointer	fd
char	output[SZ_FNAME]			# output fits file

char	fmtstr[SZ_LINE]				# format for writing table row
int	fmtstrlen
int	n
pointer	kdat

real	rnul
double	dnul
char	cnul[6]

int	i

pointer	open()
begin

# Read in parameters:

	call clgstr ("output", output, SZ_FNAME)
	fd = open (output, NEW_FILE, TEXT_FILE)

	rnul = -9.e9
	dnul = -9.e9
	call strcpy ("INDEF", cnul, 6)

# Write the primary HDU
	call ftab_init (kdat, 1, fmtstr)
	call pkwb (fd, kdat, "SIMPLE", true, "file does conform to FITS standard")
	call pkwi (fd, kdat, "BITPIX", 8, "number of bits per data pixel")
	call pkwi (fd, kdat, "NAXIS", 0, "number of data axes")
	call pkwb (fd, kdat, "EXTEND", true, "FITS dataset may contain extensions")
	call pkwspec (fd, kdat, "END")
	do n = 6, 36
		call pkwspec (fd, kdat, "")

	if (CHARCNT(kdat) != 0)
		call eprintf ("CHAR_CNT off after primary HDU\n")

	call ftab_free (fd, kdat)



	fmtstrlen = SZ_LINE
	n = 33

	call ftab_init (kdat, 4, fmtstr)

	call ftcol_defi (kdat, "index", "I6",      "", -9999,  6, fmtstr, fmtstrlen)
	call ftcol_defr (kdat, "xslit", "E12.6", "mm", rnul, 12, fmtstr, fmtstrlen)
	call ftcol_defd (kdat, "yslit", "E12.6", "mm", dnul, 12, fmtstr, fmtstrlen)
	call ftcol_defc (kdat, "id",    "A10",     "", cnul, 10, fmtstr, fmtstrlen)


	call ftab_whead (fd, kdat, n)


	do i = 1, n {
		call sprintf (BUFROW(kdat), TWUSED(kdat), fmtstr)
			call pargi (i)
			call pargr (rnul)
			call pargd (dnul)
#			call pargr (3.1415926)
#			call pargr (99.9999e9)
			call pargstr ("ID_NAME")
		call ftab_wrow (fd, kdat)
	}

	call ftab_free (fd, kdat)


# Do it again
	fmtstrlen = SZ_LINE
	n = 4

	call ftab_init (kdat, 3, fmtstr)

	call ftcol_defi (kdat, "index", "I6",      "", -999,  6, fmtstr, fmtstrlen)
	call ftcol_defr (kdat, "xslit", "E12.6", "mm", rnul, 12, fmtstr, fmtstrlen)
	call ftcol_defc (kdat, "id",    "A20",     "", "INDEF", 20, fmtstr, fmtstrlen)

	call ftab_whead (fd, kdat, n)

	do i = 1, n {
		call sprintf (BUFROW(kdat), TWUSED(kdat), fmtstr)
			call pargi (i)
			call pargr (3.1415926)
			call pargstr ("ANOTHER_NAME")
		call ftab_wrow (fd, kdat)
	}


	call ftab_free (fd, kdat)


	call close (fd)

end

#############################################################################
# 
#	ftab_init (kdat, ncol, fmtstr)
#
#	ftcol_defc (kdat, type, form, unit, nulstr, wid, fmtstr, fmtstrlen)
#	ftcol_defi (kdat, type, form, unit, nulstr, wid, fmtstr, fmtstrlen)
#	ftcol_defr (kdat, type, form, unit, nulstr, wid, fmtstr, fmtstrlen)
#	ftcol_defd (kdat, type, form, unit, nulstr, wid, fmtstr, fmtstrlen)
#
#	ftab_whead (fd, kdat, nline)
#
#	ftab_wrow (fd, kdat)
#
#	ftab_free (fd, kdat)
#
#	pkwi   (fd, kdat, name, val, comment)
#	pkwr   (fd, kdat, name, val, comment)
#	pkwstr (fd, kdat, name, val, comment)
#	pkwspec (fd, kdat, val)
#

# FTAB_INIT: initialize fits acsii table
# 
procedure	ftab_init (kdat, ncol, fmtstr)

pointer	kdat
int	ncol
char	fmtstr[ARB]

begin
	call malloc (kdat, KDATLEN, TY_STRUCT)
	call malloc (TBCOLPT(kdat), ncol, TY_INT)
	call malloc (TFORMPT(kdat), ncol*TFMTLEN, TY_CHAR)
	call malloc (TTYPEPT(kdat), ncol*TTYPLEN, TY_CHAR)
	call malloc (TUNITPT(kdat), ncol*TUNILEN, TY_CHAR)
	call malloc (TNULLPT(kdat), ncol*TNULLEN, TY_CHAR)

	TROWPT (kdat) = 0		# set to note that NOT assigned yet
	TMAXCOL(kdat) = ncol
	TFIELD(kdat) = 0
	TWUSED(kdat) = 0
	CHARCNT(kdat) = 0
	THDROK(kdat) = NO
	call strcpy ("", fmtstr, 1)
end


# FTCOL_DEF: define ascii table column, return update format string

procedure ftcol_def (kdat, type, form, unit, null, wid, fmtstr, fmtstrlen)
# procedure ftcol_def (kdat, type, form, unit, wid, fmtstr, fmtstrlen)

pointer	kdat				# pointer to struct
char	type[ARB]			# field name
char	form[ARB]			# format string for column
char	unit[ARB]			# field units
int	null				# null value (addr just gets passed)
int	wid				# width of field (TMP?)

char	fmtstr[ARB]
int	fmtstrlen

char	sppfmt[16]
char	nullstr[68]
int	ndx

int	strlen()
begin
	ndx = TFIELD(kdat)
	if (ndx + 1 > TMAXCOL(kdat))
		call fatal (0, "ftcol_def: more columns defined than requested")

	if (strlen (form) >= TFMTLEN)
		call fatal (0, "ftcol_def: format string too long")
	if (strlen (type) >= TTYPLEN)
		call fatal (0, "ftcol_def: type string too long")
	if (strlen (unit) >= TUNILEN)
		call fatal (0, "ftcol_def: unit string too long")


	TFIELD(kdat) = TFIELD(kdat) + 1
	TBCOL(kdat,ndx) = TWUSED(kdat) + 1
	TWUSED(kdat) = TWUSED(kdat) + wid

# translate format string
	if (form[1] == 'A' || form[1] == 'a') {
		call sprintf (sppfmt, 16, "%%-%ss")
			call pargstr (form[2])
		call sprintf (nullstr, 68, "%-s")
			call pargstr (null)
	} else if (form[1] == 'I' || form[1] == 'i') {
		call sprintf (sppfmt, 16, "%%%sd")
			call pargstr (form[2])
		call sprintf (nullstr, 68, sppfmt)
			call pargi (null)
	} else if (form[1] == 'F' || form[1] == 'f') {
		call sprintf (sppfmt, 16, "%%%sf")
			call pargstr (form[2])
		call sprintf (nullstr, 68, sppfmt)
			call pargr (null)
	} else if (form[1] == 'E' || form[1] == 'e') {
		call sprintf (sppfmt, 16, "%%%sg")
			call pargstr (form[2])
		call sprintf (nullstr, 68, sppfmt)
			call pargr (null)
	} else {
		call fatal (0, "ftcol_def:  no translation for format type")
	}


	if (strlen (fmtstr) + strlen (sppfmt) > fmtstrlen) {
		call fatal (0, "requested format string too long!")
	}

	call strcat (sppfmt, fmtstr, fmtstrlen)

	call strcpy (form, TFORM(kdat,ndx), TFMTLEN-1)
	call strcpy (type, TTYPE(kdat,ndx), TTYPLEN-1)
	call strcpy (unit, TUNIT(kdat,ndx), TUNILEN-1)
	call strcpy (nullstr, TNULL(kdat,ndx), TNULLEN-1)

end

# FTAB_WHEAD:  Write header including column defns
# Note: header is left open; first data write closes it

procedure	ftab_whead (fd, kdat, nline)

pointer	fd				# opened file
pointer	kdat				# tab KW struct
int	nline				# number of lines to write

char	kwd[8]				# string variable for enumerated keyword
int	n

bool	strne()
begin
	call pkwstr (fd, kdat, "XTENSION", "TABLE   ", "ASCII table extension")
	call pkwi   (fd, kdat, "BITPIX", 8, "number of bits per data pixel")
	call pkwi   (fd, kdat, "NAXIS", 2, "always 2 for a FITS table")
	call pkwi   (fd, kdat, "NAXIS1", TWUSED(kdat), "cols -- chars needed to describe fields")
	call pkwi   (fd, kdat, "NAXIS2", nline, "rows")
	call pkwi   (fd, kdat, "PCOUNT", 0, "required FITS extension keyword")
	call pkwi   (fd, kdat, "GCOUNT", 1, "required FITS extension keyword")
	call pkwi   (fd, kdat, "TFIELDS", TFIELD(kdat), "Number of columns in the table")


	if (TFIELD(kdat) > 0) {
	    do n = 0, TFIELD(kdat)-1 {
		call sprintf (kwd, 8, "TBCOL%-3d")
			call pargi (n+1)
		call pkwi (fd, kdat, kwd, TBCOL(kdat,n), "")

		call sprintf (kwd, 8, "TFORM%-3d")
			call pargi (n+1)
		call pkwstr (fd, kdat, kwd, TFORM(kdat,n), "")

		if (strne (TTYPE(kdat,n), "")) {
		    call sprintf (kwd, 8, "TTYPE%-3d")
			call pargi (n+1)
		    call pkwstr (fd, kdat, kwd, TTYPE(kdat,n), "")
		}

		if (strne (TUNIT(kdat,n), "")) {
		    call sprintf (kwd, 8, "TUNIT%-3d")
			call pargi (n+1)
		    if (strne (TUNIT(kdat,n), EOS))
			call pkwstr (fd, kdat, kwd, TUNIT(kdat,n), "")
		}

		if (strne (TNULL(kdat,n), "")) {
		    call sprintf (kwd, 8, "TNULL%-3d")
			call pargi (n+1)
		    if (strne (TNULL(kdat,n), EOS))
			call pkwstr (fd, kdat, kwd, TNULL(kdat,n), "")
		}
	    }
	}

# Finally, create the data (row) buffer
		call malloc (TROWPT(kdat), TWUSED(kdat)+1, TY_CHAR)

end	



# FTAB_WROW: write row, update count

procedure	ftab_wrow (fd, kdat)
pointer	fd
pointer	kdat

int	cnt, n

begin
	if (THDROK(kdat) == NO) {
		call pkwspec (fd, kdat, "END")

		cnt = (2880 - CHARCNT(kdat)) / 80
		if (cnt > 0) {
			do n = 1, cnt
				call pkwspec (fd, kdat, "")
		}

# and reset counter ?
		if (CHARCNT(kdat) != 0) {
			call eprintf ("Error on CHAR_CNT: %d\n")
				call pargi (CHARCNT(kdat))
			CHARCNT(kdat) = 0
		}

		THDROK(kdat) = YES
	}
	call fprintf (fd, "%s")
		call pargstr (BUFROW(kdat))
	CHARCNT(kdat) = CHARCNT(kdat) + TWUSED(kdat)	# actually char-cnt
	if (CHARCNT(kdat) > 2880)
		CHARCNT(kdat) = CHARCNT(kdat) - 2880
end



# FTAB_CLOSE: pad file, release memory

procedure	ftab_free (fd, kdat)

pointer	fd
pointer	kdat

int	i

begin
# First, pad and close the fits table  (REVIEW -- v. inefficient)
# Could be better done using salloc, TY_CHAR, then amovkc and and pargstr
# NB In that case, terminate string!!

	if (CHARCNT(kdat) > 0) {
		do i = CHARCNT(kdat)+1, 2880
			call fprintf (fd, " ")
	}
	
# Free memory
	if (TROWPT(kdat) != 0)
		call mfree (TROWPT(kdat), TY_CHAR)
	call mfree (TNULLPT(kdat), TY_CHAR)
	call mfree (TUNITPT(kdat), TY_CHAR)
	call mfree (TTYPEPT(kdat), TY_CHAR)
	call mfree (TFORMPT(kdat), TY_CHAR)
	call mfree (TBCOLPT(kdat), TY_INT)
	call mfree (kdat, TY_STRUCT)
end

procedure pkwb (fd, kdat, name, val, comment)

pointer	fd
pointer	kdat
char	name[ARB]
bool	val
char	comment[ARB]

char	card[80]

begin
	call sprintf (card, 80, "%-8s= %20s / %s")
		call pargstr (name)
		if (val)
			call pargc ("T")
		else
			call pargc ("F")
		call pargstr (comment)

	call fprintf (fd, "%-80s")
		call pargstr (card)

	CHARCNT(kdat) = CHARCNT(kdat) + 80
	if (CHARCNT(kdat) >= 2880)
		CHARCNT(kdat) = CHARCNT(kdat) - 2880
end

procedure pkwi (fd, kdat, name, val, comment)

pointer	fd
pointer	kdat
char	name[ARB]
int	val
char	comment[ARB]

char	card[80]

begin
	call sprintf (card, 80, "%-8s= %20d / %s")
		call pargstr (name)
		call pargi (val)
		call pargstr (comment)

	call fprintf (fd, "%-80s")
		call pargstr (card)

	CHARCNT(kdat) = CHARCNT(kdat) + 80
	if (CHARCNT(kdat) >= 2880)
		CHARCNT(kdat) = CHARCNT(kdat) - 2880
end


procedure pkwr (fd, kdat, name, val, comment)

pointer	fd
pointer	kdat
char	name[ARB]
real	val
char	comment[ARB]

char	card[80]

begin
	call sprintf (card, 80, "%-8s= %20f / %s")
		call pargstr (name)
		call pargr (val)
		call pargstr (comment)

	call fprintf (fd, "%-80s")
		call pargstr (card)

	CHARCNT(kdat) = CHARCNT(kdat) + 80
	if (CHARCNT(kdat) >= 2880)
		CHARCNT(kdat) = CHARCNT(kdat) - 2880
end


procedure pkwd (fd, kdat, name, val, comment)

pointer	fd
pointer	kdat
char	name[ARB]
double	val
char	comment[ARB]

char	card[80]

begin
	call sprintf (card, 80, "%-8s= %20f / %s")
		call pargstr (name)
		call pargd (val)
		call pargstr (comment)

	call fprintf (fd, "%-80s")
		call pargstr (card)

	CHARCNT(kdat) = CHARCNT(kdat) + 80
	if (CHARCNT(kdat) >= 2880)
		CHARCNT(kdat) = CHARCNT(kdat) - 2880
end

	
procedure pkwstr (fd, kdat, name, val, comment)

pointer	fd
pointer	kdat
char	name[ARB]
char	val[ARB]
char	comment[ARB]

char	card[80]
char	tval[70]

int	strlen()
begin
	call sprintf (tval, 69, "'%s")
		call pargstr (val)
	call strcat ("'", tval, 70)

	if (strlen (val) > 68)
		call eprintf ("WARNING: fits string value truncated!\n")

	if (strlen (tval) < 20)		# pad for appearance
		call strcat ("                    ", tval, 20)

	call sprintf (card, 80, "%-8s= %s / %s")
		call pargstr (name)
		call pargstr (tval)
		call pargstr (comment)

	call fprintf (fd, "%-80s")
		call pargstr (card)

	CHARCNT(kdat) = CHARCNT(kdat) + 80
	if (CHARCNT(kdat) >= 2880)
		CHARCNT(kdat) = CHARCNT(kdat) - 2880
end


procedure pkwspec (fd, kdat, val)

pointer	fd
pointer	kdat
char	val[ARB]

char	card[80]

begin
	call sprintf (card, 80, "%-s")
		call pargstr (val)

	call fprintf (fd, "%-80s")
		call pargstr (card)

	CHARCNT(kdat) = CHARCNT(kdat) + 80
	if (CHARCNT(kdat) >= 2880)
		CHARCNT(kdat) = CHARCNT(kdat) - 2880
end

#
# FTAB_WHDU0: write a simple FITS primary HDU for tables
#

procedure	ftab_whdu0 (fd)

pointer	fd				# file descriptor

char	fmtstr[SZ_LINE]
int	i
pointer	kdat

begin
	call ftab_init (kdat, 1, fmtstr)
	call pkwb (fd, kdat, "SIMPLE", true, "file does conform to FITS standard")
	call pkwi (fd, kdat, "BITPIX", 8, "number of bits per data pixel")
	call pkwi (fd, kdat, "NAXIS", 0, "number of data axes")
	call pkwb (fd, kdat, "EXTEND", true, "FITS dataset may contain extensions")
	call pkwspec (fd, kdat, "END")
	do i = 6, 36
		call pkwspec (fd, kdat, "")

	if (CHARCNT(kdat) != 0)
		call eprintf ("CHAR_CNT off after primary HDU\n")

	call mfree (TNULLPT(kdat), TY_CHAR)
	call mfree (TUNITPT(kdat), TY_CHAR)
	call mfree (TTYPEPT(kdat), TY_CHAR)
	call mfree (TFORMPT(kdat), TY_CHAR)
	call mfree (TBCOLPT(kdat), TY_INT)
	call mfree (kdat, TY_STRUCT)
end

#
# FT_DATE: (new) FITS date format
#

procedure	ft_date (date, size)

char	date[size]		# date/time string
int	size			# size

long	ltime			# seconds since 00:00:00 10-Jan-1980
int	tm[LEN_TMSTRUCT]	# broken down time structure

long	clktime()
begin
	if (size < 19)
		call eprintf ("FT_TIME: time string truncated\n")

	ltime = clktime (long (0))
	call brktime (ltime, tm)
	call sprintf (date, size, "%4d-%02d-%02dT%02d:%02d:%02d")
		call pargi (TM_YEAR(tm))
		call pargi (TM_MONTH(tm))
		call pargi (TM_MDAY(tm))
		call pargi (TM_HOUR(tm))
		call pargi (TM_MIN(tm))
		call pargi (TM_SEC(tm))
end


# FTCOL_DEFC: define text-type ascii table column, return update format string

procedure ftcol_defc (kdat, type, form, unit, null, wid, fmtstr, fmtstrlen)

pointer	kdat				# pointer to struct
char	type[ARB]			# field name
char	form[ARB]			# format string for column
char	unit[ARB]			# field units
char	null[ARB]			# null value (addr just gets passed)
int	wid				# width of field (TMP?)

char	fmtstr[ARB]
int	fmtstrlen

int	typefmt				# type of format (fixed, float, char)

char	sppfmt[16]
char	nullstr[68]
int	ndx

int	strlen()
begin
	call ftcol_setup (kdat, type, form, unit, null, wid, fmtstr, fmtstrlen,
								typefmt, sppfmt)

	if (typefmt != FT_CHAR)			# SPECIFIC TO TYPE
		call fatal (0, "ftcol_defc:  wrong format type for call")

#	call sprintf (nullstr, 68, sppfmt)
#		call pargstr (null)		# SPECIFIC TO TYPE
	call strcpy (null, nullstr, 68)

	if (strlen (nullstr) >= TNULLEN)
		call fatal (0, "ftcol_defc: null string too long")

	ndx = TFIELD(kdat) - 1
	call strcpy (form, TFORM(kdat,ndx), TFMTLEN-1)
	call strcpy (type, TTYPE(kdat,ndx), TTYPLEN-1)
	call strcpy (unit, TUNIT(kdat,ndx), TUNILEN-1)
	call strcpy (nullstr, TNULL(kdat,ndx), TNULLEN-1)

end

# FTCOL_DEFI: define int-type ascii table column, return update format string

procedure ftcol_defi (kdat, type, form, unit, null, wid, fmtstr, fmtstrlen)

pointer	kdat				# pointer to struct
char	type[ARB]			# field name
char	form[ARB]			# format string for column
char	unit[ARB]			# field units
int	null				# null value (addr just gets passed)
int	wid				# width of field (TMP?)

char	fmtstr[ARB]
int	fmtstrlen

int	typefmt				# type of format (fixed, float, char)

char	sppfmt[16]
char	nullstr[68]
int	ndx

int	strlen()
begin
	call ftcol_setup (kdat, type, form, unit, null, wid, fmtstr, fmtstrlen,
								typefmt, sppfmt)

	if (typefmt != FT_FIXED)			# SPECIFIC TO TYPE
		call fatal (0, "ftcol_defi:  wrong format type for call")

	call sprintf (nullstr, 68, sppfmt)
		call pargi (null)		# SPECIFIC TO TYPE

	if (strlen (nullstr) >= TNULLEN)
		call fatal (0, "ftcol_defi: null string too long")

	ndx = TFIELD(kdat) - 1
	call strcpy (form, TFORM(kdat,ndx), TFMTLEN-1)
	call strcpy (type, TTYPE(kdat,ndx), TTYPLEN-1)
	call strcpy (unit, TUNIT(kdat,ndx), TUNILEN-1)
	call strcpy (nullstr, TNULL(kdat,ndx), TNULLEN-1)

end

# FTCOL_DEFR: define real-type ascii table column, return update format string

procedure ftcol_defr (kdat, type, form, unit, null, wid, fmtstr, fmtstrlen)

pointer	kdat				# pointer to struct
char	type[ARB]			# field name
char	form[ARB]			# format string for column
char	unit[ARB]			# field units
real	null				# null value (addr just gets passed)
int	wid				# width of field (TMP?)

char	fmtstr[ARB]
int	fmtstrlen

int	typefmt				# type of format (fixed, float, char)

char	sppfmt[16]
char	nullstr[68]
int	ndx

int	strlen()
begin
	call ftcol_setup (kdat, type, form, unit, null, wid, fmtstr, fmtstrlen,
								typefmt, sppfmt)

	if (typefmt != FT_FLOAT)		# SPECIFIC TO TYPE
		call fatal (0, "ftcol_defr:  wrong format type for call")

	call sprintf (nullstr, 68, sppfmt)
		call pargr (null)		# SPECIFIC TO TYPE

	if (strlen (nullstr) >= TNULLEN)
		call fatal (0, "ftcol_defr: null string too long")

	ndx = TFIELD(kdat) - 1
	call strcpy (form, TFORM(kdat,ndx), TFMTLEN-1)
	call strcpy (type, TTYPE(kdat,ndx), TTYPLEN-1)
	call strcpy (unit, TUNIT(kdat,ndx), TUNILEN-1)
	call strcpy (nullstr, TNULL(kdat,ndx), TNULLEN-1)

end

# FTCOL_DEFd: define double-type ascii table column, return update format string

procedure ftcol_defd (kdat, type, form, unit, null, wid, fmtstr, fmtstrlen)

pointer	kdat				# pointer to struct
char	type[ARB]			# field name
char	form[ARB]			# format string for column
char	unit[ARB]			# field units
double	null				# null value (addr just gets passed)
int	wid				# width of field (TMP?)

char	fmtstr[ARB]
int	fmtstrlen

int	typefmt				# type of format (fixed, float, char)

char	sppfmt[16]
char	nullstr[68]
int	ndx

int	strlen()
begin
	call ftcol_setup (kdat, type, form, unit, null, wid, fmtstr, fmtstrlen,
								typefmt, sppfmt)

	if (typefmt != FT_FLOAT)		# SPECIFIC TO TYPE
		call fatal (0, "ftcol_defd:  wrong format type for call")

	call sprintf (nullstr, 68, sppfmt)
		call pargd (null)		# SPECIFIC TO TYPE

	if (strlen (nullstr) >= TNULLEN)
		call fatal (0, "ftcol_defd: null string too long")

	ndx = TFIELD(kdat) - 1
	call strcpy (form, TFORM(kdat,ndx), TFMTLEN-1)
	call strcpy (type, TTYPE(kdat,ndx), TTYPLEN-1)
	call strcpy (unit, TUNIT(kdat,ndx), TUNILEN-1)
	call strcpy (nullstr, TNULL(kdat,ndx), TNULLEN-1)

end

## FTCOL_SETUP: Generic, type-independent setup code

procedure ftcol_setup (kdat, type, form, unit, null, wid, fmtstr, fmtstrlen, typefmt, sppfmt)

pointer	kdat				# pointer to struct
char	type[ARB]			# field name
char	form[ARB]			# format string for column
char	unit[ARB]			# field units
int	null				# null value (addr just gets passed)
int	wid				# width of field (TMP?)

char	fmtstr[ARB]
int	fmtstrlen
int	typefmt				# type of format (fixed, float, char)
char	sppfmt[16]

int	ndx

int	strlen()
begin
	ndx = TFIELD(kdat)
	if (ndx + 1 > TMAXCOL(kdat))
		call fatal (0, "ftcol_def: more columns defined than requested")

	if (strlen (form) >= TFMTLEN)
		call fatal (0, "ftcol_def: format string too long")
	if (strlen (type) >= TTYPLEN)
		call fatal (0, "ftcol_def: type string too long")
	if (strlen (unit) >= TUNILEN)
		call fatal (0, "ftcol_def: unit string too long")

	TFIELD(kdat) = TFIELD(kdat) + 1
	TBCOL(kdat,ndx) = TWUSED(kdat) + 1
	TWUSED(kdat) = TWUSED(kdat) + wid

# translate format string
	if (form[1] == 'A' || form[1] == 'a') {
		call sprintf (sppfmt, 16, "%%-%ss")
			call pargstr (form[2])
		typefmt = FT_CHAR
	} else if (form[1] == 'I' || form[1] == 'i') {
		call sprintf (sppfmt, 16, "%%%sd")
			call pargstr (form[2])
		typefmt = FT_FIXED
	} else if (form[1] == 'F' || form[1] == 'f') {
		call sprintf (sppfmt, 16, "%%%sf")
			call pargstr (form[2])
		typefmt = FT_FLOAT
	} else if (form[1] == 'E' || form[1] == 'e') {
		call sprintf (sppfmt, 16, "%%%sg")
			call pargstr (form[2])
		typefmt = FT_FLOAT
	} else {
		call fatal (0, "ftcol_def:  no translation for format type")
	}

	if (strlen (fmtstr) + strlen (sppfmt) > fmtstrlen) {
		call fatal (0, "requested format string too long!")
	}

	call strcat (sppfmt, fmtstr, fmtstrlen)
end
