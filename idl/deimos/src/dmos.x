# TBD: determine amp, chip from header  -- change to use DETSEC
# There is nothing in the header to describe CCDs or the Mosaic layout.
# Thus, some values must be coded in (CCDSZ_X, CCDSZ_Y)

# TBD: limit checking in Y -- needed!!  DONE
# note that flag "any_data" in chip_sect unused

# NOTE: There is an issue of IRAF vs other conventions in terms of pixel
# location.  I have adopted the IRAF conv., which is that the center of the
# first pixel is 1.0, rather than 0.5.

# TBD: need to worry about number of im buffers? I don't thinks so, provided
# we can have as many _images_ open as we want

#
# MOS_INIT: initialize data structure mosaic.
#	Inputs: image name, ccdprop file
#	Outputs: header pointer, file of layout/properties, # CCDs
# 

define	CCDSZ_X		2048.
define	CCDSZ_Y		4096.
define	DEF_PROP_FILE	"deimos$lib/prop/redprop.txt"	# default properties

include	<imhdr.h>
include	<error.h>
include	"dmos.h"

procedure	mos_init (image, ccdprop, im0, mmap, nccd)

char	image[ARB]		# image name
char	ccdprop[ARB]		# file name of CCD properties
pointer	im0			# image pointer to primary HDU
pointer	mmap			# pointer to mosaic vectors
int	nccd

char	imagext[SZ_FNAME]	# name of image with extention
char	tchar
int	nhdu
int	chip, amp
int	i
int	i1, i2, j1, j2		# datasec ranges

real	rpix, rval, dval
int	mx, my

pointer	im
pointer	fd

int	fscan(), nscan()
int	imgeti()
real	imgetr()
pointer	open()
pointer	immap()

bool	imaccf()
int	get_datasec()

begin

# Map the initial header:
# get: nextn, numamps, precol, postcol [pre/post row?]
	call sprintf (imagext, SZ_FNAME, "%-s[0]")
		call pargstr (image)
	im0 = immap (imagext, READ_ONLY, 0)

	if (imaccf (im0, "NVIDINP")) {
		nhdu = imgeti (im0, "NVIDINP")
	} else if (imaccf (im0, "NUMAMPS")) {
		nhdu = imgeti (im0, "NUMAMPS")
	} else {
		call fatal (0, "Number of HDUs cannot be determined ...")
	}
	nccd = 8				# TMP -- FIX should be var
call eprintf ("Fix NHDU/NCCD issue in MOSINIT\n")

# Allocate the vectors
	call malloc (mmap, MMAPLEN, TY_STRUCT)
	call calloc (PTHDU1(mmap), nccd, TY_INT)	# zeroed!
	call calloc (PTHDU2(mmap), nccd, TY_INT)	# zeroed!
	call malloc (PTCPX1(mmap), nccd, TY_INT)
	call malloc (PTC0X1(mmap), nccd, TY_INT)
	call malloc (PTC1X1(mmap), nccd, TY_INT)
	call malloc (PTCPY1(mmap), nccd, TY_INT)
	call malloc (PTC0Y1(mmap), nccd, TY_INT)
	call malloc (PTC1Y1(mmap), nccd, TY_INT)
	call malloc (PTCPX2(mmap), nccd, TY_INT)
	call malloc (PTC0X2(mmap), nccd, TY_INT)
	call malloc (PTC1X2(mmap), nccd, TY_INT)
	call malloc (PTCPY2(mmap), nccd, TY_INT)
	call malloc (PTC0Y2(mmap), nccd, TY_INT)
	call malloc (PTC1Y2(mmap), nccd, TY_INT)

	call malloc (PTGF1(mmap), nccd, TY_REAL)
	call malloc (PTBL1(mmap), nccd, TY_REAL)
	call malloc (PTGF2(mmap), nccd, TY_REAL)
	call malloc (PTBL2(mmap), nccd, TY_REAL)

# These values from HDU0
	PRECOL(mmap) = imgeti (im0, "PRECOL")
	POSTCOL(mmap) = imgeti (im0, "POSTPIX")

# Initialize the flag for freeing the image buffer
	IMB_USED(mmap) = NO

# Loop through all the extensions and fill in table:
call eprintf ("[")
	do i = 1, nhdu {
# ... construct name
		call sprintf (imagext, SZ_FNAME, "%-s[%d]")
			call pargstr (image)
			call pargi (i)
call eprintf (" %d")
call pargi (i)

		im = immap (imagext, READ_ONLY, 0)
# ... identify region.  This should probably be done with DETSEC -- TMP?

		chip = imgeti (im, "CCDLOC")
		amp = imgeti (im, "AMPLOC") - 2*(chip-1)
		mx = chip
		my = 1
		if (chip > 4) {
			my = 2
			mx = mx - 4
		}

		if (get_datasec (im, i1, i2, j1, j2) != OK)
			call fatal (0, "MOS_INIT: GET_DATASEC failed")

		if (amp == 1) {
# WORK HERE
		    HDU1(mmap,chip) = im
		    rpix = imgetr (im, "CRPIX1P")
		    rval = imgetr (im, "CRVAL1P")	# - (mx-1)*CCDSZ_X
		    dval = imgeti (im, "CD1_1P")
		    CPX1(mmap,chip) = i1
		    C0X1(mmap,chip) = rval + (i1 - rpix) * dval + 0.5
		    C1X1(mmap,chip) = rval + (i2 - rpix) * dval + 0.5

		    rpix = imgetr (im, "CRPIX2P")
		    rval = imgetr (im, "CRVAL2P")	# - (my-1)*CCDSZ_Y
		    dval = imgeti (im, "CD2_2P")
		    CPY1(mmap,chip) = j1
		    C0Y1(mmap,chip) = rval + (j1 - rpix) * dval + 0.5
		    C1Y1(mmap,chip) = rval + (j2 - rpix) * dval + 0.5
		} else {
		    HDU2(mmap,chip) = im
		    rpix = imgetr (im, "CRPIX1P")
		    rval = imgetr (im, "CRVAL1P")	# - (mx-1)*CCDSZ_X
		    dval = imgeti (im, "CD1_1P")
		    CPX2(mmap,chip) = i1
		    C0X2(mmap,chip) = rval + (i1 - rpix) * dval + 0.5
		    C1X2(mmap,chip) = rval + (i2 - rpix) * dval + 0.5

		    rpix = imgetr (im, "CRPIX2P")
		    rval = imgetr (im, "CRVAL2P")	# - (my-1)*CCDSZ_Y
		    dval = imgeti (im, "CD2_2P")
		    CPY2(mmap,chip) = j1
		    C0Y2(mmap,chip) = rval + (j1 - rpix) * dval + 0.5
		    C1Y2(mmap,chip) = rval + (j2 - rpix) * dval + 0.5
		}
	}
call eprintf (" ]\n")

# Fill props with defaults
	do i = 1, nccd {
		BL1(mmap,i) = 0.
		BL2(mmap,i) = 0.
		GF1(mmap,i) = 1.
		GF2(mmap,i) = 1.
	}

# Open text file
	fd = open (ccdprop, READ_ONLY, TEXT_FILE)

# Read in Bias, Gain data
	while (fscan (fd) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#' || nscan() == 0) {
			next
		}
		call reset_scan ()

		call gargi (chip)
		call gargr (GF1(mmap,chip))
		call gargr (BL1(mmap,chip))
		call gargr (GF2(mmap,chip))
		call gargr (BL2(mmap,chip))

		if (nscan() < 5) {
			call fatal (0, "MOS_INIT: Bad format in CCDprop")
		}
	}
	call close (fd)

do i = 1, nccd {
if (HDU1(mmap,i) > 0) {
call eprintf ("DEBUG: HDU1=%-8d  X: %4d %4d %4d;  Y: %4d %4d %4d\n")
call pargi (HDU1(mmap,i))
call pargi (CPX1(mmap,i))
call pargi (C0X1(mmap,i))
call pargi (C1X1(mmap,i))
call pargi (CPY1(mmap,i))
call pargi (C0Y1(mmap,i))
call pargi (C1Y1(mmap,i))
call pargi (C1Y2(mmap,2))
}
if (HDU2(mmap,i) > 0) {
call eprintf ("DEBUG: HDU2=%-8d  X: %4d %4d %4d;  Y: %4d %4d %4d\n")
call pargi (HDU2(mmap,i))
call pargi (CPX2(mmap,i))
call pargi (C0X2(mmap,i))
call pargi (C1X2(mmap,i))
call pargi (CPY2(mmap,i))
call pargi (C0Y2(mmap,i))
call pargi (C1Y2(mmap,i))
}
}

end


#
# CHIP_SECT: get image buffer from a given chip, ICS orientation
#

pointer	procedure chip_sect (i1, i2, j1, j2, mmap, chip, biasopt)

int	i1, i2, j1, j2			# section of chip, ICS orientation
pointer	mmap				# mosaic map structure
int	chip				# CCD number
int	biasopt				# bias option [UNUSED]

int	nx, ny

int	i	# TMP
int	k1, k2, l1, l2
int	delx, dely
int	len			# length of subsection buffer
int	xoff, yoff		# offsets in output array to allow for off-image
int	add_data		# Is there data on this chip/in this HDU?
int	any_data		# is there any data in the output buffer?
pointer	bufs
pointer	bufcs

pointer	imgs2r()
begin
# If image buffer has been used before, free the existing memory
	if (IMB_USED(mmap) == YES) {
		call mfree (IMB_BUFF(mmap), TY_REAL)
		IMB_USED(mmap) = NO
	}

## FIX NEEDS WORK: put in checks, trucations on image section

	nx = i2 - i1 + 1
	ny = j2 - j1 + 1

# Allocate buffer space; fill with zeroes as a convenience
	call calloc (bufcs, nx*ny, TY_REAL)
	any_data = NO

# image buffer used, set the flag
	IMB_USED(mmap) = YES
	IMB_BUFF(mmap) = bufcs

# Figure out limits within each amplifier region
	if (HDU1(mmap,chip) > 0) {
		delx = 1
		dely = 1
		if (C1X1(mmap,chip) < C0X1(mmap,chip))
			delx = -1
		if (C1Y1(mmap,chip) < C0Y1(mmap,chip))
			dely = -1

		k1 = (i1 - C0X1(mmap,chip)) * delx
		k2 = (i2 - C0X1(mmap,chip)) * delx
		l1 = (j1 - C0Y1(mmap,chip)) * dely
		l2 = (j2 - C0Y1(mmap,chip)) * dely

# LIMIT CHECKING HERE (Need Y limits also)
		if (delx > 0) {
			k2 = min (C1X1(mmap,chip)-C0X1(mmap,chip), k2)
			xoff = -min (k1, 0)
			k1 = k1 + xoff
			if (k1 <= k2)
				add_data = YES
			else
				add_data = NO
		} else {
			xoff = -min (C0X1(mmap,chip)-C1X1(mmap,chip)-k1, 0)
			k1 = k1 - xoff
			k2 = max (k2, 0)
			if (k2 <= k1)
				add_data = YES
			else
				add_data = NO
		}
		if (dely > 0) {
			l2 = min (C1Y1(mmap,chip)-C0Y1(mmap,chip), l2)
			yoff = -min (l1, 0)
			l1 = l1 + yoff
			if (l1 <= l2 && add_data == YES)
				add_data = YES
			else
				add_data = NO
		} else {
			yoff = -min (C0Y1(mmap,chip)-C1Y1(mmap,chip)-l1, 0)
			l1 = l1 - yoff
			l2 = max (l2, 0)
			if (l2 <= l1 && add_data == YES)
				add_data = YES
			else
				add_data = NO
		}
	} else {
		add_data = NO
	}

	if (add_data == YES) {
# adjust for datasec (skip over non-imaging pixels)
		k1 = k1 + CPX1(mmap,chip)
		k2 = k2 + CPX1(mmap,chip)
		l1 = l1 + CPY1(mmap,chip)
		l2 = l2 + CPY1(mmap,chip)

		bufs = imgs2r (HDU1(mmap,chip), k1, k2, l1, l2)
		len = (abs (k2-k1) + 1) * (abs (l2-l1) + 1)

		if (BIASOPT == YES) {
		    call aaddkr (Memr[bufs], -BL1(mmap,chip), Memr[bufs], len)
		    call amulkr (Memr[bufs], 1./GF1(mmap,chip), Memr[bufs], len)
		}

# Transfer into array
		len = abs (k2-k1) + 1
		do i = 0, ny-1 {
		    call amovr (Memr[bufs+i*len], Memr[bufcs+i*nx+xoff], len)
		}
		any_data = YES
	}

## Check the second HDU and process if appropriate:
#
	if (HDU2(mmap,chip) > 0) {
		delx = 1
		dely = 1
		if (C1X2(mmap,chip) < C0X2(mmap,chip))
			delx = -1
		if (C1Y2(mmap,chip) < C0Y2(mmap,chip))
			dely = -1

		k1 = (i1 - C0X2(mmap,chip)) * delx
		k2 = (i2 - C0X2(mmap,chip)) * delx
		l1 = (j1 - C0Y2(mmap,chip)) * dely
		l2 = (j2 - C0Y2(mmap,chip)) * dely

# LIMIT CHECKING HERE (Need Y limits also)
		if (delx > 0) {
			k2 = min (C1X2(mmap,chip)-C0X2(mmap,chip), k2)
			xoff = -min (k1, 0)
			k1 = k1 + xoff
			if (k1 <= k2)
				add_data = YES
			else
				add_data = NO
		} else {
			xoff = -min (C0X2(mmap,chip)-C1X2(mmap,chip)-k1, 0)
			k1 = k1 - xoff
			k2 = max (k2, 0)
			if (k2 <= k1)
				add_data = YES
			else
				add_data = NO
		}
		if (dely > 0) {
			l2 = min (C1Y2(mmap,chip)-C0Y2(mmap,chip), l2)
			yoff = -min (l1, 0)
			l1 = l1 + yoff
			if (l1 <= l2 && add_data == YES)
				add_data = YES
			else
				add_data = NO
		} else {
			yoff = -min (C0Y2(mmap,chip)-C1Y2(mmap,chip)-l1, 0)
			l1 = l1 - yoff
			l2 = max (l2, 0)
			if (l2 <= l1 && add_data == YES)
				add_data = YES
			else
				add_data = NO
		}
	} else {
		add_data = NO
	}

	if (add_data == YES) {

# adjust for datasec (skip over non-imaging pixels)
		k1 = k1 + CPX2(mmap,chip)
		k2 = k2 + CPX2(mmap,chip)
		l1 = l1 + CPY2(mmap,chip)
		l2 = l2 + CPY2(mmap,chip)

		bufs = imgs2r (HDU2(mmap,chip), k1, k2, l1, l2)
		len = (abs (k2-k1) + 1) * (abs (l2-l1) + 1)

		if (BIASOPT == YES) {
		    call aaddkr (Memr[bufs], -BL2(mmap,chip), Memr[bufs], len)
		    call amulkr (Memr[bufs], 1./GF2(mmap,chip), Memr[bufs], len)
		}

# Transfer into array
		len = abs (k2-k1) + 1
		do i = 0, ny-1 {
		    call amovr (Memr[bufs+i*len], Memr[bufcs+i*nx+xoff], len)
		}
		any_data = YES
	}
		
	return (bufcs)

end
	

#
# GET_DATASEC: decode datasec.  Currently assumes NAXIS=2
#

define	SZ_FITS_STR	72

int	procedure get_datasec (im, i1, i2, j1, j2)

pointer	im		# image descriptor
int	i1, i2		# datasec limits in x
int	j1, j2		# datasec limits in y

char	kwval[SZ_FITS_STR]
char	wkstr[SZ_FITS_STR]
char	tchar
int	i, n, ia, ib, io
int	naxis

int	stridx(), sscan(), nscan()
begin
	naxis = IM_NDIM(im)
	if (naxis != 2)
		call fatal (0, "GET_DATASEC: naxis != 2!")

	call imgstr (im, "DATASEC", kwval, SZ_FITS_STR)

	tchar = '['
	ia = stridx (tchar, kwval)
	tchar = ']'
	ib = stridx (tchar, kwval)

	if (ia < 1 || ib < ia)
		call fatal (0, "GET_DATASEC: no bracket pair")

	tchar = ' '
	call amovkc (tchar, wkstr, SZ_FITS_STR)
	i = ia + 1
	io = 1
	tchar = ':'
	n = stridx (tchar, kwval[i]) - 1
	if (n < 1)
		call fatal (0, "GET_DATASEC: poor datasec format")
	call amovc (kwval[i], wkstr[io], n)

	i = i + n + 1
	io = io + n + 1
	tchar = ','
	n = stridx (tchar, kwval[i]) - 1
	if (n < 1)
		call fatal (0, "GET_DATASEC: poor datasec format")
	call amovc (kwval[i], wkstr[io], n)

	i = i + n + 1
	io = io + n + 1
	tchar = ':'
	n = stridx (tchar, kwval[i]) - 1
	if (n < 1)
		call fatal (0, "GET_DATASEC: poor datasec format")
	call amovc (kwval[i], wkstr[io], n)

	i = i + n + 1
	io = io + n + 1
	tchar = ']'
	n = stridx (tchar, kwval[i]) - 1
	if (n < 1)
		call fatal (0, "GET_DATASEC: poor datasec format")
	call amovc (kwval[i], wkstr[io], n)

	io = io + n
	wkstr[io] = EOS

	i = sscan (wkstr)
		call gargi (i1)
		call gargi (i2)
		call gargi (j1)
		call gargi (j2)
	if (nscan() < 4)
		call fatal (0, "GET_DATASEC: decode failed; bad format?")
	else
		return (OK)
end


procedure	t_tmos ()

char	image[SZ_FNAME]		# image name
# char	ccdprop[SZ_FNAME]	# file name of CCD properties
pointer	im0			# image pointer to primary HDU
pointer	mmap			# pointer to mosaic vectors
int	nccd

int	nx, ny, i1, i2, j1, j2
int	chip
pointer	im2

pointer	immap(), imps2r()
pointer	chip_sect()
int	det_chip()

begin
	call clgstr ("image", image, SZ_FNAME)
#	call clgstr ("ccdprop", ccdprop, SZ_FNAME)
#	call mos_init (image, ccdprop, im0, mmap, nccd)

 	call mos_init (image, DEF_PROP_FILE, im0, mmap, nccd)

	nx = 70
	ny = 70
	i1 = 1000

	i2 = i1 + nx - 1
	j1 = 1191
	j2 = j1 + ny - 1

	im2 = immap ("testim", NEW_IMAGE, 0)
	IM_NDIM(im2) = 2
	IM_LEN(im2,1) = nx
	IM_LEN(im2,2) = ny
	IM_PIXTYPE(im2) = TY_REAL

call eprintf ("DEBUG: calling Chip_Sect, %d:%d,%d:%d\n")
call pargi (i1)
call pargi (i2)
call pargi (j1)
call pargi (j2)
	chip = det_chip (real (i1), real(j1))
	call amovr (Memr[chip_sect (i1, i2, j1, j2, mmap, chip, YES)],
			Memr[imps2r (im2, 1, nx, 1, ny)], nx*ny)
	call imunmap (im2)
end


procedure	mos_free (im0, mmap)

pointer	im0			# image pointer to primary HDU
pointer	mmap			# pointer to mosaic vectors

int	i

begin
	do i = 1, 8 {			# HARDCODE HERE ...
		if (HDU1(mmap,i) > 0)
			call imunmap (HDU1(mmap,i))
		if (HDU2(mmap,i) > 0)
			call imunmap (HDU2(mmap,i))
	}

	call imunmap (im0)
call eprintf ("image unmapped\n")
end

#
# DET_CHIP: like ident_ccd, except uses "detsec" coordinates
#

int	procedure det_chip (xim, yim)

real	xim, yim			# x,y in DETSEC coords

int	mx, my

begin
	
	if (yim > 4096.)
		my = 2
	else
		my = 1

# find appropriate x-division
	if (xim < 2048.)
		mx = 1
	else if (xim < 4096.)
		mx = 2
	else if (xim < 6144.)
		mx = 3
	else
		mx = 4

	return ((my-1)*4+mx)
end
