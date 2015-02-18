#
# These are the graphical routines for DSIMULATOR
#

include	<math.h>
include	<gset.h>
include	<gim.h>
include	"dsimulator.h"
include	"deimos.h"

define	FAST	10.
define	MEDIUM	1.
define	SLOW	0.1
define	LBL1	"-X"
define	LBL2	"+Y"
define	LBL3	"Seen From Front"
define	DEF_FLD_WIDTH	1100.
define	DEF_FLD_CENX	0.
define	DEF_FLD_CENY	180.


procedure	gselect_targ (tdat, ntarg, sdat, nslit, indat)

pointer	tdat
int	ntarg
pointer	sdat
int	nslit
pointer	indat

pointer	mdat			# pointer to data structure for mask outlines

char	objname[SZ_ID]
int	i, j
int	nlist			# number of the active list
int	psum

real	anginc, xyinc, speed
real	xoff, yoff, aoff, doff
real	cosp, sinp		# could be double
real	cosdec			# could be double
real	zfact

char	command[32]			# not sure if 32 is good
int	wcs, key
real	wx, wy
pointer	gp

pointer	gopen()
int	clgcur(), get_sw_nearest0()
int	get_valr()
int	sscan()
begin
	nlist = PRIMARY		## XXX TMP DEBUG

# Read in the mask outlines (also defines window params)
	call get_outlines (mdat)

	FPWIN_X(mdat) = DEF_FLD_CENX
	FPWIN_Y(mdat) = DEF_FLD_CENY
	FPWIN_W(mdat) = DEF_FLD_WIDTH


# Open the graphics stream
	gp = gopen ("stdgraph", NEW_FILE, STDGRAPH)

	call fp_layout (gp, mdat, indat, tdat, ntarg, sdat, nslit, nlist, LBL1, LBL2, LBL3)
	if (sdat != 0 && nslit > 0) {
		do i = 0, nslit-1
			call mark_slit (gp, sdat, i)
	}

# start the interactive stuff

	speed = FAST
	anginc = DEGTORAD(1.) * speed		# TMP_HARDCODE (was PA_INCR)
	xyinc = 6. * speed		# TMP_HARDCODE (was TR_INCR)

	cosdec = cos (DEC_FLD(indat))
	cosp = cos (PA_ROT(indat))
	sinp = sin (PA_ROT(indat))

	while ( clgcur("coord", wx, wy, wcs, key, command, 32) != EOF ) {

	    if (key == 'q') {
		break
	    }

	    switch (key) {			# NB note TWO switch cases!

# x = cosp
# y = sinp
# translations  XXX there's a lot of excess calc's going on here

		case 'h':
			xoff = DEGTORAD(xyinc/3600.)
			yoff = 0.
			doff =  xoff * cosp - yoff * sinp
			aoff =  xoff * sinp + yoff * cosp
			RA_FLD(indat) = RA_FLD(indat) + aoff / cosdec
			DEC_FLD(indat) = DEC_FLD(indat) + doff
			
		case 'j':
			xoff = 0.
			yoff = DEGTORAD(xyinc/3600.)
			doff =  xoff * cosp - yoff * sinp
			aoff =  xoff * sinp + yoff * cosp
			RA_FLD(indat) = RA_FLD(indat) + aoff / cosdec
			DEC_FLD(indat) = DEC_FLD(indat) + doff
			
		case 'k':
			xoff = 0.
			yoff = -DEGTORAD(xyinc/3600.)
			doff =  xoff * cosp - yoff * sinp
			aoff =  xoff * sinp + yoff * cosp
			RA_FLD(indat) = RA_FLD(indat) + aoff / cosdec
			DEC_FLD(indat) = DEC_FLD(indat) + doff
			
		case 'l':
			xoff = -DEGTORAD(xyinc/3600.)
			yoff = 0.
			doff =  xoff * cosp - yoff * sinp
			aoff =  xoff * sinp + yoff * cosp
			RA_FLD(indat) = RA_FLD(indat) + aoff / cosdec
			DEC_FLD(indat) = DEC_FLD(indat) + doff
			
		case 'p':
			PA_ROT(indat) = PA_ROT(indat) + anginc
			if (PA_ROT(indat) > PI)
				PA_ROT(indat) = PA_ROT(indat) - 2.*PI

		case 'n':
			PA_ROT(indat) = PA_ROT(indat) - anginc
			if (PA_ROT(indat) < -PI)
				PA_ROT(indat) = PA_ROT(indat) + 2.*PI

		case 's':
			do nlist = 1, NLIST_MAX
				call selector (indat, tdat, ntarg, nlist, 2.*(DEF_HLEN(indat)+SLIT_GAP(indat)), psum)	# XXX
			nlist = PRIMARY
			do i = 0, ntarg-1
				call mark_obj (gp, tdat, i, nlist)

			call printf ("Total priorities: %d\n")
				call pargi (psum)

		case 'a':
			i = get_sw_nearest0 (gp, XARCS(tdat,0), YARCS(tdat,0),
				SEL(tdat,0), ntarg, NO, wx, wy, wcs)
			SEL(tdat,i) = YES
			call mark_obj (gp, tdat, i, SAMPL(tdat,i))

		case 'd':
			i = get_sw_nearest0 (gp, XARCS(tdat,0), YARCS(tdat,0),
				SEL(tdat,0), ntarg, YES, wx, wy, wcs)
			SEL(tdat,i) = NO
			call mark_obj (gp, tdat, i, nlist)

		case 'i':
			call amovkr (NO, SEL(tdat,0), ntarg)
			nslit = 0

		case 'g':
call eprintf ("NO projection\n")
			if (sdat != 0)
				call slit_free (sdat, nslit)
			call gen_slits (tdat, ntarg, sdat, nslit, indat)

		case 'z':
			if (get_valr ("zoom_factor:", "%3f", zfact, zfact, -5., 80.) == ERR)
				next
			if (zfact == 0.) {
				FPWIN_X(mdat) = DEF_FLD_CENX
				FPWIN_Y(mdat) = DEF_FLD_CENY
				FPWIN_W(mdat) = DEF_FLD_WIDTH
			} else if (zfact == -1.) {
				FPWIN_X(mdat) = FLIP*-375.-80.
				FPWIN_Y(mdat) = 340.
				FPWIN_W(mdat) = DEF_FLD_WIDTH / 2.4
			} else if (zfact == -2.) {
				FPWIN_X(mdat) = FLIP*-125.-80.
				FPWIN_Y(mdat) = 340.
				FPWIN_W(mdat) = DEF_FLD_WIDTH / 2.4
			} else if (zfact == -3.) {
				FPWIN_X(mdat) = FLIP*125.-80.
				FPWIN_Y(mdat) = 340.
				FPWIN_W(mdat) = DEF_FLD_WIDTH / 2.4
			} else if (zfact == -4.) {
				FPWIN_X(mdat) = FLIP*375.-80.
				FPWIN_Y(mdat) = 340.
				FPWIN_W(mdat) = DEF_FLD_WIDTH / 2.4
			} else if (zfact == -5.) {
				FPWIN_X(mdat) = FLIP*100.-50.
				FPWIN_Y(mdat) = 195.
				FPWIN_W(mdat) = DEF_FLD_WIDTH / 3.
			} else if (zfact > 0.) {
				FPWIN_X(mdat) = wx
				FPWIN_Y(mdat) = wy
				FPWIN_W(mdat) = DEF_FLD_WIDTH / zfact
			}
				

		case ' ':
			i = get_sw_nearest0 (gp, XARCS(tdat,0), YARCS(tdat,0),
				SEL(tdat,0), ntarg, INDEFI, wx, wy, wcs)

# print full data line on text screen, partial on graphics
			call gdeactivate (gp, 0)
			call printf ("%s\n")
				call pargstr (DATLINE(tdat,i))
			call greactivate (gp, 0)
			j = sscan (DATLINE(tdat,i))
			call gargwrd (objname, SZ_ID)
			call printf ("%-s   code=%d   %5.2f_%1s  (x,y)= %.1f %.1f arcsec\n")
				call pargstr (objname)
				call pargi (PCODE(tdat,i))
				call pargr (MAG(tdat,i))
				call pargc (PBAND(tdat,i))
				call pargr (XARCS(tdat,i))
				call pargr (YARCS(tdat,i))


		case 't':
			i = get_sw_nearest0 (gp, XARCS(tdat,0), YARCS(tdat,0),
				SEL(tdat,0), ntarg, INDEFI, wx, wy, wcs)
			j = sscan (DATLINE(tdat,i))
			call gargwrd (objname, SZ_ID)
			if (PCODE(tdat,i) == CODE_GS) {
				PCODE(tdat,i) = CODE_AS
				call mark_obj (gp, tdat, i, nlist)
				call printf ("%-s changed to Alignment\n")
					call pargstr (objname)
			} else if (PCODE(tdat,i) == CODE_AS) {
				PCODE(tdat,i) = CODE_GS
				call mark_obj (gp, tdat, i, nlist)
				call printf ("%-s changed to Guide\n")
					call pargstr (objname)
			} else {
				call printf ("%-s is not a star\n")
					call pargstr (objname)
			}

	    }

	    switch (key) {

		case 'h','j','k','l','p','n':
			cosdec = cos (DEC_FLD(indat))
			cosp = cos (PA_ROT(indat))
			sinp = sin (PA_ROT(indat))
		    call fld2telax (indat)
		    call tel_coords (tdat, ntarg, indat, YES)
#  things have moved; remove any slits, if present:
		    if (sdat != 0)
			call slit_free (sdat, nslit)

	    }

	    switch (key) {

		case 'r','c','h','j','k','l','p','n','f','i','w','g','z':
			call fp_layout (gp, mdat, indat, tdat, ntarg, sdat, nslit, nlist, LBL1, LBL2, LBL3)

		case '.':
			if (speed == FAST) {
				speed = MEDIUM
				call printf (" medium")
			} else if (speed == MEDIUM) {
				speed = SLOW
				call printf (" slow")
			} else if (speed == SLOW) {
				speed = FAST
				call printf (" fast")
			}
			anginc = DEGTORAD(1.) * speed	# TMP -- see above
			xyinc = 6. * speed	# TMP -- see above


		case '?':
		    call gpagefile (gp, KEYSFILE, "simulator cursor commands")

		case 'I':
			call fatal (0, "INTERRUPT")
	    }
	}

	call gclose (gp)
end

#
# GET_OUTLINES: read in the outlines of mask, TV in Focal plane
# 

procedure	get_outlines (mdat)

pointer	mdat				# pointer to Foc Plane struct

char	tchar
int	nact
int	ndx
pointer	fdf

int	line_count()
int	fscan(), nscan()
pointer	open()
begin
	call malloc (mdat, MDATLEN, TY_STRUCT)

	fdf = open (FP_FILE, READ_ONLY, TEXT_FILE)

	nact = line_count (fdf)

	NFPDES(mdat) = nact
	call malloc (PTFPLX(mdat), nact, TY_REAL)
	call malloc (PTFPLY(mdat), nact, TY_REAL)
	call malloc (PTFPLZ(mdat), nact, TY_INT)
	call malloc (PTFPWIN(mdat), NFPWPAR, TY_REAL)

# read in the descriptions (x,y,pen)
	ndx = 0
	while (fscan (fdf) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#' || nscan() == 0) {
			next
		}
		call reset_scan()
		call gargr (FPLX(mdat,ndx))
		call gargr (FPLY(mdat,ndx))
		call gargi (FPLZ(mdat,ndx))
		if (nscan() < 3)
			FPLZ(mdat,ndx) = -1
		else if (nscan() < 2)
			call fatal (0, "Problem with foc_plane file!")
		ndx = ndx + 1
	}

# do ndx = 0, nact-1 {
# call eprintf (" %6.0f %6.0f  %2d\n")
# call pargr (FPLX(mdat,ndx))
# call pargr (FPLY(mdat,ndx))
# call pargi (FPLZ(mdat,ndx))
# }

	call close (fdf)
end

#
# FP_LAYOUT: graphically describe the Focal Plane layout
#
	

procedure	fp_layout (gp, mdat, indat, tdat, ntarg, sdat, nslit, nlist, lab1, lab2, lab3)

pointer	gp
pointer	mdat
pointer	indat
pointer	tdat
int	ntarg
pointer	sdat
int	nslit
int	nlist			# number of the "active" list
char	lab1[ARB], lab2[ARB], lab3[ARB]

int	i
int	pen				# code for "pen" action

char	pafld[SZ_ID]
real	gx1, gx2, gy1, gy2
real	gszx, gszy
real	xarea[5], yarea[5]
real	hbox				# box half-width (arcsec)
real	x0, y0, x, y
real	cosp, sinp
real	szasec, swid			# sz of 1 arcsec, slitwid

real	xback[5], yback[5]		# TMP!! XXX

real	aspect
real	axlen

real	ggetr()
# pointer	gopen()
begin
#	hbox = PL_FWID(limit) / 2.
#	hbox = 1100. / 2.
# get the aspect ratio for the device, so that squares are square.
	aspect = ggetr (gp, "ar")

	hbox = FPWIN_W(mdat) / 2.

	gx1 = FPWIN_X(mdat) - hbox	# PL_XMAG(limit)
	gx2 = FPWIN_X(mdat) + hbox	# PL_XMAG(limit)
	gy1 = FPWIN_Y(mdat) - hbox * aspect * 1.2
	gy2 = FPWIN_Y(mdat) + hbox * aspect * 0.8

	gszx = gx2 - gx1
	gszy = gx2 - gx1
	
# set windows
# XXX All windows, viewports need to be reviewed; aspect currently used poorly
	call gclear (gp)
#	call gsview (gp, (0.5-0.49*aspect), (0.5+0.49*aspect), 0.01, 0.99)
	call gsview (gp, 0.005, 0.995, 0.005, 0.995)
	call gswind (gp, gx1, gx2, gy1, gy2)

# TMP! blacken background: XXX
xback[1] = gx1; yback[1] = gy1
xback[2] = gx1; yback[2] = gy2
xback[3] = gx2; yback[3] = gy2
xback[4] = gx2; yback[4] = gy1
xback[5] = gx1; yback[5] = gy1
	call gseti (gp, G_FACOLOR, BLACK)		# BACKGROUND)
	call gfill (gp, xback, yback, 4, GF_SOLID)


# Draw the FP boundaries (X refers to DEIMOS X-coords)
	call gseti (gp, G_PLCOLOR, YELLOW)
	do i = 0, NFPDES(mdat)-1 {
		pen = FPLZ(mdat,i)
		if (pen == 0) {
			call gamove (gp, FLIP*FPLX(mdat,i), FPLY(mdat,i))
			next
		}
		if (pen > 0) {
			if (pen == 1)
				call gseti (gp, G_PLTYPE, GL_SOLID)
			else if (pen == 2)
				call gseti (gp, G_PLTYPE, GL_DASHED)
			else if (pen == 3)
				call gseti (gp, G_PLTYPE, GL_DOTTED)
		}
		call gadraw (gp, FLIP*FPLX(mdat,i), FPLY(mdat,i))
	}

# Now mark the objects/slits
	call gseti (gp, G_PLTYPE, GL_SOLID)

	do i = 0, ntarg-1
		call mark_obj (gp, tdat, i, nlist)

	if (nslit > 0) {
	    do i = 0, nslit-1
		call mark_slit (gp, sdat, i)
	}



# Auxilliary plots:
# clear area ...
	xarea[1] = 0.01 * gszx + gx1;   yarea[1] = 0.01 * gszy + gy1
	xarea[2] = xarea[1]         ;   yarea[2] = 0.23 * gszy + gy1
	xarea[3] = 0.33 * gszx + gx1;   yarea[3] = yarea[2]
	xarea[4] = xarea[3]         ;   yarea[4] = yarea[1]
	xarea[5] = xarea[1]         ;   yarea[5] = yarea[1]
	call gseti (gp, G_FACOLOR, AUX_COLOR_2)		# BACKGROUND)
	call gfill (gp, xarea, yarea, 4, GF_SOLID)

	call gseti (gp, G_PLTYPE, GL_SOLID)
	call gseti (gp, G_PLCOLOR, CYAN)	# make text match (5)
	call gpline (gp, xarea, yarea, 5)

# put in axes and labels (X refers to display X)
	x0 = 0.05 * gszx + gx1
	y0 = 0.07 * gszy + gy1
	axlen = 0.05 * gszy
	call gamove (gp, x0, y0)
	x =  axlen
	y =  0
	call grdraw (gp, x, y)
	call gtext (gp, x0+1.2*x, y0, lab1, "h=l;v=c;q=h;s=0.8;col=5")
	call gamove (gp, x0, y0)
	x =  0
	y =  axlen
	call grdraw (gp, x, y)
	call gtext (gp, x0, y0+1.2*y, lab2, "h=c;v=c;q=h;s=0.8;col=5")

	call gtext (gp, x0+1.0*axlen, y0-0.3*axlen, lab3, "h=c;v=t;q=h;s=0.9;col=5")


# Draw the compass rose (X currently refers to display X):
	x0 = 0.26 * gszx + gx1
	y0 = 0.10 * gszy + gy1
	axlen = 0.05 * gszy
	call gseti (gp, G_PLCOLOR, RED)
	call gamove (gp, x0, y0)
	cosp = cos (PI-PAR_ANG(indat)+PA_ROT(indat)) 
	sinp = sin (-PAR_ANG(indat)+PA_ROT(indat))
	x =  axlen * cosp / 1.		# PL_XMAG(limit)
	y =  axlen * sinp
	call grdraw (gp, x, y)

	call gseti (gp, G_PLCOLOR, CYAN)
	call gamove (gp, x0, y0)
	cosp = cos (PI-PA_ROT(indat)) 		# PI needed bc this is wrt -X
	sinp = sin (PA_ROT(indat))
	x =  axlen * cosp / 1.		# PL_XMAG(limit)
	y =  axlen * sinp
	call grdraw (gp, x, y)
	call gtext (gp, x0+1.2*x, y0+1.2*y, "N", "h=c;v=c;q=h;s=0.7;col=5")
	call gamove (gp, x0, y0)
	x = -axlen * sinp / 1.		# PL_XMAG(limit)
	y =  axlen * cosp
	call grdraw (gp, x, y)
	call sprintf (pafld, SZ_ID, "PA=%6.1f")
		call pargd (RADTODEG(PA_ROT(indat)))
	call gtext (gp, x0, y0-1.4*axlen, pafld, "h=c;v=t;q=h;s=0.8;col=5")

# Dispersion Diagram
	call gseti (gp, G_PLCOLOR, CYAN)
	szasec = 0.05 * gszy
	swid = szasec * DEF_SLWID(indat)
	y0 = yarea[2] - 0.04 * gszy
	x0 = xarea[1] + 0.02 * gszy
	call gamove (gp, x0, y0+0.5*swid)
	x0 = xarea[3] - 0.02 * gszy
	call gadraw (gp, x0, y0+0.5*swid)
	x0 = xarea[1] + 0.02 * gszy
	call gamove (gp, x0, y0-0.5*swid)
	x0 = xarea[3] - 0.02 * gszy
	call gadraw (gp, x0, y0-0.5*swid)
	x0 = 0.5 * (xarea[1] + xarea[3])
# get parallactic angle, set disk to slit wid
	cosp = cos (PI-PAR_ANG(indat)+PA_ROT(indat)) 
	sinp = sin (-PAR_ANG(indat)+PA_ROT(indat))
	call gseti (gp, G_FACOLOR, BLUE)
	x =  x0 - AD1(indat) * szasec * cosp / 1.		# PL_XMAG(limit)
	y =  y0 - AD1(indat) * szasec * sinp
	call gmark (gp, x, y, GM_CIRCLE+GM_FILL, -swid, -swid)
	call gseti (gp, G_FACOLOR, RED)
	x =  x0 - AD2(indat) * szasec * cosp / 1.		# PL_XMAG(limit)
	y =  y0 - AD2(indat) * szasec * sinp
	call gmark (gp, x, y, GM_CIRCLE+GM_FILL, -swid, -swid)
	call gseti (gp, G_FACOLOR, GREEN)
	call gmark (gp, x0, y0, GM_CIRCLE+GM_FILL, -swid, -swid)
# ... label slit
	x0 = xarea[1] + 0.03 * gszy
	call gtext (gp, x0, y0, "Slit", "h=l;v=c;q=h;s=1.0;col=5")
	x0 = xarea[3] - 0.03 * gszy
	call sprintf (pafld, SZ_ID, "%5.2f''")
		call pargr (DEF_SLWID(indat))
	call gtext (gp, x0, y0, pafld, "h=r;v=c;q=h;s=1.0;col=5")


end

#
# MARK_OBJ: mark an object with relevant details
#

procedure	mark_obj (gp, tdat, i, nlist)

pointer	gp
pointer	tdat
int	i
int	nlist			# number of the "active" list

real	x, y, xa, xb, ya, yb
real	szbx
begin
# Work out the location/size
	x = FLIP * XARCS(tdat,i)
	y = YARCS(tdat,i)
	xa = FLIP * X1(tdat,i)
	xb = FLIP * X2(tdat,i)
	ya = Y1(tdat,i)
	yb = Y2(tdat,i)

# assign appropriate color
	if (STAT(tdat,i) == YES) {
	    if (SAMPL(tdat,i) == nlist)
		call gseti (gp, G_PLCOLOR, GREEN)
	    else if (PCODE(tdat,i) == CODE_AS)		# TMP! XXX
		call gseti (gp, G_PLCOLOR, MAGENTA)	# TMP! XXX
	    else if (PCODE(tdat,i) == CODE_GS)		# TMP! XXX
		call gseti (gp, G_PLCOLOR, MAGENTA)	# TMP! XXX
	    else
		call gseti (gp, G_PLCOLOR, BLUE)
	} else {
		call gseti (gp, G_PLCOLOR, RED)
	}


	if (SEL(tdat,i) == YES)
		call gseti (gp, G_PLCOLOR, WHITE)	# TMP XXX TEST


	if (PCODE(tdat,i) == CODE_AS) {
		szbx = -1. * abs (xa-xb)		# note negative
		call gmark (gp, x, y, GM_BOX, szbx, szbx)
	} else if (PCODE(tdat,i) == CODE_GS) {
		call gmark (gp, x, y, GM_CROSS+GM_PLUS, -5., -5.)
	} else {
		call gamove (gp, xa, ya)
		call gadraw (gp, xb, yb)
		if (SEL(tdat,i) == YES)
		    call gmark (gp, x, y, GM_BOX, -1.*SLWID(tdat,i), -1.*SLWID(tdat,i))	# TMP
		else
		    call gmark (gp, x, y, GM_PLUS, -1.*SLWID(tdat,i), -1.*SLWID(tdat,i))	# TMP
	}
end

#
# MARK_SLIT: mark an object with relevant details
#

procedure	mark_slit (gp, sdat, i)

pointer	gp
pointer	sdat
int	i

real	x, y
real	yoff
begin
	call gseti (gp, G_PLCOLOR, CYAN)

	yoff = 0.5*SLWID(sdat,i)

	x = FLIP * X1(sdat,i)
	y = Y1(sdat,i)
	call gamove (gp, x, y+yoff)
	call gadraw (gp, x, y-yoff)

	x = FLIP * X2(sdat,i)
	y = Y2(sdat,i)
	call gadraw (gp, x, y-yoff)
	call gadraw (gp, x, y+yoff)

	x = FLIP * X1(sdat,i)
	y = Y1(sdat,i)
	call gadraw (gp, x, y+yoff)

end

#
# GET_SW_NEAREST0: get nearest point with appropriate value of an int switch
# "0" indicates zero-indexed. NB!! CAVE!! in this version, x is flipped
# since +X(display) is -X(DEIMOS)
#

int procedure get_sw_nearest0 (gp, xdata, ydata, isw, ndata, sw_val, wx, wy, wcs)

pointer	gp
real 	xdata[ARB], ydata[ARB]
int	isw[ARB]
int	ndata
int	sw_val
real	wx, wy
int	wcs

int	nearest, i
real	ycorr
real	rsq, rsq_min
real	xndc, yndc, xgcndc, ygcndc

real	ggetr()

begin

# need to put in INDEF check 

# Get aspect ratio 
	ycorr = ggetr (gp, "ar")
	if (ycorr == 0.)
		ycorr = 1.

	rsq_min = 2.			# by def'n larger than NDC possible
	nearest = 0

	call gctran (gp, wx, wy, xgcndc, ygcndc, wcs, 0)

	if (sw_val == INDEFI) {
	    do i = 1, ndata {
		call gctran (gp, FLIP*xdata[i], ydata[i], xndc, yndc, wcs, 0)
		rsq = (xndc - xgcndc) ** 2 + ( (yndc - ygcndc) * ycorr) ** 2
		if (rsq < rsq_min) {
			rsq_min = rsq
			nearest = i
		}
	    }
	} else {
	    do i = 1, ndata {
		if (isw[i] != sw_val)
			next
		call gctran (gp, FLIP*xdata[i], ydata[i], xndc, yndc, wcs, 0)
		rsq = (xndc - xgcndc) ** 2 + ( (yndc - ygcndc) * ycorr) ** 2
		if (rsq < rsq_min) {
			rsq_min = rsq
			nearest = i
		}
	    }
	}

	if (nearest != 0)
		call gscur (gp, FLIP*xdata[nearest], ydata[nearest])
	else
		call eprintf ("no appropriate points")
	
	return (nearest - 1)		# Zero indexed
end
