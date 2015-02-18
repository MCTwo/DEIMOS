#
# WRITE_MONGO: produces a MONGO input file
#

include	<math.h>
include	"deimos.h"
include	"dsimulator.h"

define	TVMINMAG 20.5

procedure	write_mongo (fdp, tdat, nt, sdat, nslit, indat, fdtv1, fdtv2, gsndx, ngs, delpa, plotfile)

pointer	fdp
pointer	tdat
int	nt
pointer	sdat
int	nslit
pointer	indat
pointer	fdtv1, fdtv2		# pointers to input files for TV surf map
int	gsndx[ngs]		# Object indices of TV stars, sorted mag
int	ngs
real	delpa
char	plotfile[ARB]		# plotfile name

int	i, j
int	pen
pointer	mdat

char	rootname[SZ_FNAME]	# attempt to contruct meaningful filename

real	xtv, ytv		# Coordinates in the TV system
real	xas, yas		# x,y in arcsec
int	ibrit, imirr		# counter for bright, mirror stars
real	xlab, ylab		# position for TV star info
pointer	tvsf1x, tvsf1y, tvsf2x, tvsf2y	# pointers to surface fits for TV

int	strldx()
real	gseval()
long	clktime()
begin
	i = strldx (".", plotfile)
	if (i == 0)
		i = SZ_FNAME
	call strcpy (plotfile, rootname, i-1)
	
	call fprintf (fdp, "psland %-s.ps\n")
		call pargstr (rootname)
	call fprintf (fdp, "color 1 0 0 0 \n")
	call fprintf (fdp, "margin 0.6 0.1 1.6 0.8\n")
	call fprintf (fdp, "submargin 0. 0.\n")
	call fprintf (fdp, "limits -375.4 375.4  240. -160.\n")
	call fprintf (fdp, "square -1 -1 -1 -1\n")
	do i = 0, nslit-1 {
		call fprintf (fdp, "relocate %8.4f %8.4f\n")
			call pargd (XMM1(sdat,i))
			call pargd (YMM1(sdat,i))
		call fprintf (fdp, "draw   %8.4f %8.4f\n")
			call pargd (XMM2(sdat,i))
			call pargd (YMM2(sdat,i))
		call fprintf (fdp, "draw   %8.4f %8.4f\n")
			call pargd (XMM3(sdat,i))
			call pargd (YMM3(sdat,i))
		call fprintf (fdp, "draw   %8.4f %8.4f\n")
			call pargd (XMM4(sdat,i))
			call pargd (YMM4(sdat,i))
		call fprintf (fdp, "draw   %8.4f %8.4f\n")
			call pargd (XMM1(sdat,i))
			call pargd (YMM1(sdat,i))
	}

	call get_outlines (mdat)
	call fprintf (fdp, "color 2 0 1 0 \n")
	do i = 0, NFPDES(mdat)-1 {
		pen = FPLZ(mdat,i)
		if (pen == 0) {
			call fprintf (fdp, "relocate %8.4f %8.4f\n")
				call pargr (FPLX(mdat,i)*0.7277)
				call pargr (FPLY(mdat,i)*0.7277-ZPT_YM)
			next
		}
		if (pen > 0) {
			if (pen == 1)
				call fprintf (fdp, "ltype 0\n")
			else if (pen == 2)
				call fprintf (fdp, "ltype 2\n")
			else if (pen == 3)
				call fprintf (fdp, "ltype 1\n")
		}
		call fprintf (fdp, "draw   %8.4f %8.4f\n")
			call pargr (FPLX(mdat,i)*0.7277)
			call pargr (FPLY(mdat,i)*0.7277-ZPT_YM)
	}

# TV guide star info -- draw labels first so they lie underneath
	call fprintf (fdp, "ltype 0\n")
	call fprintf (fdp, "ltype 0\n")
	call fprintf (fdp, "expand 0.75\n")

# Now get TV star info:
	call gs_ingest (fdtv1, tvsf1x, tvsf1y)
	call gs_ingest (fdtv2, tvsf2x, tvsf2y)

	xlab = 200.
	ylab = -105.
	ibrit = 0
	imirr = 0
	call fprintf (fdp, "angle 180.\n")	# label upside down
	do j = 1, ngs {
		i = gsndx[j]
		xas = XARCS(tdat,i)
		yas = YARCS(tdat,i)
		if (yas < 176.) {
			xtv = gseval (tvsf1x, xas, yas)
			ytv = gseval (tvsf1y, xas, yas)
			imirr = imirr + 1
		} else {
			xtv = gseval (tvsf2x, xas, yas)
			ytv = gseval (tvsf2y, xas, yas)
		}
		ibrit = ibrit + 1

# make sure at least 2 mirror stars are listed
		if (ibrit > 3 && imirr < 2)
			next

# put label near star
		call fprintf (fdp, "color 4 0.7 0.7 0.7 \n")
		call fprintf (fdp, "rel %6.1f %5.1f \n")
			call pargr (xas*0.7277)
			call pargr (yas*0.7277-ZPT_YM)
		call fprintf (fdp, "draw %6.1f %5.1f \n")
			call pargr (xlab)
			call pargr (ylab)
		call fprintf (fdp, "color 1 0 0 0 \n")
		call fprintf (fdp, "putlabel 4 tv: %4.0f %4.0f\n")
			call pargr (xtv)
			call pargr (ytv)

		ylab = ylab + 18.

# if we reach 5, with at least 2 mirror stars, quit
		if (ibrit >= 5 && imirr >=2)
			break
	}
	call fprintf (fdp, "angle 0.\n")	# return

# now plot the stars themselves
	call fprintf (fdp, "ltype 0\n")
	call fprintf (fdp, "color 3 0 0 1 \n")

	do i = 0, nt-1 {
		if (SEL(tdat,i) == NO)
			next
		if (PCODE(tdat,i) != CODE_GS)
			next
		if (MAG(tdat,i) > TVMINMAG)		# XXX kludge
			next

		call fprintf (fdp, "expand %4f \n ptype 12 0 \n")
			call pargr (0.16*(TVMINMAG-MAG(tdat,i)))

		call fprintf (fdp, "rel %6.1f %5.1f \n dot\n")
			call pargr (XARCS(tdat,i)*0.7277)
			call pargr (YARCS(tdat,i)*0.7277-ZPT_YM)
	}
#call fprintf (fdp, "limits 0 1 0 1 \n")
#call fprintf (fdp, "rel 0 0\ndraw 0 1\ndraw 1 1\ndraw 1 0\ndraw 0 0\n")

	call fprintf (fdp, "margin 1.0 0.9 1.6 1.0\n")
	call fprintf (fdp, "limits 0 1 0 1 \n")
	call fprintf (fdp, "color 1 0 0 0 \n")
#call fprintf (fdp, "rel 0 0\ndraw 0 1\ndraw 1 1\ndraw 1 0\ndraw 0 0\n")

	call fprintf (fdp, "expand 1.4 \n")
	call fprintf (fdp, "rel 0.02 0.95 \n")
	call fprintf (fdp, "putlabel 6 %s\n")
		call pargstr (DESNAME(indat))
	call fprintf (fdp, "rel 0.98 0.95 \n")
	call fprintf (fdp, "putlabel 4 %s\n")
		call pargstr (GUINAME(indat))

	call fprintf (fdp, "expand 1.05 \n")
	call fprintf (fdp, "rel 0.02 0.88 \n")
	call fprintf (fdp, "putlabel 6 RA= %011.2h   Dec= %011.1h \n")
		call pargd (RADTODEG(RA0_FLD(indat))/15.0D0)
		call pargd (RADTODEG(DEC0_FLD(indat)))

	call fprintf (fdp, "rel 0.02 0.84 \n")
	call fprintf (fdp, "putlabel 6 PA(%8.3f)= %8.3f; PA(%8.3f)= %8.3f\n")
		call pargd (STD_EQX(indat))
		call pargd (RADTODEG(PA_ROT(indat)))
		call pargd (EPOCH(indat))
		call pargd (RADTODEG(PA_ROT(indat)+delpa))

	call fprintf (fdp, "expand 0.7 \n")
	call fprintf (fdp, "rel 0.02 0.79 \n")
	call cnvtime (clktime (0), rootname, SZ_FNAME)
	call fprintf (fdp, "putlabel 6 (%s    %s) \n")
		call pargstr (plotfile)
		call pargstr (rootname)

	call fprintf (fdp, "hard \n")

end

