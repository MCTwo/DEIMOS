####### Definitions for DSIMULATOR -- DEIMOS slitmask design

define  FLIP	-1.	# Fudge factor for flipping X axis in graphics

define	ASEC_RAD	206264.8D0

# May want to replace STAT with "OK"; then can use YES/NO unambiguously

define	SZ_ID		16	## TMP? reconcile elsewhere
define	SZ_PBAND	12	# Size of passband string (XXX not well impl.)
define	PRIMARY		 1	# TMP? resolve how to name lists
define	NLIST_MAX	4	# TMP?  maximum number of lists

# Define the struct for the target data
define	TDATLEN		24
define	PTINDEX		Memi[$1]
define	PTRA0		Memi[$1+1]
define	PTDEC0		Memi[$1+2]
define	PTRA		Memi[$1+3]
define	PTDEC		Memi[$1+4]
define	PTPA		Memi[$1+5]
define	PTLEN1		Memi[$1+6]
define	PTLEN2		Memi[$1+7]
define	PTWID		Memi[$1+8]
define	PTPCODE		Memi[$1+9]
define	PTXARCS		Memi[$1+10]
define	PTYARCS		Memi[$1+11]
define	PTRELPA		Memi[$1+12]
define	PTX1		Memi[$1+13]
define	PTY1		Memi[$1+14]
define	PTX2		Memi[$1+15]
define	PTY2		Memi[$1+16]
define	PTSTAT		Memi[$1+17]
define	PTMAG		Memi[$1+18]
define	PTPBAND		Memi[$1+19]
define	PTSAMPL		Memi[$1+20]
define	PTSEL		Memi[$1+21]
define	PTSLNDX		Memi[$1+22]
define	PTLINE		Memi[$1+23]

# Note that the number of targets can also be described here.

define	INDEX		Memi[PTINDEX($1)+$2]		# index
define	RA0		Memd[PTRA0($1)+$2]		# unrefracted RA (rad)
define	DEC0		Memd[PTDEC0($1)+$2]		# unrefracted Dec (rad)
define	RA		Memd[PTRA($1)+$2]		# refracted RA (rad)
define	DEC		Memd[PTDEC($1)+$2]		# refracted Dec (rad)
define	PA		Memr[PTPA($1)+$2]		# PA on sky (rad)
define	LEN1		Memr[PTLEN1($1)+$2]		# length 1
define	LEN2		Memr[PTLEN2($1)+$2]		# length 2
define	SLWID		Memr[PTWID($1)+$2]		# width
define	PCODE		Memi[PTPCODE($1)+$2]		# pcode
define	XARCS		Memr[PTXARCS($1)+$2]		# X in arcsec
define	YARCS		Memr[PTYARCS($1)+$2]		# Y in arcsec
define	RELPA		Memr[PTRELPA($1)+$2]		# Rel PA wrt tel (rad)
define	X1		Memr[PTX1($1)+$2]		# X1 in arcsec
define	Y1		Memr[PTY1($1)+$2]		# Y1 in arcsec
define	X2		Memr[PTX2($1)+$2]		# X2 in arcsec
define	Y2		Memr[PTY2($1)+$2]		# Y2 in arcsec
define	STAT		Memi[PTSTAT($1)+$2]		# OK for selection?
define	MAG		Memr[PTMAG($1)+$2]		# Magnitude
define	PBAND		Memc[PTPBAND($1)+$2]		# Passband for magnitude
define	SAMPL		Memi[PTSAMPL($1)+$2]		# sample code
define	SEL		Memi[PTSEL($1)+$2]		# Selected?
define	SLNDX		Memi[PTSLNDX($1)+$2]		# index to slit
define	DATLINE		Memc[PTLINE($1)+$2*SZ_LINE]	# full data line


# We want to bundle other items here:
# ra_tel, dec_tel, eqnx_std, pa_rot, ha0, temp, pressure

define	NINDAT		3
define	PTTELDAT	Memi[$1]
define	PTDEFDAT	Memi[$1+1]
define	PTMSKDAT	Memi[$1+2]

define	NTELPAR		21
define	RA_TEL		Memd[PTTELDAT($1)]	# RA of tel. axis (rad)
define	DEC_TEL		Memd[PTTELDAT($1)+1]	# Dec of tel. axis (rad)
define	PA_ROT		Memd[PTTELDAT($1)+2]	# PA of rotator (rad)
define	RA0_FLD		Memd[PTTELDAT($1)+3]	# RA of fld (rad)
define	DEC0_FLD	Memd[PTTELDAT($1)+4]	# dec of fld (rad)
define	RA_FLD		Memd[PTTELDAT($1)+5]	# refracted RA of fld (rad)
define	DEC_FLD		Memd[PTTELDAT($1)+6]	# refracted dec of fld (rad)
define	HA_FLD		Memd[PTTELDAT($1)+7]	# HA of fld (rad)
define	STD_EQX		Memd[PTTELDAT($1)+8]	# standard equinox (yr)
define	TEMP		Memd[PTTELDAT($1)+9]	# temp (C)
define	PRES		Memd[PTTELDAT($1)+10]	# atm. press (mm Hg)
define	PAR_ANG		Memd[PTTELDAT($1)+11]	# parallactic angle (rad)
define	WAVER		Memd[PTTELDAT($1)+12]	# wavel. for refract. (microns)
define	REF1		Memd[PTTELDAT($1)+13]	# Refract coeff, 1st order
define	REF3		Memd[PTTELDAT($1)+14]	# Refract coeff, 3rd order
define	WAVEMN		Memd[PTTELDAT($1)+15]	# min wavel. for refr (um)
define	WAVEMX		Memd[PTTELDAT($1)+16]	# max wavel. for refr (um)
define	AD1		Memd[PTTELDAT($1)+17]	# Atm disp. for WAVEMN (asec)
define	AD2		Memd[PTTELDAT($1)+18]	# Atm disp. for WAVEMN (asec)
define	AMASS		Memd[PTTELDAT($1)+19]	# Airmass
define	EPOCH		Memd[PTTELDAT($1)+20]	# Epoch of planned observation

define	NDEFPAR		6
define	SLIT_GAP	Memr[PTDEFDAT($1)]	# slit separation (asec)
define	DEF_HLEN	Memr[PTDEFDAT($1)+1]	# min slit 1/2-length (asec)
define	DEF_BOXR	Memr[PTDEFDAT($1)+2]	# box radius (asec)
define	DEF_SLWID	Memr[PTDEFDAT($1)+3]	# default slit width (asec)
define	PROJ_LEN 	Memi[PTDEFDAT($1)+4]	# proj. to preserve slit length
define	ADJ_LEN 	Memi[PTDEFDAT($1)+5]	# Adj. len for no overlap
						# NB LAST FUDGED TY_REAL?=TY_INT

define	NMSKPAR		9
define	DESNAME		Memc[PTMSKDAT($1)]		# Design name
define	DESAUTH		Memc[PTMSKDAT($1)+1*SZ_LINE]	# Design Author
define	DESCREAT	Memc[PTMSKDAT($1)+2*SZ_LINE]	# Design Creation Tool
define	PROJNAME	Memc[PTMSKDAT($1)+3*SZ_LINE]	# Project name
define	INSTRUME	Memc[PTMSKDAT($1)+4*SZ_LINE]	# Instrument
define	BLUNAME		Memc[PTMSKDAT($1)+5*SZ_LINE]	# Blueprint name
define	GUINAME		Memc[PTMSKDAT($1)+6*SZ_LINE]	# Name to appear on GUI
define	BLUOBSR		Memc[PTMSKDAT($1)+7*SZ_LINE]	# Observer
define	USEDATE		Memc[PTMSKDAT($1)+8*SZ_LINE]	# Date of intended use

# NB We ASSUME that PA_ROT is the on-axis PA of the rotator, || +X axis

# Definitions for slits: Note that many are same as for targets and use th
# same definition

define	SDATLEN		19
#define	PTINDEX		Memi[$1]	# defined above with tdat
#define	PTRA0		Memi[$1+1]
#define	PTDEC0		Memi[$1+2]
#define	PTRA		Memi[$1+3]	# defined above with tdat
#define	PTDEC		Memi[$1+4]	# defined above with tdat
#define	PTPA		Memi[$1+5]	# defined above with tdat
#define	PTLEN1		Memi[$1+6]	# defined above with tdat
#define	PTLEN2		Memi[$1+7]	# defined above with tdat
#define	PTWID		Memi[$1+8]	# defined above with tdat
#define	PTPCODE		Memi[$1+9]	# defined above with tdat
#define	PTXARCS		Memi[$1+10]	# defined above with tdat
#define	PTYARCS		Memi[$1+11]	# defined above with tdat
#define	PTRELPA		Memi[$1+12]	# defined above with tdat
#define	PTX1		Memi[$1+13]	# defined above with tdat
#define	PTY1		Memi[$1+14]	# defined above with tdat
#define	PTX2		Memi[$1+15]	# defined above with tdat
#define	PTY2		Memi[$1+16]	# defined above with tdat
#define	PTSTAT		Memi[$1+17]	# defined above with tdat
define	PTSCOOR		Memi[$1+18]


define	NSCOOR		12			# number of bundled elements

#define	X1		Memd[PTSCOOR($1)+$2*NSCOOR]	# double? # In asec
#define	Y1		Memd[PTSCOOR($1)+$2*NSCOOR+1]	# double? # In asec
#define	X2		Memd[PTSCOOR($1)+$2*NSCOOR+2]	# double? # In asec
#define	Y2		Memd[PTSCOOR($1)+$2*NSCOOR+3]	# double? # In asec
define	XMM1		Memd[PTSCOOR($1)+$2*NSCOOR+4]	# double? # In mm
define	YMM1		Memd[PTSCOOR($1)+$2*NSCOOR+5]	# double? # In mm
define	XMM2		Memd[PTSCOOR($1)+$2*NSCOOR+6]	# double? # In mm
define	YMM2		Memd[PTSCOOR($1)+$2*NSCOOR+7]	# double? # In mm
define	XMM3		Memd[PTSCOOR($1)+$2*NSCOOR+8]	# double? # In mm
define	YMM3		Memd[PTSCOOR($1)+$2*NSCOOR+9]	# double? # In mm
define	XMM4		Memd[PTSCOOR($1)+$2*NSCOOR+10]	# double? # In mm
define	YMM4		Memd[PTSCOOR($1)+$2*NSCOOR+11]	# double? # In mm



define	KEYSFILE	"deimos$lib/keys/dsimulator.keys"

define	RDB_MAP_FILE	"deimos$lib/RDBmap.hdu"

define	TV_MIRR_MAP	"deimos$mappings/tvmirr.map"	# Arcsec to TV pix
define	TV_MASK_MAP	"deimos$mappings/tvmask.map"	# Arcsec to TV pix
define	XPROJ_MAP	"deimos$mappings/xproj.map"	# "X-distortion"


#############################################################################
# FOCAL PLANE OUTLINES
# This is pure drawing only, and is described by FP_FILE below.
# Need file of (x,y,pen) where pen 0=move,1=solid,2=dashed,3=dotted
# In principle, by simply changing the file we should end up with
# a new instrument description


define	MDATLEN		5
define	NFPDES		Memi[$1]
define	PTFPLX		Memi[$1+1]
define	PTFPLY		Memi[$1+2]
define	PTFPLZ		Memi[$1+3]		# Z is the pen action
define	PTFPWIN		Memi[$1+4]		# bundle windowing params here

define	FPLX		Memr[PTFPLX($1)+$2]
define	FPLY		Memr[PTFPLY($1)+$2]
define	FPLZ		Memi[PTFPLZ($1)+$2]	# Z is the pen action
define	NFPWPAR	3
define	FPWIN_X		Memr[PTFPWIN($1)]	# Window center, X
define	FPWIN_Y		Memr[PTFPWIN($1)+1]	# Window center, Y
define	FPWIN_W		Memr[PTFPWIN($1)+2]	# Window width (X)

#############################################################################

# INSTRUMENT SPECIFIC DEFINITIONS
# Adapting to new instrument includes not only changes here but in check_stat()
# several items needed here for general case:
# - inversion? Attempt to do with SENSX,SENSY
# - axis of dispersion?
# -- we could simply adopt that X will always be slit direction, y-dispersion

# This is the focal-plane outline in arcsec on sky. NB follows same convention
# of N=x-axis, E=y-axis

define	FP_FILE		"deimos$src/foc_plane.dat"

define	SENSX		-1.		# either 1 or -1	# XXX unused?
define	SENSY		 1.		# either 1 or -1	# XXX unused?
define	FLDCEN_X	0.		# Offset (asec) between tel & fld
define	FLDCEN_Y	270.		# Offset (asec) between tel & fld
define	XLOW_LIM	-498.		# lowest X value for slits  XXX
define	XUPP_LIM	 498.		# highest X value for slits XXX
define	YMSKMIN		187.3		# 0.5 arcsec more than "soft" limit (XX)
define	YMSKMAX		479.1		# 0.5 arcsec less than "soft" limit (XX)
define	YCAMCEN		700.		# camera center, in arcsec (11.66')
define	RADVIGN		302.		# radius of cam vign., asec (8.62")

define	GAP1CEN		-254.4		# CCD gap in arcsec -- center
define	GAP1HWD		   5.2		# CCD gap in arcsec -- half-wid
define	GAP2CEN		   0.0		# CCD gap in arcsec -- center
define	GAP2HWD		   5.2		# CCD gap in arcsec -- half-wid
define	GAP3CEN		 254.4		# CCD gap in arcsec -- center
define	GAP3HWD		   5.2		# CCD gap in arcsec -- half-wid

#############################################################################
