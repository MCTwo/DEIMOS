####### Definitions for DEIMOS MOSAIC

define	SZ_ID		16	## TMP? reconcile elsewhere

# Define the struct for the mosaic data
define	MMAPLEN		22
define	PTHDU1		Memi[$1]
define	PTHDU2		Memi[$1+1]

define	PTCPX1		Memi[$1+2]	# 1st X data pix, HDU1
define	PTC0X1		Memi[$1+3]	# X Val at CP, HDU1
define  PTC1X1		Memi[$1+4]	# Fin X for HDU1
define	PTCPY1		Memi[$1+5]	# 1st Y data pix, HDU1
define	PTC0Y1		Memi[$1+6]	# Y Val at CP, HDU1
define  PTC1Y1		Memi[$1+7]	# Fin Y for HDU1
define	PTCPX2		Memi[$1+8]	# 1st X data pix, HDU2
define	PTC0X2		Memi[$1+9]	# X Val at CP, HDU2
define  PTC1X2		Memi[$1+10]	# Fin X for HDU2
define	PTCPY2		Memi[$1+11]	# 1st Y data pix, HDU2
define	PTC0Y2		Memi[$1+12]	# Y Val at CP, HDU2
define  PTC1Y2		Memi[$1+13]	# Fin Y for HDU2
define  PTGF1		Memi[$1+14]	# Gain factor for HDU1
define  PTBL1		Memi[$1+15]	# Typ. bias lev for HDU1
define  PTGF2		Memi[$1+16]	# Gain factor for HDU2
define  PTBL2		Memi[$1+17]	# Typ. bias lev for HDU2

define	PRECOL		Memi[$1+18]	# Number of Pre columns
define	POSTCOL		Memi[$1+19]	# Number of Post (overscan) Columns
define	IMB_USED	Memi[$1+20]	# Has ImBuf been allocated (used)?
define	IMB_BUFF	Memi[$1+21]	# Pointer to ImBuf

# Note that -1 accounts for 1-indexed CCDs

define	HDU1		Memi[PTHDU1($1)+$2-1]		# image pointer HDU1
define	HDU2		Memi[PTHDU2($1)+$2-1]		# image pointer HDU2
define	CPX1		Memi[PTCPX1($1)+$2-1]		# 1st X data pix, HDU1
define	C0X1		Memi[PTC0X1($1)+$2-1]		# X Val at CP, HDU1
define  C1X1		Memi[PTC1X1($1)+$2-1]		# Fin X for HDU1
define	CPY1		Memi[PTCPY1($1)+$2-1]		# 1st Y data pix, HDU1
define	C0Y1		Memi[PTC0Y1($1)+$2-1]		# Y Val at CP, HDU1
define  C1Y1		Memi[PTC1Y1($1)+$2-1]		# Fin Y for HDU1
define	CPX2		Memi[PTCPX2($1)+$2-1]		# 1st X data pix, HDU2
define	C0X2		Memi[PTC0X2($1)+$2-1]		# X Val at CP, HDU2
define  C1X2		Memi[PTC1X2($1)+$2-1]		# Fin X for HDU2
define	CPY2		Memi[PTCPY2($1)+$2-1]		# 1st Y data pix, HDU2
define	C0Y2		Memi[PTC0Y2($1)+$2-1]		# Y Val at CP, HDU2
define  C1Y2		Memi[PTC1Y2($1)+$2-1]		# Fin Y for HDU2
define  GF1		Memr[PTGF1($1)+$2-1]		# Gain factor for HDU1
define  BL1		Memr[PTBL1($1)+$2-1]		# Typ. bias lev for HDU1
define  GF2		Memr[PTGF2($1)+$2-1]		# Gain factor for HDU2
define  BL2		Memr[PTBL2($1)+$2-1]		# Typ. bias lev for HDU2
