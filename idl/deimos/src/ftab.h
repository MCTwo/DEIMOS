####### Definitions for FITS table writing utils
define	FT_CHAR		1	# character-style format
define	FT_FIXED	2	# fixed-style format (int)
define	FT_FLOAT	3	# floating-style format (real, double)

define	TFMTLEN		7
define	TTYPLEN		23
define	TUNILEN		23
define	TNULLEN		69

define	KDATLEN		11
define	TFIELD		Memi[$1]
define	TWUSED		Memi[$1+1]
define	CHARCNT		Memi[$1+2]
define	TBCOLPT		Memi[$1+3]
define	TFORMPT		Memi[$1+4]
define	TTYPEPT		Memi[$1+5]
define	TUNITPT		Memi[$1+6]
define	TNULLPT		Memi[$1+7]
define	TROWPT		Memi[$1+8]
define	THDROK		Memi[$1+9]
define	TMAXCOL		Memi[$1+10]
define	TBCOL		Memi[TBCOLPT($1)+$2]
define	TFORM		Memc[TFORMPT($1)+$2*TFMTLEN]
define	TTYPE		Memc[TTYPEPT($1)+$2*TTYPLEN]
define	TUNIT		Memc[TUNITPT($1)+$2*TUNILEN]
define	TNULL		Memc[TNULLPT($1)+$2*TNULLEN]
define	BUFROW		Memc[TROWPT($1)]

#############################################################################
