#{ Package DEIMOS is Drew Phillips' deimos software


print ("\n DEIMOS (ver.0b1 for IRAF 2.11.3) -- Software in development -- User assumes risk\n")

cl < "deimos$lib/zzsetenv.def"
package	deimos, bin = deimosbin$

task	dsimulator,
	refl,
	trace,
	qtrmap,
	qmodel,
	qrevmod 	= "deimos$src/x_deimos.e"


clbye()
