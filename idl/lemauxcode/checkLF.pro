pro checkLF

	readcol, 'deep2/Nzlog4.65_0.1_F.dat', format='D,D,D,D', lstar, phistar, N, extra	;B,D,& F are alpha=-1.6, B is a different limiting luminosity, D & F are two different instances of the same limiting luminosity

	readcol, 'BL16a4.65_0.1.dat', format='D,D,D,D', lstar2, phistar2, N2, extra2

	quotient = N2/N
	print, N2[0], N[0], N2[0]/N[0], quotient[0]

	print, 'Mean difference in predicted galaxies between run 1 and run 2 is', mean(quotient)
	print, 'The rms is', stddev(quotient)
end
