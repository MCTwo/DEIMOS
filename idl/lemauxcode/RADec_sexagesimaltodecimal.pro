pro RADec_sexagesimaltodecimal, RAh, RAm, RAs, dech, decm, decs

	RAdeg = ((((RAs/60. )+ RAm)/60.)+RAh)/24.*360.

	Decdeg = ((((decs/60. )+ decm)/60.)+dech)

	if decm lt 10 then decm = '0' + string(decm)
	if decs lt 10 then decs = '0' + string(decs)
	if RAm lt 10 then RAm = '0' + string(RAm)
	if RAs lt 10 then RAs = '0' + string(RAs)

	print, 'The original Right Ascension of ', strcompress(string(RAh),/REMOVE_ALL), ':', strcompress(string(RAm),/REMOVE_ALL), ':', strcompress(string(RAs),/REMOVE_ALL),  ' transforms to ', strcompress(string(RADeg), /REMOVE_ALL)
	print, 'The original Declination of ', strcompress(string(dech),/REMOVE_ALL), ':', strcompress(string(decm),/REMOVE_ALL), ':', strcompress(string(Decs),/REMOVE_ALL),  ' transforms to ', strcompress(string(Decdeg), /REMOVE_ALL) 

end
