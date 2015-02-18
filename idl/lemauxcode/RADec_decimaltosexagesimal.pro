pro RADec_decimaltosexagesimal, RA, Dec

	RAsexh = floor(RA/360.*24)
	RAsexm = floor((RA/360*24 - RAsexh)*60)
	RAsexs = ((RA/360*24 - floor(RAsexh))*60 - RAsexm)*60

	Decsexdeg = floor(Dec)
	Decsexam = floor((Dec-Decsexdeg)*60)
	Decsexas = ((Dec-Decsexdeg)*60 - Decsexam)*60

	print, 'The original Right Ascension of ', RA, 'transforms to ', RAsexh, ':', RAsexm, ':', RAsexs, format='(a32, d7.3, a15, i2.2, a1, i2.2, a1, d6.3)'
	print, 'The original Declination of ', dec, ' transforms to ', Decsexdeg, ':', Decsexam, ':', Decsexas, format='(a28, d7.4, a15, i2.2, a1, i2.2, a1, d6.3)'

end
