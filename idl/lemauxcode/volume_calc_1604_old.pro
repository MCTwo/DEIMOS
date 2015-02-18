pro volume_calc_1604_old, incat

d2 = getenv('D2_RESULTS')	

readcol, incat, format='A,A,A, D, D, D', LFCID, mask, slit, rband, iband, zband
ncol = n_elements(LFCID)
	
lambdaarr = findgen(11000)*0.33 + 6000.
numberarray = fltarr(11000)

	for i=0, ncol-1 do begin

	spec = d2 + 'sc1604/' + mask[i] + '/*/spec1d.' + mask[i] + '.' + slit[i] + '.' + LFCID[i] + '.fits'
	specblue = mrdfits(spec,3)
	specred = mrdfits(spec,4)
	print, spec	
	bluegood = where(specblue.spec ne 0)
	redgood = where(specred.spec ne 0) 

	if bluegood[0] ne -1 then specbluegood = specblue.lambda[bluegood]
	if redgood[0] ne -1 then specredgood = specred.lambda[redgood]
	
		for j=0, n_elements(lambdaarr) - 2 do begin
	
			blue = where(specbluegood ge lambdaarr[j] and specbluegood le lambdaarr[j+1])
			red = where(specredgood ge lambdaarr[j] and specredgood le lambdaarr[j+1])
			if red[0] eq -1 then numberarray[j] = numberarray[j] else numberarray[j] = n_elements(red) + numberarray[j]
			if blue[0] eq -1 then numberarray[j] = numberarray[j] else numberarray[j] = n_elements(blue) + numberarray[j] 
		endfor

	endfor
	
openw, lun, 'lambdahist.dat', /get_lun

	for i=0, n_elements(lambdaarr)-1 do begin
	printf, lun, lambdaarr[i], '    ', numberarray[i]
	endfor

free_lun, lun

end
	
