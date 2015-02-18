pro volume_calc_1604, incat

t1 = systime(1) 
d2 = getenv('D2_RESULTS')	

readcol, incat, format='A,A,A, D, D, D', LFCID, mask, slit, rband, iband, zband
ncol = n_elements(LFCID)
	
lambdaarr = findgen(6600)*0.5 + 6100.
numberarray = fltarr(6600)

	for i=0, ncol-1 do begin

	spec = d2 + 'sc1604/' + mask[i] + '/*/spec1d.' + mask[i] + '.' + slit[i] + '.fits'
	bluefits = d2 + 'sc1604/' + mask[i] + '/*/slit.' + mask[i] + '.' + slit[i] + 'B.fits.gz'
	redfits = d2 + 'sc1604/' + mask[i] + '/*/slit.' + mask[i] + '.' + slit[i] + 'R.fits.gz'
	print, spec	

	specblue = mrdfits(bluefits,1, /silent)
	specred = mrdfits(redfits,1, /silent)

	; OLD, not using this anymore, using the 2d wavelength solution instead
	bluehead = headfits(bluefits, ext=1)
	redhead = headfits(redfits,ext=1)
	bluelaminit = sxpar(bluehead, 'CRVAL1')	;getting wavelength information of the 2d flux array, this is waveleneght of the 1st pixel in the blue slit
	redlamlimit = sxpar(redhead, 'CRVAL1')
	dlamblue = sxpar(bluehead, 'CD1_1')	;the wavelength spacing on the blue side
	dlamred = sxpar(redhead, 'CD1_1')
		
	;removing any pixels that have zero flux, mostly to mask the ends of the slit in the dispersion direction, since they don't really contribute anything to the volume. This also removes any bad pixels and chgarge traps, hopefully does not have any wavelength dependance, but maybe important if it does. In any case these pixels don't contribute anything to the volume and should be taken out. The lambda0 tag on the slit files is the initial wavelength from the bottom of the slit and dlambda0 is a 2d array of corrections for the pixel in the kth row on the jth lambda0 value (column)
	bluelam = specblue.flux*0
	redlam = specred.flux*0
	
		for j=0, 4096-1 do begin
		bluelam[j,*] = specblue.dlambda[j,*] + specblue.lambda0[j]
		redlam[j,*] = specred.dlambda[j,*] + specred.lambda0[j]
		endfor	
	badblue = where(specblue.flux eq 0)
	badred = where(specred.flux eq 0)
	
	bluelam[badblue] = 0.	;trivializing the pixels in the lambda array where the flux array is zero	
	redlam[badred] = 0.	

	blueslitsize = size(specblue.flux, /dimensions)
	redslitsize = size(specred.flux, /dimensions)	

	bluespatialsize = blueslitsize[1]	;only getting the spatial sizes, spectral dimensions is always 4096 pixls for both
	redspatialsize = redslitsize[1]	
	
	
		for j=0, n_elements(lambdaarr) - 2 do begin
	
			blue = where(bluelam ge lambdaarr[j] and bluelam le lambdaarr[j+1])
			red = where(redlam ge lambdaarr[j] and redlam le lambdaarr[j+1])
			if red[0] eq -1 then numberarray[j] = numberarray[j] else numberarray[j] = n_elements(red) + numberarray[j]
			if blue[0] eq -1 then numberarray[j] = numberarray[j] else numberarray[j] = n_elements(blue) + numberarray[j] 
		endfor

	endfor
	
openw, lun, 'lambdahist.dat', /get_lun

	for i=0, n_elements(lambdaarr)-1 do begin
	printf, lun, lambdaarr[i], '    ', numberarray[i]
	endfor

free_lun, lun

print
print, systime()+'  Done with calculation!'
print
print, 'Total time elapsed: ', (systime(1)-t1)/3600., ' hours.'

end


	
