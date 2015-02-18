pro k_correct, filter, filter2, incat, outcat, slitloss
	
d2 = getenv('D2_RESULTS')
;incat = d2 + 'sc1604/LFCfluxcalib/' + incat
;outcat = d2 + 'sc1604/LFCfluxcalib/' + outcat
readcol, incat, format='A,A,A,A,D,D,D,D,D,D,D,D,D', LFCID, specID, mask, slit, RA, Dec, rband, iband, zband, F606W, F814W, z, q 
ncol = n_elements(LFCID)

readcol, filter + '.dat', format = 'D,D', lam1, through, /SILENT
minlam = min(lam1)       ;min wavelength of the filter transmission curve
maxlam = max(lam1)	

readcol, filter2 + '.dat', format = 'D,D', lam2, through2, /SILENT 

openw, lun, outcat, /get_lun 
	
	for i=0, ncol-1 do begin

	spec = d2 + 'sc1604/' + mask[i] + '/*/spec1d.' + mask[i] + '.' + slit[i] + '.' + specID[i] + '.fits'
	h = headfits(spec, ext=1)
        exptime = sxpar(h,'EXPTIME')
	
	ss1d=fill_gap(spec,/tweak,/horne,/telluric,/silent,header=header) 	; do the fill_gap stuff, remove A/B band and interpolate over the CCD gap to make one spectrum that contains both the red and blue CCDs with Horne extraction, changed to boxcar extraction (rather than Horne) to minimize aperture effects, response correction is also done here, check fix_response.pro to make sure it is coming out in proper form 
	wave=ss1d.lambda
       	airtovac,wave	; tweak the wavlenegths for the index of refraction of air
	ss1d.lambda=wave
	
	spec_in=ss1d.spec
     	lambda_in=ss1d.lambda
     	ivar_in=ss1d.ivar


	lam2obsframe = lam2*(1+z[i])
	minlam2 = min(lam2obsframe)
	maxlam2 = max(lam2obsframe)

	;print, minlam2, maxlam2, min(lambda_in), max(lambda_in)

	print, spec	
	constfluxcorr = 2e-08/(3600.*!pi*(998./2)^2) 	;2*10^-8 erg*Ang/e-, 3600 s per frame (flux units in spec1d are in counts per hour from headers), pi*498 cm ^2 effective Keck II aperture
	;ALREADY DONE W/ fill_gap corrfactor = deimos_correction(ss1d.lambda)	;deimos response correction as a function of wavelength, multiplicative factor
	corr_spec_in = fltarr(n_elements(ss1d.spec))

		for j=0, n_elements(ss1d.spec)-1 do begin
		corr_spec_in[j] = spec_in[j]*constfluxcorr/((1-slitloss)*lambda_in[j])	;get a flux in units of ergs/s/pixel/cm^2 see http://www.ucolick.org/~npk/dokuwiki/doku.php?id=deep2_flux 
		endfor
	print, minlam, maxlam

	dummy = where(lambda_in ge minlam and lambda_in le maxlam)      ;check where the actual spectrum overlaps the filter transmission curve
        arraylength = n_elements(dummy)
        maxlambda = lambda_in[dummy[arraylength-1]]
        minlambda = lambda_in[dummy[0]]
	througharray = generate_interpol_filtercurve_forKcorr(filter, arraylength, minlambda, maxlambda)	;interpolate values so that the interpolated grid matches up with the spectrum grid, both in intial wavelength values and spacing

	deltalamband = fltarr(n_elements(dummy))
	meanlam = fltarr(n_elements(dummy))
        
	        for j=0, n_elements(dummy)-2 do begin
                deltalamband[j] = abs(lambda_in[dummy[j]]-lambda_in[dummy[j+1]])   ;the frequency interval that the filter throughput covers, needed to get a flux density for converting to AB magnitudes 
                meanlam[j] = (lambda_in[dummy[j]] + lambda_in[dummy[j+1]])/2	;get the effective wavelength in each bin
		endfor
	
	dumdums = where(deltalamband eq 0)               ;clean up last element
        deltalamband[dumdums] = mean(deltalamband)
        dumdumss = where(meanlam eq 0)
        meanlam[dumdumss] = max(meanlam)

	;now do this all again for the 2nd band
	dummy2 = where(lambda_in ge minlam2 and lambda_in le maxlam2)
	arraylength2 = n_elements(dummy2)
	maxlambda2 = lambda_in[dummy2[arraylength2-1]]
	minlambda2 = lambda_in[dummy2[0]]
	througharray2 = generate_interpol_filtercurve_forKcorr_z(filter2, z[i], arraylength2, minlambda2, maxlambda2)
	deltalamband2 = fltarr(n_elements(dummy2))
	meanlam2 = fltarr(n_elements(dummy2))

		for j=0, n_elements(dummy2)-2 do begin
		deltalamband2[j] = abs(lambda_in[dummy2[j]]-lambda_in[dummy2[j+1]])
		meanlam2[j] = (lambda_in[dummy2[j]] + lambda_in[dummy2[j+1]])/2
		endfor

	dumdums2 = where(deltalamband2 eq 0)               ;clean up last element
        deltalamband2[dumdums2] = mean(deltalamband2)
        dumdumss2 = where(meanlam2 eq 0)
        meanlam2[dumdumss2] = max(meanlam2)
	
	bandfluxnum = fltarr(n_elements(dummy))
	bandfluxdenom = fltarr(n_elements(dummy))

	        for j=0, n_elements(dummy)-1 do begin
                bandfluxnum[j] = througharray[j]*corr_spec_in[dummy[j]]*deltalamband[j]/deltalamband[j]*meanlam[j]/(3.0D+18)	;this is pedantic see Fukugita 1995 eqn 9, F_nu is corr_spec_in/deltafband
		bandfluxdenom[j] = througharray[j]*deltalamband[j]/meanlam[j]
                endfor

	;and again for the 2nd band
	
	bandfluxnum2 = fltarr(n_elements(dummy2))
        bandfluxdenom2 = fltarr(n_elements(dummy2))

                for j=0, n_elements(dummy2)-1 do begin
                bandfluxnum2[j] = througharray2[j]*corr_spec_in[dummy2[j]]*deltalamband2[j]/deltalamband2[j]*meanlam2[j]/(3.0D+18)    ;this is pedantic see Fukugita 1995 eqn 9, F_nu is corr_spec_in/deltafband
                bandfluxdenom2[j] = througharray2[j]*deltalamband2[j]/meanlam2[j]
                endfor

	totbandfluxdens = total(bandfluxnum)/total(bandfluxdenom)
	totbandfluxdens2 = total(bandfluxnum2)/total(bandfluxdenom2)

	;splot, lambda_in, spec_in	

	
	ZP = 48.60	;photometric zero point for AB system
	extinction = -0.1	;approximate extinction in magnitudes/airmass at Mauna Kea, see http://www.cfht.hawaii.edu/Instruments/ObservatoryManual/CFHT_ObservatoryManual_(Sec_2).html	
	airmass = sxpar(h,'AIRMASS')
	filterABmag = -ZP - 2.5*alog10(totbandfluxdens) + airmass*extinction
	filter2ABmag = -ZP - 2.5*alog10(totbandfluxdens2) + airmass*extinction
	print, iband[i], filterABmag, filter2ABmag

	;help, LFCID, specID, RA, Dec, mask, slit, rband, iband, zband, F606W, F814W, filterABmag, filter2ABmag, /struct	
	printf, lun, LFCID[i], '   ', RA[i], '   ', Dec[i], '   ', z[i], '   ', mask[i], '   ', slit[i], '   ', rband[i], '   ', iband[i], '   ', zband[i], '   ', F606W[i], '   ', F814W[i], '   ', filterABmag, '   ', filter2ABmag, format ='(a15, a3, d11.7, a3, d10.7, a3, d10.7, a3, a6, a3, a4, a3, d8.5, a3, d8.5, a3, d8.5, a3, d8.5, a3, d8.5, a3, d8.5, a3, d8.5, a3, d8.5)'
	endfor
	
free_lun, lun
end
		
