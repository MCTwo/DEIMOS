pro fluxcalibcheck, filter, incat, outcat
	
d2 = getenv('D2_RESULTS')
;incat = d2 + 'sc1604/LFCfluxcalib/' + incat
;outcat = d2 + 'sc1604/LFCfluxcalib/' + outcat
readcol, incat, format='A,A,A, D, D, D', LFCID, mask, slit, rband, iband, zband 
ncol = n_elements(LFCID)

readcol, filter + '.dat', format = 'D,D,D', lam, pointthrough, extendedthrough, /SILENT
minlam = min(lam)       ;min wavelength of the filter transmission curve
maxlam = max(lam)	
openw, lun, outcat, /get_lun 
	
	for i=0, ncol-1 do begin

	spec = d2 + 'sc1604/' + mask[i] + '/*/spec1d.' + mask[i] + '.' + slit[i] + '.' + LFCID[i] + '.fits'
	;h = headfits(spec, ext=1)
        ;exptime = sxpar(h,'EXPTIME')
	
	ss1d=fill_gap(spec,/tweak,/horne,/telluric,/silent,header=header) 	; do the fill_gap stuff, remove A/B band and interpolate over the CCD gap to make one spectrum that contains both the red and blue CCDs with Horne extraction, changed to boxcar extraction (rather than Horne) to minimize aperture effects, response correction is also done here, check fix_response.pro to make sure it is coming out in proper form 
       	wave=ss1d.lambda
       	airtovac,wave	; tweak the wavlenegths for the index of refraction of air
	ss1d.lambda=wave
	
	spec_in=ss1d.spec
     	lambda_in=ss1d.lambda
     	ivar_in=ss1d.ivar

	print, spec	
	constfluxcorr = 2e-08/(3600.*!pi*(996./2)^2) 	;2*10^-8 erg*Ang/e-, 3600 s per frame (units are in counts per hour) , pi*498 cm ^2 effective Keck II aperture
	corrfactor = deimos_correction(ss1d.lambda)	;deimos response correction as a function of wavelength, multiplicative factor
	corr_spec_in = fltarr(n_elements(ss1d.spec))
	
		for j=0, n_elements(corrfactor)-1 do begin
		corr_spec_in[j] = spec_in[j]*constfluxcorr/(lambda_in[j])	;get a flux in units of ergs/s/pixel/cm^2 see http://www.ucolick.org/~npk/dokuwiki/doku.php?id=deep2_flux 
		endfor

	dummy = where(lambda_in ge minlam and lambda_in le maxlam)      ;check where the actual spectrum overlaps the filter transmission curve
        arraylength = n_elements(dummy)
        maxlambda = lambda_in[dummy[arraylength-1]]
        minlambda = lambda_in[dummy[0]]
	througharray = generate_interpol_filtercurve(filter, arraylength, minlambda, maxlambda)	;interpolate values so that the interpolated grid matches up with the spectrum grid, both in intial wavelength values and spacing

	deltafband = fltarr(n_elements(dummy))
	meanf = fltarr(n_elements(dummy))
                for j=0, n_elements(dummy)-2 do begin
                deltafband[j] = abs(3.0D+8/(lambda_in[dummy[j]]*1.0D-10)-3.0D+8/(lambda_in[dummy[j+1]]*1.0D-10))   ;the frequency interval that the filter throughput covers, needed to get a flux density for converting to AB magnitudes 
                meanf[j] = (3.0D+8/(lambda_in[dummy[j]]*1.0D-10) + 3.0D+8/(lambda_in[dummy[j+1]]*1.0D-10))/2	;get the effective frequency in each bin
		endfor


	dummy2 = where(deltafband eq 0)		;clean up last element
	deltafband[dummy2] = mean(deltafband)
	dummy3 = where(meanf eq 0)
	meanf[dummy3] = max(meanf)
	bandfluxnum = fltarr(n_elements(dummy))
	bandfluxdenom = fltarr(n_elements(dummy))

	        for j=0, n_elements(dummy)-1 do begin
                bandfluxnum[j] = througharray[j]*corr_spec_in[dummy[j]]/deltafband[j]*deltafband[j]/meanf[j]	;this is pedantic see Fukugita 1995 eqn 9, F_nu is corr_spec_in/deltafband
		bandfluxdenom[j] = througharray[j]*deltafband[j]/meanf[j]
                endfor

	totbandfluxdens = total(bandfluxnum)/total(bandfluxdenom)

	;splot, lambda_in, spec_in	

	
	ZP = 48.60	;photometric zero point for AB system
	extinction = -0.1	;approximate extinction in magnitudes/airmass at Mauna Kea, see http://www.cfht.hawaii.edu/Instruments/ObservatoryManual/CFHT_ObservatoryManual_(Sec_2).html	
	airmass = sxpar(h,'AIRMASS')
	filterABmag = -ZP - 2.5*alog10(totbandfluxdens) + airmass*extinction
	print, iband[i], filterABmag
	printf, lun, LFCID[i], '   ', mask[i], '   ',slit[i], '   ', rband[i], '   ', iband[i], '   ', zband[i], '   ', filterABmag, format ='(a15, a3, a6, a3, a4, a3, d8.5, a3, d8.5, a3, d8.5, a3, d8.5)'
	endfor
	
free_lun, lun
end
		
