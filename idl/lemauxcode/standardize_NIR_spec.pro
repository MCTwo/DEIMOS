pro standardize_NIR_spec, cat, standard

	readcol, cat, format='A,A,A,D', mask, slit, spectrum, z, airmass
	
	m = mrdfits(standard,1)
	stdflux = m.spec
	stdlam = m.lambda

	extinction = 0.114 	;+/- 0.007, extinction per airmass on Mauna Kea in J band
	

	for i=0, n_elements(spectrum)-1 do begin
	
		NIR_spec = mrdfits(spectrum[i],1)
		NIR_flux = NIR_spec.spec*10^(airmass[i]*extinction/2.5)/900.		;convert to per second (15 minute exposures) and airmass correct
		NIR_lam = NIR_spec.lambda*10^4		;convert from microns to Angstroms
		NIR_err = sqrt(1/NIR_spec.ivar)/900.	;should use the error in the standard here as well, only using the spectral uncertainties
		rest_NIR_lam = NIR_lam;/(1.+z[i])

		shifted_stdflux = interpol(stdflux,stdlam,NIR_lam, /quadratic)
		print, mask[i], '   ', slit[i]
		std_NIR_spectrum = NIR_flux/shifted_stdflux
		std_NIR_ivar = 1/(NIR_err/shifted_stdflux)^2
		result = {lambda:rest_NIR_lam, spec:std_NIR_spectrum, ivar:std_NIR_ivar}
		mwrfits, result, 'standardized.'+spectrum[i]
		;mwrfits, result, 'standardized.restframe.'+spectrum[i]
	endfor

	
end

	
