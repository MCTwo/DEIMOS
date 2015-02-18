function Ez, x

   distance = 1/sqrt(0.3*(1+x)^3+0.7) ;concordance cosmology assumtion omega_M = 0.3, omega_lam = 0.7, to be consistant with DEEP2
   return, distance

end

pro density_error, cat, n
	
	readcol, cat, format='A,A,D,D,D,D,D,D,D', mask, slit, z, flux, fneg, fpos, lum, lneg, lpos
	
	flux = flux/10.^(18)
	fpos = fpos/10.^(18)
	fneg = fneg/10.^(18)
	lum = lum*10.d^42
	lpos = lpos*10.d^42
	lneg = lneg*10.d^42
	binL = lum[sort(lum)]	;sort luminosities in ascending order
	
	;openw, lun, 'luminosity_bins.dat', /get_lun, width=600
	;printf, lun, binL	
	;free_lun, lun

	openw, lun, 'simulated_LAE_luminosity_bins.dat', /get_lun, width=600
	
	for i=0, n-1 do begin
	
		random = randomn(seed, n_elements(mask))	;generate a number of Gaussian random variables to equal the number of LAE candidates
		lumsim = dblarr(n_elements(lum))		;make the luminosity and flux arrays for this iteration
		fluxsim = dblarr(n_elements(flux))
		number = fltarr(n_elements(binL))

		for j=0, n_elements(random)-1 do begin
	
			if random[j] le 0 then fluxsim[j] = flux[j]+random[j]*fneg[j]	;if random variable is negative use the negative error in the flux to generate a new flux
			if random[j] gt 0 then fluxsim[j] = flux[j]+random[j]*fpos[j]	;same for positive random variables
			if fluxsim[j] lt 1.9d-18 then fluxsim[j] = 0			;impose flux cut
			lumsim[j] = (qromb('Ez', 0.0, z[j])*1.323d28*(1.+z[j]))^2*4*!PI*fluxsim[j]
		endfor
		
		for j=0, n_elements(lumsim)-1 do begin
	
			number[j] = n_elements(where(lumsim gt binL[j]))
			if number[j] lt 0 then number[j] = 0
		endfor
	printf, lun, number
	endfor
	free_lun, lun	
end	
		

	
	
		
							
		
