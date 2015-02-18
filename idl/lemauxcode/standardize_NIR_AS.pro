; notes:
; - I assume that the spectrum is normalized to counts per second, all NIRSPEC standards we took were 4 coadds of 0.25 so good for those
; - Magnitudes should be in J in the 2MASS system since im scaling a Vega spectrum based on Vega's 2MASS J magnitude
; - Standard should be an A0 star, since I am scaling Vega's spectrum

pro standardize_NIR, spec, magnitude, airmass


	m = mrdfits(spec,1)
	
	newspec = dblarr(n_elements(m.spec))
	newivar =  dblarr(n_elements(m.ivar)) 
	deltalam = deriv(m.lambda)
	
	extinction = 0.114	;+/- 0.007, extinction per airmass on Mauna Kea in J band

	print, n_elements(m.spec), n_elements(newspec), n_elements(m.lambda), n_elements(deltalam)
	
	for i=0, n_elements(m.spec)-1 do begin

		newspec[i] = 10^(airmass*extinction/2.5)*m.spec[i]	;h*c/(A_eff*lambda*delta_lambda), is a flux density
		newivar[i] = (1/(10^(airmass*extinction/2.5)*1/sqrt(m.ivar[i])*6.63e-27*3e18/(6.33e05*m.lambda[i]*deltalam[i])))^2	;same 

	endfor

	delta_m = magnitude 		;magnitude difference between the standard and vega in J (m_J_vega = -0.177 +/- 0.206 from 2MASS)
	readcol, 'vega.005.fnu', format='D, D', vega_temp_lam, vega_fnu
	vega_jansks = vega_fnu*1e23
	vega_lam = vega_temp_lam/10000.

	roi = where(vega_lam ge m.lambda[0] and vega_lam le m.lambda[n_elements(m.lambda)-1])
	
	std_fluxdens = 10^(-delta_m/2.5)*vega_jansks[roi]		; I should put in errors here from the error in Vega's J band measure
	
	interpol_spec = interpol(newspec,m.lambda,vega_lam[roi], /quadratic)
	
	loadcolors

	throughput = interpol_spec/std_fluxdens
	splot, vega_lam[roi], throughput

	result = {lambda:vega_lam[roi], spec:throughput}
	mwrfits, result, 'HD203856.countspersperJy.fits'
end	
