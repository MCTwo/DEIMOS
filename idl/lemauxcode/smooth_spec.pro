pro smooth_spec, spec, sigma, sigstring
	
	spec2 = mrdfits(spec,1)

	halfwidth = floor(1.5*sigma)
        kernel=findgen(2*halfwidth+1)-halfwidth
        kernel=exp(-kernel^2/2/sigma^2)
        kernel=kernel/total(kernel)
        stemp_flux =  convol(spec2.spec,  kernel,  /center)
	stemp_ivar = convol(spec2.ivar, kernel, /center)
	
	result = {lambda:spec2.lambda, spec: stemp_flux, ivar: stemp_ivar}
	
	sigstring = strcompress(string(sigma), /REMOVE_ALL)
	mwrfits, result, 'smoothed.10sigma.' + spec
	splot, result.lambda, result.spec, yrange=[-200,600]	
end
