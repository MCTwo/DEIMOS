pro combine_NIR_spec, name

	name2 = 'spec1d.NIR.' + name  
	m = mrdfits(name2+'_spec.fits')
	n = mrdfits(name2+'_lam.fits')
	o = mrdfits(name2+'_err.fits')
		
	oivar = 1/o^2
	
	result = {lambda:n, spec:m, ivar:oivar}
	mwrfits, result, 'spec1d.NIR.' + name + '.fits'
end
