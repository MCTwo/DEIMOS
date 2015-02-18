pro makeNIR2dfits, num, setup

	sub = mrdfits('may22s00' + num +'_sub.fits')
	var = mrdfits('may22s00' + num +'_var.fits')
	lam = mrdfits('may22s00' + num +'_wav.fits')
	sky = mrdfits('may22s00' + num +'_sky.fits')
	
	ivar = 1./var

	result = {flux:sub, ivar:ivar, lambda:lam, sky:sky}
	mwrfits, result, 'slit.NIR.' + setup + '.00' + num + '.fits'
	;$rm 'may22s00' + num +'_sub.fits'
	;$rm 'may22s00' + num +'_var.fits'
	;$rm 'may22s00' + num +'_sky.fits'
	;$rm 'may22s00' + num +'_wav.fits'

end
