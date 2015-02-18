pro makeNIR1d, name, outname, ypos, fwhm

	mm = mrdfits(name)
	
	spec = fltarr(1024)
	lambda = fltarr(1024)

	for i=0, 1024-1 do begin
	spec[i] = total(mm[ypos-fwhm/2:ypos+fwhm/2,i])
	lambda[i] = 1.1299+i*0.0002677325 ;very rough estimate of the wavelength solution this is bullshit for real
	endfor 

	ss = {ss, lambda:lambda, spec:spec}
	mwrfits, ss, outname+'.1d.fits'
	device, retain=2
	splot, ss.lambda, ss.spec
end
