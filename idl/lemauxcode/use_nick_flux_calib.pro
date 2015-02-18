pro use_nick_flux_calib, spec, mag, z

	s = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/'
	spec = s + spec
	h = headfits(spec, ext=1)
        exptime = sxpar(h,'EXPTIME')
	spec = fill_gap(spec, /tweak,/horne,/telluric,/silent,header=header)
	fltr = mrdfits('i.fits',1)	

	result = flux_calib_spec(spec, fltr, mag, z)
	
	openw, lun, 'SC1NM1.041.cat', /get_lun
	
	for i=0, n_elements(result.spec)-1 do begin
	printf, lun, result.lambda[i], '    ', result.spec[i]
	endfor
	
	free_lun, lun

	spec.spec = 2.16e-17/(exptime/3600)*spec.spec/spec.lambda

	openw, lun, 'SC1NM1.041.fake.cat', /get_lun
	
	for i=0, n_elements(spec.spec)-1 do begin
	printf, lun, spec.lambda[i], '   ', spec.spec[i]
	endfor

	free_lun, lun

	print, mean(result.spec), mean(spec.spec)
end
	
	
