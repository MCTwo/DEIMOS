pro spec1dtoascii_1604, ID, mask, slit, horne=horne, boxcar=boxcar
		
	spec = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/' + mask + '/*/spec1d.' + mask + '.' + slit + '.' + ID + '.fits'
	if n_elements(boxcar) then boxcar = boxcar[0] ge 1 else boxcar = 0
	if n_elements(horne) then horne = horne[0] ge 1 else horne = 0
	if boxcar eq 0 then horne = 1

	if horne eq 1 then ss1d = fill_gap(spec,/tweak, /telluric,/silent,header=header)
	
	if boxcar eq 1 then ss1d = fill_gap(spec,/tweak, /boxsprof, /telluric,/silent,header=header)
	
	wave=ss1d.lambda
        airtovac,wave   ; tweak the wavlenegths for the index of refraction of air
        ss1d.lambda=wave
	
	spec_in=ss1d.spec
        lambda_in=ss1d.lambda
        ivar_in=ss1d.ivar
	; could search for redshift in the catalog by ID and make a restfram grid option	
	result = {lambda:lambda_in, spec:spec_in, ivar:ivar_in}
	mwrfits, result, 'spec1d.' + mask + '.' + slit + '.' + ID + '.fillgapped.fits'
		
	openw, lun, mask + '.' + slit + '.' + ID + '.new.dat', /get_lun

		for i=0, n_elements(spec_in)-1 do begin
		printf, lun, lambda_in[i], '   ', spec_in[i], format='(d10.5, 3a, d9.5)'
		endfor
	
	free_lun, lun

end	
		
