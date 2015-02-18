pro useinterpol, cat


	;numread = 2
	readcol, cat, format='A,A,A,F,A,F', mask, slit, obj, Iband, q, z;, skipline=2,  numline = numread
	numread = n_elements(mask)
	list = replicate({list}, numread)
	
	for i=0, numread-1 do begin
   	list[i] = {list, maskname:mask[i], slitname:slit[i], objectno:obj[i],Iband:Iband[i], zquality:q[i], z:z[i]}
   	endfor
	for i=0, numread-1 do begin
	
	spec = ('/Volumes/Data2/orelse/lemaux/deimos/1322/' + list[i].maskname + '/*/spec1d.' + list[i].maskname + '.' + list[i].slitname +'.' + list[i].objectno + '.fits')[0]
	z_in = list[i].z
	;ss1d=fill_gap(spec,/tweak,/telluric,/silent,header=header)
        ;wave=ss1d.lambda
        ;airtovac,wave
        ;ss1d.lambda=wave	
	ss1db = mrdfits(spec,3)
	ss1dr = mrdfits(spec,4)
	tempspec = fltarr(8192)
	templam = fltarr(8192)
	tempivar = fltarr(8192)
	tempspec[0:4095] = ss1db.spec
	tempspec[4096:8191] = ss1dr.spec
	templam[0:4095] = ss1db.lambda
	templam[4096:8191] = ss1dr.lambda
	tempivar[0:4095] = ss1db.ivar
	tempivar[4096:8191] = ss1dr.ivar
	ss1d = {spec: tempspec, lambda:templam, ivar:tempivar}

	spec_in=ss1d.spec
     	lambda_in=ss1d.lambda
     	ivar_in=ss1d.ivar

	rlamlimit = [3300., 5600.]
	dlambda = .195
	npoints = long((rlamlimit[1]-rlamlimit[0])/dlambda)
	lambda = findgen(npoints)*dlambda +rlamlimit[0]
	ivar = lambda
	
	restframe_spec = interpol(ss1d.spec, ss1d.lambda/(1+z_in), lambda)
	restframe_ivar = interpol(ss1d.ivar, ss1d.lambda/(1+z_in), lambda)
	splot, lambda, restframe_spec
	result = { spec: restframe_spec, lambda: lambda, ivar: restframe_ivar}
	mwrfits, result, './spec/spec1d.' + list[i].maskname + '.' + list[i].slitname + '.' + list[i].objectno + '.restframe.fits'
	device, retain=2
	;splot, result.lambda, result.spec
	endfor

end
