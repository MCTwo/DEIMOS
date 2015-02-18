pro useinterpol, cat

	numread = 222
	readcol, cat, format='A,A,F,F,A,A', mask, slit, Iband, z, q, file, skipline=0,  numline = numread
	list = replicate({list}, numread)
	
	for i=0, numread-1 do begin
   	list[i] = {maskname:mask[i], slitname:slit[i], Iband:Iband[i], z:z[i], zquality:q[i], file:file[i]}
   	endfor
	
	for i=0, numread-1 do begin
	
	spec = ('/Volumes/Data2/orelse/lemaux/deimos/sc1604/' + list[i].file)[0]
	z_in = list[i].z
		
	ss1d=fill_gap(spec,/tweak,/telluric,/silent,header=header)
        wave=ss1d.lambda
        airtovac,wave
        ss1d.lambda=wave	

	spec_in=ss1d.spec
     	lambda_in=ss1d.lambda
     	ivar_in=ss1d.ivar
	
	rlamlimit = [3200., 4900.]
	llimit = alog10(rlamlimit)
	dloglam = 1.95e-5
	npoints = long((llimit[1]-llimit[0])/dloglam)
	llambda = findgen(npoints)*dloglam +llimit[0]
	lambda = 10.^llambda
	ivar = lambda
	
	restframe_spec = interpol(ss1d.spec, ss1d.lambda/(1+z_in), lambda)
	restframe_ivar = interpol(ss1d.ivar, ss1d.lambda/(1+z_in), lambda)
	result = { spec: restframe_spec, lambda: lambda, ivar: restframe_ivar}
	mwrfits, result, 'spec1d.' + list[i].maskname + '.' + list[i].slitname + '.restframe.fits'
	device, retain=2
	splot, result.lambda, result.spec
	endfor

end
