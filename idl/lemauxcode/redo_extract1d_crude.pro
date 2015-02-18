pro redo_extract1d_crude, blueslit, redslit, ypos, fwhm 
	
	mm = mrdfits(blueslit,1)
	ss = mrdfits(redslit,1)
	a = fltarr(4096)
	b = fltarr(4096)
	c = fltarr(4096)
	d = fltarr(4096)
	s1d = {lambda:fltarr(8192), spec:fltarr(8192), ivar: fltarr(8192)} 
	f = 1/mm.ivar
	g = 1/ss.ivar
	for i=0, 4095-1 do begin
	a[i] = total(mm.flux[i, ypos-fwhm/2:ypos+fwhm/2])
	endfor
	
	for i=0, 4095-1 do begin
        b[i] = total(ss.flux[i, ypos-fwhm/2:ypos+fwhm/2])
        endfor
	
	for i=0, 4095-1 do begin
        c[i] = 1/total(f[i, ypos-fwhm/2:ypos+fwhm/2])
        endfor

	for i=0, 4095-1 do begin
        d[i] = 1/total(g[i, ypos-fwhm/2:ypos+fwhm/2])
        endfor


	s1d.lambda[0:4095] = mm.lambda0
	s1d.lambda[4096:8191] = ss.lambda0
	s1d.spec[0:4095] = a
	s1d.spec[4096:8191] = b
	s1d.ivar[0:4095] = c
	s1d.ivar[4096:8191] = d
	mwrfits, s1d, 'test.fits'
	
end
