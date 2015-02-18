pro maketargetlist, mask
	path = getenv('D2_RESULTS')
	slitfile = path + 'sc1604/' + mask + '/*/' + mask + '.bintabs.fits'
	s2nfile = path + 'sc1604/' + mask + '/*/obj_info.' + mask + '.fits'
	a = mrdfits(slitfile,1)
	b = mrdfits(slitfile,3)
	c = mrdfits(s2nfile,1)
	objname = strcompress(a.object, /remove_all)
	RA = a.RA_OBJ
	dec = a.DEC_OBJ
	mag = a.mag
	class = a.objclass
	slit = strcompress(b.slitname, /remove_all)
	slitlength = b.slitlen
	slitPA = b.slitlPA
	
	N = n_elements(objname) 
	s2nB = fltarr(N)
	s2nR = fltarr(N)
	objposB = fltarr(N)
	objposR = fltarr(N)
	slitno = fltarr(N)

	for i=0, N-1 do begin
	indexB = where(objname[i] eq c.objno and c.color eq 'B')
	if indexB ne -1 then s2nB[i] = c[indexB].s2n_fwhm else s2nB[i] = 0.	
	indexR = where(objname[i] eq c.objno and c.color eq 'R')
	if indexR ne -1 then s2nR[i] = c[indexR].s2n_fwhm else s2nR[i] = 0.
	indexopB = where(objname[i] eq c.objno and c.color eq 'B')
	if indexopB ne -1 then objposB[i] = c[indexopB].objpos else objposB[i] = 0.
	indexopR = where(objname[i] eq c.objno and c.color eq 'R')
        if indexopR ne -1 then objposR[i] = c[indexopR].objpos else objposR[i] = 0.		
	index = where(objname[i] eq c.objno and c.color eq 'R')
	if index ne -1 then slitno[i] = c[index].slitno else slitno[i] = 0.
	endfor 

	target = replicate({target}, N)
	
	for i=0, N-1 do begin  	
	target[i] = {target, objname:objname[i], slitlen:slitlength[i], PA:slitPA[i], s2nB:s2nB[i], s2nR:s2nR[i], objposB:objposB[i], objposR:objposR[i], slit:slit[i], RA:RA[i], dec:dec[i], imag:mag[i]}	
	endfor 

	mwrfits, target, 'target.'+ mask + '.fits'

end	
