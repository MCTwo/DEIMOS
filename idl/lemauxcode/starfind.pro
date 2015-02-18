pro starfind, cat1, cat2

	readcol, cat1, format='D,D,A,A', RA, dec, mask, slit
	
	readcol, cat2, format='A,D,D,D,D,D,A', LFCID, RAcat, deccat, rband, iband, zband

	RArad = RA*2.*!PI/360.
	decrad = dec*2.*!PI/360.
	RAcatrad = RAcat*2.*!PI/360.
	deccatrad = deccat*2.*!PI/360.

	N = n_elements(RA)-1
	Ncat = n_elements(RAcat)-1
	
		

	for i=0, N do begin
	k=0.	
	d = dblarr(Ncat)
	dgood = dblarr(2000)
	p = dblarr(2000)
	s = dblarr(2000)
	galaxy1 = strarr(2000)
	galaxy2 = strarr(2000)
	rcat = dblarr(2000)
	icat = dblarr(2000)
	zcat = dblarr(2000)
	openw, lun, './starcats/' + mask[i] + '.' + slit[i] + '.staroutput.cat', /get_lun

		for j=0, Ncat-1 do begin
		d[j] = acos(cos(RArad[i]-RAcatrad[j])*cos(decrad[i])*cos(deccatrad[j]) + sin(decrad[i])*sin(deccatrad[j]))*360./(2.*!PI)*3600
		
		if d[j] lt 110. and rband[j] le 20 then begin
		dgood[k] = d[j]
		p[k] = atan(-cos(deccatrad[j])*sin(RArad[i]-RAcatrad[j])/( cos(decrad[i])*sin(deccatrad[j])-sin(decrad[i])*cos(deccatrad[j])*cos(RArad[i]-RAcatrad[j])))*360./(2.*!pi)
		if p[k] lt 0 and RArad[i] gt RAcatrad[j] then p[k] = 360. + p[k] else if p[k] lt 0 and RArad[i] lt RAcatrad[j] then p[k] = 180. + p[k] else if p[k] gt 0 and RArad[i] gt RAcatrad[j] then p[k] = 180. + p[k] else if p[k] gt 0 and RArad[i] lt RAcatrad[j] then p[k] = p[k]
		galaxy1[k] = mask[i]+'.'+slit[i]
		galaxy2[k] = LFCID[j]
		rcat[k] = rband[j]
		icat[k] = iband[j]
		zcat[k] = zband[j]		
		printf, lun, galaxy1[k], galaxy2[k], dgood[k], p[k], rcat[k], icat[k], zcat[k], format='(a10, 3x, a14, 3x, d12.6, 3x, d10.6, 3x, d10.6, 3x, d10.6, 3x, d10.6)'
		k = k+1.
		endif
		endfor	
		
	;hmm = n_elements(d_good)	
		
	;openw, lun, mask[i] + '.' + slit[i] + '.staroutput.cat', /get_lun
	;	for j=0, hmm-1. do begin
	;	if d[j] ne 0 then begin
	;	endif
	;	endfor
	free_lun, lun
	endfor
end 
