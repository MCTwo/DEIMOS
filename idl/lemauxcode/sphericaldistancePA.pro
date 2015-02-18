pro sphericaldistancePA, cat

	readcol, cat, format='D,D,A,A,A', RA, dec, mask, slit, priority

	RArad = RA*2.*!PI/360.
	decrad = dec*2.*!PI/360.
	
	N = n_elements(RA)-1
	k = 0.	
	totalsize = (N^2+N)/2

	d = dblarr((N^2+N)/2)
	galaxy1 = strarr((N^2+N)/2)
	galaxy2 = strarr((N^2+N)/2)
	p = dblarr((N^2+N)/2)
	priority1 = strarr((N^2+N)/2)
	priority2 = strarr((N^2+N)/2)	

	for i=0, N-1 do begin
	
	j=1.
		while j le N-i do begin
		p[k] = 360/(2*!pi)*atan( (-cos(decrad[i+j])*sin(RArad[i]-RArad[i+j]))/(cos(decrad[i])*sin(decrad[i+j])-sin(decrad[i])*cos(decrad[i+j])*cos(RArad[i] - RArad[i+j]))) 
		if p[k] lt 0 and RArad[i] gt RArad[i+j] then p[k] = 360. + p[k] else if p[k] lt 0 and RArad[i] lt RArad[i+j] then p[k] = 180. + p[k] else if p[k] gt 0 and RArad[i] gt RArad[i+j] then p[k] = 180. + p[k] else if p[k] gt 0 and RArad[i] lt RArad[i+j] then p[k] = p[k]
		d[k] = 360/(2*!pi)*3600*acos(sin(decrad[i])*sin(decrad[i+j])+cos(decrad[i])*cos(decrad[i+j])*cos(RArad[i]-RArad[i+j]))
		galaxy1[k] = mask[i]+'.'+slit[i]
		galaxy2[k] = mask[i+j]+'.'+slit[i+j]
		priority1[k] = priority[i]
		priority2[k] = priority[i+j]
		j = j+1.
		k = k+1.
		endwhile 
	endfor

	a = n_elements(d)
	print, a
	openw, lun, 'output.dat', /get_lun
	for i=0, totalsize-1 do begin
	printf, lun, galaxy1[i], galaxy2[i], priority1[i], priority2[i], d[i], p[i], format='(a10, 3x, a10, 3x, a3, 3x, a3, 3x, d12.6, 3x, d10.6)'		
	endfor
	free_lun, lun
end

