pro makebins2, deltay, dlogr, numy
	
	ybegin = fltarr(numy)
	yend = fltarr(numy)
	
	for i=0, numy-1 do begin
	ybegin[i] = 1573. - deltay*numy/2 + deltay*i
	yend[i] = 1573. - deltay*(numy-2)/2 + (deltay*i-1.)
	endfor
	
	;for i=0, numy-1 do begin
	;ybegin[i] = 1573. - (deltay*numy/2 - deltay*i)
	;yend[i] = 1573. - (deltay*(numy-2.)/2 - (deltay*i-1.))
	;endfor

	if n_elements(dlogr) eq 0 then dlogr = 4e-3
	rposlimit = [3282, 5000]
	lrposlimit = alog10(rposlimit)
	npoints = long((lrposlimit[1]-lrposlimit[0])/dlogr)
	r1 = findgen(npoints)*dlogr + lrposlimit[0]
	rinput = 10.^r1
	print, rinput
	
	rposbegin = fltarr(npoints)
        rposend = fltarr(npoints)
	
	for i=0, npoints-3 do begin
	rposbegin[i] = rinput[i]
	rposend[i] = (rinput[i+1])-1 
	endfor
	
	

	openw, lun, 'r_output.dat', /get_lun, width = 400
   	
	for i=0, numy-1 do begin
		for j=0, npoints-1 do begin
   		printf, lun, rposbegin[j], rposend[j], ybegin[i], yend[i], format='(i4, 2x, i4, 2x, i4, 2x, i4)'
   		endfor
	endfor
   	free_lun, lun
end
	
