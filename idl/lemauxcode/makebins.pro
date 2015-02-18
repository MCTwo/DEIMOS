pro makebins, deltax, dlogz, numx
	
	xbegin = fltarr(numx)
	xend = fltarr(numx)
	
	for i=0, numx-1 do begin
	xbegin[i] = 3282. - deltax*numx/2 + deltax*i
	xend[i] = 3282. - deltax*(numx-2)/2 + (deltax*i-1.)
	endfor
	
	;for i=0, numy-1 do begin
	;ybegin[i] = 1573. - (deltay*numy/2 - deltay*i)
	;yend[i] = 1573. - (deltay*(numy-2.)/2 - (deltay*i-1.))
	;endfor

	if n_elements(dlogz) eq 0 then dlogz = 4e-3
	zposlimit = [1573, 2323]
	lzposlimit = alog10(zposlimit)
	npoints = long((lzposlimit[1]-lzposlimit[0])/dlogz)
	z1 = findgen(npoints)*dlogz + lzposlimit[0]
	zinput = 10.^z1
	print, zinput
	
	zposbegin = fltarr(npoints)
        zposend = fltarr(npoints)
        znegbegin = fltarr(npoints)
        znegend = fltarr(npoints) 
	
	for i=0, npoints-3 do begin
	zposbegin[i] = zinput[i]
	zposend[i] = (zinput[i+1])-1 
	znegbegin[i] = 1573 - ((zinput[i]) - 1573)
	znegend[i] = 1573 - (zinput[i+1] - 1573)+1
	endfor
	
	

	openw, lun, 'output.dat', /get_lun, width = 400
   	
	for i=0, numx-1 do begin
		for j=0, npoints-1 do begin
   		printf, lun, xbegin[i], xend[i], zposbegin[j], zposend[j], format='(i4, 2x, i4, 2x, i4, 2x, i4)'
		printf, lun, xbegin[i], xend[i], znegend[j], znegbegin[j], format='(i4, 2x, i4, 2x, i4, 2x, i4)'
   		endfor
	endfor
   	free_lun, lun
end
	
