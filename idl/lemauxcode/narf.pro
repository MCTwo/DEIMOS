pro narf
	mm = mrdfits('V.fits',1)
	openw, lun, 'V.dat', /get_lun
	
	for i=0, n_elements(mm.lambda)-1 do begin
	printf, lun, mm.lambda[i], '    ', mm.spec[i]
	endfor
	
	free_lun, lun
end
