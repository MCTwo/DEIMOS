function generate_interpol_filtercurve_forKcorr_z, filter, z, arraylength, minlambda, maxlambda

	result = -1
	readcol, filter + '.dat', format = 'D,D', lam, through, /SILENT
	lam = lam*(1+z)
	mocklam = minlambda + findgen(arraylength)*(maxlambda-minlambda)/arraylength
	result = interpol(through, lam, mocklam, /quadratic)		;quadratic significantly better than spline and lsquadratic by eye, slightly better than linear	

	;openw, lun, outcat, width=100, /get_lun
	
	;for i=0, n_elements(results)-1 do begin
		;print, result[i]
		;mockthrough[i] = result[i]*mocklam^[i] + mockthrough[i]
	;endfor
	;free_lun, lun
	
	;splot, mocklam, mockthrough, color=2, psym=2
	;oplot, lam, avgthrough, psym=4
	;return, mocklam, mockthrough
	return, result
end
		
	
