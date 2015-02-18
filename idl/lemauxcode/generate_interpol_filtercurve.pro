function generate_interpol_filtercurve, filter, arraylength, minlambda, maxlambda

	result = -1
	readcol, filter + '.dat', format = 'D,D', lam, through
	;avgthrough = (pointthrough + extendedthrough)/2		; averaging the extended source and point source transmissionfor lack of a better idea
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
		
	
