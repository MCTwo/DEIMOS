function smooth_spec_kernel, x
	
	halfwidth = floor(1.5*x)
	kernel=findgen(2*halfwidth+1)-halfwidth
	kernel=exp(-kernel^2/2/x^2)
	return, kernel/total(kernel)
end

pro fitlineprofile, spec
	
	t1 = systime(1)
	spectrum = mrdfits(spec,1)
	range = where(spectrum.lambda gt 1210 and spectrum.lambda lt 1222)	;limiting the range so things don't get fucked up by edge effects, used to be 1210 to 1222	

	specerr = 1/sqrt(spectrum.ivar[range])
	flux = spectrum.spec[range]	
	lambda = spectrum.lambda[range]

	;ok going to do the fitting myself inputting the parameters as a bunch of arrays
	meanlam = findgen(17)*0.1 + 1214
	sigma = findgen(14)*0.05 + 0.55
	amplitude = findgen(10)*25 + 500	;used to be 24 element starting at 400
	background = findgen(6)*10 - 20
	convolsigma = findgen(16)*0.1 + 2.0
	;trunclam = findgen(11)*0.02 - 0.1
	x = lambda	
	;blueside = where(x lt 1215.5)	;not doing it right at 1215.7 cuz my centering sucks
	;redside = where(x ge 1215.5) 

	chisqr = 1000000	
	
	for i=0, n_elements(meanlam)-1 do begin
	
		for j=0, n_elements(sigma)-1 do begin
		
			for k=0, n_elements(amplitude)-1 do begin
			
				for m=0, n_elements(background)-1 do begin
	
					for n=0, n_elements(convolsigma)-1 do begin	

						;for h=0, n_elements(trunclam)-1 do begin
	
						blueside = where(x lt meanlam[i]); + trunclam[h])   ;not doing it right at 1215.7 cuz my centering sucks
					        redside = where(x ge meanlam[i]); + trunclam[h])

						origgauss = background[m] + amplitude[k]*exp(-(x-meanlam[i])^2/(2*sigma[j]^2))
						gauss = origgauss	
						gauss[redside] = gauss[redside]
						gauss[blueside] = 0.
	    					gauss = convol(gauss, smooth_spec_kernel(convolsigma[n]), /center)
						chisqrtemp = total( (gauss-flux)^2/specerr^2 ) 		
						
						if chisqrtemp lt chisqr then begin
							
							bestogauss = origgauss
							bestgauss = gauss	
							chisqr = chisqrtemp
							redchisqr = chisqr/n_elements(flux)
							;besttrunc = trunclam[h]+meanlam[i]
							bestlam = meanlam[i]
							bestsig = sigma[j]
							bestamp = amplitude[k]
							bestback = background[m]
							bestconvol = convolsigma[n]
							iiter = float(i)
							jiter = float(j)
							kiter = float(k)
							miter = float(m)
							niter = float(n)	
							hiter = 0 ;float(h)
							count = ((iiter)*(jiter+1.)*(kiter+1.)*(miter+1.)*(niter+1.))+((jiter)*(kiter+1.)*(miter+1.)*(niter+1.))+ ((kiter)*(miter+1.)*(niter+1.)) + (miter)*(niter+1.) + niter + ((hiter)*(iiter+1.)*(jiter+1.)*(kiter+1.)*(miter+1.)*(niter+1.))
							print, 'At iteration number', count
							print, 'The chi-squared is', chisqr, ' and a reduced chi-squared of', redchisqr
							;print, 'Truncation wavelength=', besttrunc
							print, 'mean=', bestlam
							print, 'sigma=', bestsig
							print, 'A =', bestamp
							print, 'C =', bestback
							print, 'smooth sigma=', bestconvol
						endif	
					;endfor
				endfor
			endfor
		endfor
	endfor
	endfor
loadcolors

print, systime()+'  Done with minimization!'
print
print, 'Total time elapsed: ', (systime(1)-t1)/60., ' minutes.'

sigv = (bestsig/bestlam)*3e+05

openw, lun, 'original_gauss_bestfit.dat', /get_lun
	
	printf, lun, bestsig, '  ', bestlam, '  ', bestamp, '   ', bestback, '   ', bestconvol
	for i=0, n_elements(bestogauss)-1 do begin
	printf, lun, x[i], '  ', bestogauss[i]
	endfor 
free_lun, lun

openw, lun, 'truncated_gauss_bestfit.dat'

printf, lun, bestsig, '  ', bestlam, '  ', bestamp, '   ', bestback, '   ', bestconvol
        for i=0, n_elements(bestgauss)-1 do begin
        printf, lun, x[i], '  ', bestgauss[i]
        endfor 
free_lun, lun



;plotting crap

;splot, lambda, flux, psym=7, yrange=[-50, 800]
;oplot, x, bestogauss, color=5, linestyle=2
;oplot, x, bestgauss, color=6, linestyle=0

end

;gauss = p[0] + p[1]*exp(-(x-p[2])^2/(2*p[3]^2))
	
;	if x lt 1215.7 then gauss = 0 else gauss = gauss	
;	gauss = smooth_spec_forfit, gauss, lambda, err, p[4]		

;	halfwidth = floor(1.5*p[4])
 ;       kernel=findgen(2*halfwidth+1)-halfwidth
  ;      kernel=exp(-kernel^2/2/p[4]^2)
   ;     kernel=kernel/total(kernel)
    ;    gauss =  convol(gauss,  kernel,  /center)

;	start = [0.D, 500., 1214., 4.];, 5.]  
;	result = mpfitfun('gaussconvol', lambda, flux, err, start)

;	print, result
;	splot, lambda, flux, yrange=[-200,500], psym=3
	;oplot, lambda, smooth_spec_forfit, result[0] + result[1]*exp(-(x-result[2])^2/(2*result[3]^2)), lambda, err, result[4]
	
;end
