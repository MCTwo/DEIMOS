pro measure_1d_significance, spec, centralpix

	;centralpix is the x values that is read out by gauss2dfakeslit, alternatively it is (emission wavlength-6500)/0.33
	mm = mrdfits(spec,1)
	centrallam = centralpix*0.33 + 6500
	start = [0, centrallam, 4*0.33, 4] ;constant, mean, sigma, area; sigma is 4*0.33 as a consequencee of the dispersion sigma input in gauss2dfakeslit
	expr = 'P[0] + GAUSS1(X, P[1:3])'
	badpix = where(mm.ivar ge 10)
	mm.ivar[badpix] = 1/mean(mm.spec)
	plotrange = where(mm.lambda ge centrallam-20)
	result = MPFITEXPR(expr,mm.lambda[plotrange[0]:plotrange[0]+120], mm.spec[plotrange[0]:plotrange[0]+120], 1/sqrt(abs(mm.ivar)), start)
	print, result
	plotrange = where(mm.lambda ge centrallam-10)
	plot, mm.lambda[plotrange[0]:plotrange[0]+61], mm.spec[plotrange[0]:plotrange[0]+61]
	oplot, mm.lambda, result(0)+gauss1(mm.lambda,result(1:3)), color=50, thick=5
	;narfgauss = -0.03756 + 2.75*exp(-(mm.lambda-7185.97)^2/(2*0.7)^2)
	;oplot, mm.lambda, narfgauss, color=50, thick=5 
	specrms = stddev(mm.spec)
	print, stddev(mm.spec)
	totarea = result[3]/specrms
	print, 'The total area of the Gaussian is', totarea, 'sigma'
end
