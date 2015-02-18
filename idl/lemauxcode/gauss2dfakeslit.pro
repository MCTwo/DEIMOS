pro gauss2dfakeslit, num
	
	ngauss = round(randomn(seed,1, /uniform)*3)+3	;generate between 0 and 3 Gaussians with equal probability of each number 
	print, ngauss
	
	if ngauss eq 0 then begin
		slitarray = fltarr(8192,55) ;make a slit which is 2700 Angstroms by 6.5"
        	slitarray = slitarray + 100.      ;background eventually 0 counts, pretty realistic in the DEIMOS rsponse corrected case, for now set to 100 to generate noise

        for i=0, 8192-1 do begin
                for j=0, 55-1 do begin
                        slitarray[i,j] = randomu(seed, poisson=0.1>slitarray[i,j], /double)-100.  ; add poisson noise and subtract off the background
                endfor
        endfor
	endif
	
	if ngauss ge 1 then begin
		sigx = fltarr(ngauss) + 4. 	;average dispersion in spectral dimension for 1604 sngl non-LYA in observed frame, make ngauss of em
		ranuniform = randomn(seed, 4, ngauss[0], /uniform)
		
		sigy = ranuniform[0,0:ngauss-1] + 4	;size in the spatial dimension, 0.85" (PSF) +/- 2*0.118"
		;help, sigy, /struct
		meanx = ranuniform[1, 0:ngauss-1]*7800 + 100	;avoid ends of slits by 33 Angstroms
		meany = ranuniform[2, 0:ngauss-1]*43 + 6 ;avoid the edge of the slit by at least 0.71"
		height = ranuniform[3, 0:ngauss-1]*5 + 7.  ;the last number n this expression should be ~ rms of the slit array

	;print, sigx, sigy, meanx, meany, height
		slitarray = fltarr(8192,55) ;make a slit which is 2700 Angstroms by 6.5"
		slitarray = slitarray + 100.	;background eventually 0 counts, somewhat realistic in the DEIMOS rsponse corrected case, for now set to 100 to make noise. 100 is based on average rms of 2d slit files of 1604 data, with rms~10
		
		for i=0, 8192-1 do begin
			for j=0, 55-1 do begin
				slitarray[i,j] = randomu(seed, poisson=0.1>slitarray[i,j], /double)-100. 	; add poisson noise and subtract off the background
			endfor
		endfor
	
		for k=0, n_elements(meanx)-1 do begin ;if multiple peaks are desired
	
		if sigy[k] ge 42./3. then begin
			message, "Gaussian too large in spatial dimension"
		endif
			
		xdimen = float(floor(3.*sigx[k]))		;array goes out to 2 sigma
		ydimen = float(floor(3.*sigy[k]))
		ellipsearray = fltarr(xdimen, ydimen)
		gaussarray = fltarr(xdimen, ydimen) 
		x = findgen(xdimen)		;indexed arraay to create x contributions
		xellipse = ((x-xdimen/2.)/sigx[k])^2 	;force the Gaussian to peak in the middle of the array, put in mean later
		y = findgen(ydimen)
		yellipse = ((y-ydimen/2)/sigy[k])^2
		
		for i=0, xdimen-1 do begin
			for j=0, ydimen-1 do begin
				ellipsearray[i,j] = xellipse[i] + yellipse[j]
			endfor
		endfor
	
		gaussarray = mean(slitarray) + height[k]*exp(-ellipsearray/2)
		gaussarray = 100. + gaussarray	;want the rms to be from background poisson fluctuations + gaussian count poisson fluctuations, original background was 9	

		;adding noise to the Gaussian which will then add in quadrature to the noise originally in slit
		for i=0, xdimen-1 do begin
			for j=0, ydimen-1 do begin
				gaussarray[i,j] = randomu(seed, poisson=0.1>gaussarray[i,j], /double)
			endfor	
		endfor
	
		gaussarray = gaussarray - 100. 		;subtract off the background again

		; now incorperate the mean and imbed the gaussian array to the slit
		
		beginx = meanx[k]-floor(xdimen/2)
		endx = meanx[k]+floor(xdimen/2)
		beginy = meany[k]-floor(ydimen/2)
		endy = meany[k]+floor(ydimen/2)

		;make sure the Gaussian array is the same size as the array that it will be imbedded to in the slit
		if (endx-beginx) ne (xdimen-1.) then endx = endx - (endx-beginx - (xdimen-1.))
		if (endy-beginy) ne (ydimen-1.) then endy = endy - (endy-beginy - (ydimen-1.))
	
		if beginx lt 0 or endx gt 8191. then begin
			message, "Feature could not be placed on slit, change x mean or sigma", /info
		endif
		
		if beginy lt 0 or endy gt 55. then begin	
			message, "Feature could not be placed on slit, change y mean or sigma", /info
		endif
	
		slitarray[beginx:endx, beginy:endy] = gaussarray ;+ slitarray[beginx:endx, beginy:endy], don't add it to original array, would increase noise by sqrt(2), replace
		area = height[k]*2*!PI*sigx[k]*sigy[k]
	
		print, 'The number of the Gaussians in the slit are', ngauss	
		print, 'Area of Gaussian number', k, 'at y =', meany[k], ' and x = ', meanx[k], 'is', area, 'counts or', area/stddev(slitarray), ' sigma' 
	
		endfor	
	endif
	
	atv, slitarray
	mwrfits, slitarray, 'fakeslit.' + num + '.fits'

end		
