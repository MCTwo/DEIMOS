; NAME:
;	simulateslit_loss
;
; INPUTS:
;	rhlow = smallest half light radius that is simulated (in arcseconds)
;	rhhigh = largest half light radius that is simulated (in arcseconds)
;	seeing = Gaussian seeing FWHM (in arcseonds)
;
;  KEYWORD_PARAMETERS
;	deV = specify if a de Vaucouleurs is desired (this is the default)
;	exponential = specify if an exponential profile is desired

pro simulate_slitloss_serendip, rh, seeing, outcat, deV=deV, exponential=exponential
	
	if n_elements(deV) then deV = deV[0] else deV = 0 
	if n_elements(exponential) then exponential = exponential[0] else exponential=0
	if exponential eq 0 then deV = 1
	pixelscale = 0.03	;arcseconds/pixel
	x = findgen(400)*pixelscale
	y = findgen(400)*pixelscale
	seeingFWHM = seeing/pixelscale ;want the FWHM not sigma ;/2.3548/pixelscale	;standard deviation from FWHM in pixels
	r= fltarr(n_elements(x), n_elements(y))
	slitbound = fltarr(n_elements(x), n_elements(y))

	posarray = findgen(82)*0.05 - 2	

	openw, lun, outcat, /get_lun

	for i=0, n_elements(x)-1 do begin
		for j=0, n_elements(y)-1 do begin
		if abs(x[i]-200*pixelscale) le 0.5 then slitbound[i,j] = 1 else slitbound[i,j] = 0 	;between -0.5" and 0.5" along direction perpendicular to the slit (x direction), so that the slit area extends, realistically, many " in the y direction
		endfor
	endfor	

	for k=0, n_elements(posarray)-1 do begin	
		for i=0, n_elements(x)-1 do begin
			for j=0, n_elements(y)-1 do begin
			r[i,j] = sqrt((x[i]-(200*pixelscale-posarray[k]))^2+(y[j]-(200*pixelscale))^2)
			endfor
		endfor


		Ic = 1*10^4 	;central intensity of the galaxy, taking ratios in the end so kinda irrelevant	
	
		if deV eq 1 then begin
	        fluxarr = Ic*exp(-7.67*((r/rh)^(0.25)-1))
		convolarr = filter_image(fluxarr, smooth=(seeingFWHM/2)^2, /ITER)	;this is equivalent to convolving with a Gaussian of FWHM = seeingFWHM, see the notes in filter_image 
		endif

		if exponential eq 1 then begin
		fluxarr = Ic*exp(-r/rh)
		convolarr = filter_image(fluxarr, smooth=(seeingFWHM/2)^2, /ITER)     ;this is equivalent to convolving with a Gaussian of FWHM = seeingFWHM, see the notes in filter_image
		endif
		
		;atv, fluxarr
		slitsize = dblarr(2)    ;2 array with slit width and slit length 1" x 42"
                slitsize[0] = 1./0.18
                slitsize[1] = 42./0.18
               
                ;tvbox, slitsize, 499., 499., color='0000FF'XL
		
		;window, /free, xsize=1000, ysize=1000
		;device, decomposed=1
		atv, convolarr
		;slitbound = where(1. - r gt 0)	;slit is 1", 0.5 in either direction from the central point of the galaxy, by flooring maximizing slit loss
		slitsize = dblarr(2)    ;2 array with slit width and slit length 1" x 42"
        	slitsize[0] = 1./0.18
        	slitsize[1] = 42./0.18  

			
		;slitboundy = [ceil(500-rhy[i]/(2*pixelscale)), ceil(500+rhy[i]/(2*pixelscale))] 	; extraction window is typically no more than 2.2", on average is ~ 14*.118 
		;rhbound = where(rh[i] - r gt 0)
		;halflight = total(fluxarr[rhbound])
		slitconvolfluxarr = convolarr*slitbound
		;atv, slitconvolfluxarr
		slitfluxarr = fluxarr*slitbound;slitboundx[0]:slitboundx[1],slitboundy[0]:slitboundy[1]]
		slitconvolflux = total(slitconvolfluxarr)
		totflux = total(fluxarr)
		totconvolflux = total(convolarr)
		slitflux = total(slitfluxarr)

		a = fluxarr-convolarr
		;atv, a
			
		totslitloss = 1-slitconvolflux/totconvolflux
		;fluxrat = halflight/totflux
		;print, 'This should be .50,', fluxrat 
		print, 'The total slit loss for an object position off set', posarray[k], ' from the central slit position is ', totslitloss
		;print, 'These two should be the same. The total flux in the galaxy is:',totflux, 'and convolved with the seeing the total flux is:', totconvolflux 
		printf, lun, posarray[k], '   ', totslitloss, format = '(d7.4, a3, d7.4)'
	endfor
	free_lun, lun

end		
			
