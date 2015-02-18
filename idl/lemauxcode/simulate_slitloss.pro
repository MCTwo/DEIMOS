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

pro simulate_slitloss, rhlow, rhhigh, seeing, outcat, deV=deV, exponential=exponential
	
	if n_elements(deV) then deV = deV[0] else deV = 0 
	if n_elements(exponential) then exponential = exponential[0] else exponential=0
	if exponential eq 0 then deV = 1
	pixelscale = 0.03	;arcseconds/pixel
	x = findgen(1000)*pixelscale
	y = findgen(1000)*pixelscale
	seeingFWHM = seeing/pixelscale ;want the FWHM not sigma ;/2.3548/pixelscale	;standard deviation from FWHM in pixels
	r= fltarr(n_elements(x), n_elements(y))
	slitbound = fltarr(n_elements(x), n_elements(y))
	
	for i=0, n_elements(x)-1 do begin
		for j=0, n_elements(y)-1 do begin
		r[i,j] = sqrt((x[i]-500*pixelscale)^2+(y[j]-500*pixelscale)^2)
		if abs(x[i]-500*pixelscale) le 0.5 then slitbound[i,j] = 1 else slitbound[i,j] = 0	;between -0.5" and 0.5" along direction perpendicular to the slit (x direction), so that the slit area extends, realistically, many " in the y direction
		endfor
	endfor


	;modx = abs(500*pixelscale-x)
	;mody = abs(500*pixelscale-y)
	
	;seeingxarray =  1/(seeingsig*sqrt(2*!PI))*exp(-(seex-499*pixelscale)^2/(2*seeingsig))
	;seeingyarray =  1/(seeingsig*sqrt(2*!PI))*exp(-(seey-499*pixelscale)^2/(2*seeingsig))
	Ic = 1*10^4 	;central intensity of the galaxy, taking ratios in the end so kinda irrelevant	
	
	rh = fltarr(12)
	for i=0, 12-1 do begin
		rh[i] = rhlow + (rhhigh-rhlow)/12*i	;make 12 bins of half-light radii with values equally spaced between the low and high values
	endfor

	;OLD, doing convolution in 2d now	
	;rhx = rharray/sqrt(2)		;can only do convolution in 1d so have to project r assuming the profile is symmetric in x and y
	;rhy = rharray/sqrt(2)

	openw, lun, outcat, /get_lun
	for i=0, 12-1 do begin 

		if deV eq 1 then begin

		;fluxarr = fltarr(n_elements(modx), n_elements(mody))
		;posarr = fltarr(n_elements(modx), n_elements(mody))
		;seeingarr = fltarr(n_elements(modx), n_elements(mody))
		;convolarr = fltarr(n_elements(modx), n_elements(mody))
	
                        ;for j=0, n_elements(modx)-1 do begin
                                ;for k=0, n_elements(mody)-1 do begin
                fluxarr = Ic*exp(-7.67*((r/rh[i])^(0.25)-1))
                                ;seeingarr[j,k] = 1/(seeingsig^2*2*!PI)*exp(-modx[j]^2/(2*seeingsig^2) - mody[k]^2/(2*seeingsig^2))
				;posarr = sqrt(modx[j]^2+mody[k]^2)
				;endfor
                        ;endfor
		convolarr = filter_image(fluxarr, smooth=(seeingFWHM/2)^2, /ITER)	;this is equivalent to convolving with a Gaussian of FWHM = seeingFWHM, see the notes in filter_image 
		;gfilter, fluxarr, seeingsig*2, [seeingsig/pixelscale, seeingsig/pixelscale], convolarr	;convolve a Gaussian array of dimensions 1.5 seeing sigma x 1.5 seeing sigma with sigma defined by the seeing input
		;xflux = 1*10^4*exp(-7.67*((modx/rhx[i])^(0.25)-1))
        	;yflux = 1*10^4*exp(-7.67*((mody/rhy[i])^(0.25)-1))
		endif

		if exponential eq 1 then begin
		;xflux = 1*10^4*exp(xmod/rhx[i])
                ;yflux = 1*10^4*exp(ymod/rhy[i])
		
		;fluxarr = fltarr(n_elements(xmod), n_elements(ymod))
		;posarr = fltarr(n_elements(modx), n_elements(mody))

		fluxarr = Ic*exp(-r/rh[i])
                        ;for j=0, n_elements(xflux)-1 do begin
                                ;for k=0, n_elements(yflux)-1 do begin
                                ;fluxarr[j,k] = Ic*exp((modx[j]^2+mody[j]^2)/rharray[i])
 				;seeingarr[j,k] = 1/(seeingsig^2*2*!PI)*exp(-modx[j]^2/(2*seeingsig^2) - mody[k]^2/(2*seeingsig^2))
                                ;convolarr[j,k] = fluxarr[j,k]*seeingarr[j,k]
	                        ;posarr[j,k] = sqrt(modx[j]^2+mody[k]^2)
				;endfor
			;endfor
		;gfilter, fluxarr, 11, seeingsig, convolarr
		convolarr = filter_image(fluxarr, smooth=(seeingFWHM/2)^2, /ITER)     ;this is equivalent to convolving with a Gaussian of FWHM = seeingFWHM, see the notes in filter_image
		endif
		
	
		;convolx = blk_con(seeingxarray, fluxarr, B_length=201)
		;convoly = blk_con(seeingyarray, fluxarr, B_length=201)
		
		;convolfluxarr = fltarr(n_elements(convolx), n_elements(convoly))
	
		;	for j=0, n_elements(convolx)-1 do begin
		;		for k=0, n_elements(convoly)-1 do begin
		;		convolfluxarr[j,k] = convolx[j] + convoly[k]
		;		endfor
		;	endfor	
		;window, /free, xsize=1000, ysize=1000
		;tv, fluxarr
		;slitbound = where(r le 1)     ;slit is 1", 0.5 in either direction from the central point of the galaxy, 
		slitsize = dblarr(2)    ;2 array with slit width and slit length 1" x 42"
                slitsize[0] = 1./0.18
                slitsize[1] = 42./0.18
               
                ;tvbox, slitsize, 499., 499., color='0000FF'XL
		
		;window, /free, xsize=1000, ysize=1000
		;device, decomposed=1
		;tv, convolarr
		;slitbound = where(1. - r gt 0)	;slit is 1", 0.5 in either direction from the central point of the galaxy, by flooring maximizing slit loss
		slitsize = dblarr(2)    ;2 array with slit width and slit length 1" x 42"
        	slitsize[0] = 1./0.18
        	slitsize[1] = 42./0.18  

        	;tvbox, slitsize, 499., 499., color='0000FF'XL
		
			
		;slitboundy = [ceil(500-rhy[i]/(2*pixelscale)), ceil(500+rhy[i]/(2*pixelscale))] 	; extraction window is typically no more than 2.2", on average is ~ 14*.118 
		;rhbound = where(rh[i] - r gt 0)
		;halflight = total(fluxarr[rhbound])
		slitconvolfluxarr = convolarr*slitbound
		slitfluxarr = fluxarr*slitbound;slitboundx[0]:slitboundx[1],slitboundy[0]:slitboundy[1]]
		slitconvolflux = total(slitconvolfluxarr)
		totflux = total(fluxarr)
		totconvolflux = total(convolarr)
		slitflux = total(slitfluxarr)

		a = fluxarr-convolarr
		;atv, a
			
		totslitloss = 1-slitconvolflux/totflux
		seeingloss = 1 - slitconvolflux/slitflux
		;fluxrat = halflight/totflux
		;print, 'This should be .50,', fluxrat 
		print, 'The total slit loss for an effective radius of', rh[i], 'is', totslitloss, 'and the loss due to seeing is', seeingloss
		;print, 'These two should be the same. The total flux in the galaxy is:',totflux, 'and convolved with the seeing the total flux is:', totconvolflux 
		printf, lun, rh[i], '   ', totslitloss, format = '(d7.4, a3, d7.4)'
	endfor
	free_lun, lun

end		
			
