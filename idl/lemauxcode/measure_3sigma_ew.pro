
;function biweight, xi
;	
;	M = 10^double(median(alog10(xi)))
;	MAD = median(abs(xi-M))
;	ui = (xi-M)/(MAD*6d)
;
;	r = where(abs(ui) lt 1d, count)
;	if count lt 1 then return, M
;	Cbi = M + total((xi[r]-M)*(1-ui[r]^2)^2)/total(1-ui[r]^2)^2
;	return, cbi
;end

pro loadcolors, bottom=bottom, names=names

        if (n_elements(bottom) eq 0) then bottom=0

        ;Load Graphics Colors
        red = [ 0, 255, 0, 255, 0, 255, 0, 255, $
                0, 255, 255, 112, 219, 127, 0, 255]
        grn = [ 0, 0, 255, 255, 255, 0, 0, 255, $
                0, 187, 127, 219, 112, 127, 163, 171]
        blu = [ 0, 255, 255, 0, 0, 0, 255, 255, $
                115, 0, 127, 147, 219, 127, 255, 127]
        tvlct, red, grn, blu, bottom

        ; Set color names

        names = ['Black', 'Magenta', 'Cyan', 'Yellow', $
                 'Green', 'Red', 'Blue', 'White', $
                 'Navy', 'Gold', 'Pink', 'Aquamarine', $
                 'Orchid', 'Gray', 'Sky', 'Beige' ]
end

; EW = feature flux/continuum flux x 
; 	Delta lambda flux/delta lambda continuum x 
;	delta lambda flux
pro measure_3sigma_ew, $
	spec, $		; The spectrum structure
	continuum, $	; The wavelength range for the continuum [[l0,l1],[l2,l3]]
	feature,$		; The wavelength range for the feature
	quiet=quiet,$
	deimos=deimos

	;add in stuff here stitching together spec2d created 1d spectra with fill_gap n airtovac if deimos keyword set

	fail = [!values.F_NAN, !values.F_NAN]
	roiF = where(spec.lambda ge feature[0] and spec.lambda le feature[1], count)
	
	if keyword_set(quiet) then q=1 else q=0
	if count le 2 then begin
		if q eq 0 then $
		message, "Wavelength range of feature is non existent", /info
		;return, fail
	endif
	dlf = deriv(spec.lambda[roiF])
	
	roi = where(spec.lambda ge continuum[0,0] and spec.lambda le continuum[1,1])
	
	roiB = where(spec.lambda ge continuum[0,0] and spec.lambda le continuum[1,0], count1) 
	roiR = where(spec.lambda ge continuum[0,1] and spec.lambda le continuum[1,1], count2)

	if count1 le 2 or count2 le 2 then begin	
		if q eq 0 then $
		message, "Wavelength range of continuum is non existent", /info
		;return, fail
	endif

	; apparently ivar is calculated from corrected counts because it's not returning correct values if I deimos_correct 1/sqrt(ivar)

	ivarB = spec.ivar[roiB]
	ivarR = spec.ivar[roiR]
	ivarF = spec.ivar[roiF]
	roivar = spec.ivar[roi]	

	bluelam = spec.lambda[roiB]
	redlam = spec.lambda[roiR]
	featurelam = spec.lambda[roiF]
	roilam = spec.lambda[roi]

	bluespec = spec.spec[roiB]
	redspec = spec.spec[roiR]
	featurespec = spec.spec[roiF]	
	roispec = spec.spec[roi]	

	Nx = n_elements(featurespec)
	Ny = n_elements(roi)
	Nb = n_elements(bluespec)
	Nr = n_elements(redspec)

	if n_elements(deimos) then begin
		dcorrF = fltarr(Nx)
		dcorrB = fltarr(Nb)
		dcorrR = fltarr(Nr)
		dcorrROI = fltarr(Ny)

		for x=0, Nx-1 do begin
        	dcorrF[x] = deimos_correction(featurelam[x])*featurespec[x]
        	endfor

		for x=0, Nb-1 do begin
		dcorrB[x] = deimos_correction(bluelam[x])*bluespec[x]
		endfor
	
		for x=0, Nr-1 do begin
		dcorrR[x] = deimos_correction(redlam[x])*redspec[x]
		endfor

		for x=0, Ny-1 do begin
		dcorrROI[x] = deimos_correction(roilam[x])*roispec[x]
		endfor	

	endif else begin

	dcorrF = featurespec
	dcorrROI = roispec
	dcorrB = bluespec
	dcorrR = redspec
	sigma = (mean(sqrt(1/ivarB)) + mean(sqrt(1/ivarR)) + mean(sqrt(1/ivarF)))/3

	endelse
	
	S = total(ivarB) + total(ivarR)
        Sx = total(bluelam*ivarB) + total(redlam*ivarR)
        Sxx = total(bluelam^2*ivarB) + total(redlam^2*ivarR)
        Sy = total(dcorrB*ivarB) + total(dcorrR*ivarR)
        Sxy = total(dcorrB*bluelam*ivarB) + total(dcorrR*redlam*ivarR)
        delta = Sxx*S - Sx^2

	y0 = (Sxx*Sy - Sx*Sxy)/delta
        m = (Sxy*S - Sx*Sy)/delta
        vary0 = Sxx/delta/(n_elements(bluelam)+n_elements(redlam))
	varm = S/delta/(n_elements(bluelam)+n_elements(redlam))
	
	errlinebackarr = dblarr(Nx)

        for k=0, Nx-1 do begin
        errlinebackarr[k] = vary0 + (featurelam[k]*sqrt(varm))^2    ; error elements
        endfor

	cnt = m*spec.lambda[roiF]+y0
        plotcnt = m*spec.lambda[roi] + y0

	;print, sigma
	;print, mean(cnt)
	;print, mean(errlinebackarr)

	ewerrarr = dblarr(Nx)

	for k=0, Nx-1 do begin
	ewerrarr[k] = dcorrF[k]/cnt[k]*sqrt( (1/ivarF[k])/dcorrF[k]^2 + errlinebackarr[k]/cnt[k]^2 )*dlf[k]
	endfor	

	fakefluxdens = fltarr(n_elements(dcorrF)) + 3*sigma + cnt			;for every pixel in the defined feature bandpass populate the pixels with flux density equivalent to 3 * rms of the spectrum in addition to the continuum flux density, rms is calculated from the whole bandpass (feature plus background bandpasses)	

	fakeflux = (fltarr(n_elements(dcorrF)) + 3*sigma)*deriv(featurelam)*3e18/featurelam^2*1e-23
        fakelineflux = total(fakeflux)

	ew = -fakelineflux/(mean(cnt)*3e18/mean(featurelam)^2*1e-23)
	ewerr = sqrt(total(ewerrarr))

	print, mean(mean(cnt)*3e18/mean(featurelam)^2*1e-23), fakelineflux

	;if plotew ne 0 then begin
		splot,roilam,dcorrROI*3e18/mean(featurelam)^2*1e-23,xrange=[continuum[0,0],continuum[1,1]],$
			yrange=[min(dcorrROI),max(dcorrROI)],psym=10
		oplot, [bluelam[0],redlam[Nr-1]], [plotcnt[0], plotcnt[Ny-1]],color=1
		oplot, [continuum[0,0],continuum[0,0]],[0,200],color=4
		oplot, [continuum[1,0],continuum[1,0]],[0,200],color=4
		oplot, [continuum[0,1],continuum[0,1]],[0,200],color=4
		oplot, [continuum[1,1],continuum[1,1]],[0,200],color=4
		oplot, [feature[0],feature[0]],[0,200],color=2
		oplot, [feature[1],feature[1]],[0,200],color=2
		oplot, featurelam, fakefluxdens, color=5
	;endif

	avglam = mean(featurelam)
	print, 'The 3 sigma equivalent width of the feature centered at ' + strcompress(string(avglam), /remove_all) + ' Angstroms is ' + strcompress(string(ew), /remove_all) + ' Angstroms'
	print, 'The 3 sigma line flux of the feature centered at ' + strcompress(string(avglam), /remove_all) + ' Angstroms is ' + strcompress(string(fakelineflux), /remove_all) + ' ergs/s/cm2'
	;return, [ew, 0d]
end
;
;
;	continuua is in format [ [[a,b], [c,d]], [[e,f],[g,h]] ]
;	feature is in format [ [a,b], [c,d] ]
