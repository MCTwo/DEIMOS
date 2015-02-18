
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
pro measure_bandpass_ew, $
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
	
	goodpixB = where(spec.ivar[roiB] gt 1e-08)              ;changed 4/7 to remove fill gap areas where the ivar is set to 1e-30, could set this higher if wanted to exclude bad sky areas, CHANGE IF NOT USING DEIMOS DATA!!
        goodpixR = where(spec.ivar[roiR] gt 1e-08)
        
	goodpixF = where(spec.spec[roiF] ne 0.0000)             ;I don't cut the feature based on cuz that seems not right, that uncertainty is real. I am just removing the pixels on the ends of the spectrum which were not fill gapped properly. I select in the blueward and redward bandpasses because the continuum is really not that uncertain. 

        totivarB = spec.ivar[roiB]                              
        totivarR = spec.ivar[roiR]
        if goodpixB[0] ge 0 then ivarB = totivarB[goodpixB] else ivarR = totivarR             ;only taking the pixels with reasonable variances, also checks to make sure there is some good portion of the spectrum, if not then it just calculates it with the high variance parts included
        if goodpixR[0] ge 0 then ivarR = totivarR[goodpixR] else ivarR = totivarR
        totivarF = spec.ivar[roiF]
        ivarF = totivarF[goodpixF]

        totbluelam = spec.lambda[roiB]
        totredlam = spec.lambda[roiR]
        if goodpixB[0] ge 0 then bluelam = totbluelam[goodpixB] else bluelam = totbluelam
        if goodpixR[0] ge 0 then redlam = totredlam[goodpixR] else redlam = totredlam
        totfeaturelam = spec.lambda[roiF]
        featurelam = totfeaturelam[goodpixF]
        roilam = spec.lambda[roi]

        totbluespec = spec.spec[roiB]
        totredspec = spec.spec[roiR]
        if goodpixB[0] ge 0 then bluespec = totbluespec[goodpixB] else bluespec = totbluespec
       	if goodpixR[0] ge 0 then redspec = totredspec[goodpixR] else redspec = totredspec
        totfeaturespec = spec.spec[roiF]   
        featurespec = totfeaturespec[goodpixF]
        roispec = spec.spec[roi]     

          ; What I am going to do is set the ivar to Poissonian where the fill_gap is IF it is in the feature bandpass, it's not right to remove the pixels entirely, but this is the next best thing. Otherwise, any EW measurement where the gap was filled has infinite error and realistically, it's known better than that.

        gapped_feature = where(spec.ivar[roiF] lt 1e-10)

        if gapped_feature[0] ge 0. then begin

                ivarF[gapped_feature] = 1/featurespec[gapped_feature]^2
	endif

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

	ewerrarr = dblarr(Nx)
	
	print, n_elements(cnt)
	print, Nx
	
	for k=0, Nx-1 do begin
	ewerrarr[k] = dcorrF[k]/cnt[k]*sqrt( (1/ivarF[k])/dcorrF[k]^2 + errlinebackarr[k]/cnt[k]^2 )*dlf[k]
	endfor	
	
	ew = total((1d - dcorrF/cnt)*dlf,/nan,/double)
	ewerr = sqrt(total(ewerrarr))


	;if plotew ne 0 then begin
		splot,roilam,dcorrROI,xrange=[continuum[0,0],continuum[1,1]],$
			yrange=[0,max(dcorrROI)],psym=10
		oplot, [bluelam[0],redlam[Nr-1]], [plotcnt[0], plotcnt[Ny-1]],color=1
		oplot, [continuum[0,0],continuum[0,0]],[0,200],color=4
		oplot, [continuum[1,0],continuum[1,0]],[0,200],color=4
		oplot, [continuum[0,1],continuum[0,1]],[0,200],color=4
		oplot, [continuum[1,1],continuum[1,1]],[0,200],color=4
		oplot, [feature[0],feature[0]],[0,200],color=2
		oplot, [feature[1],feature[1]],[0,200],color=2
	;endif

	avglam = (feature[0] +  feature[1])/2
	print, 'The Equivalent Width of the feature centered at ' + strcompress(string(avglam), /remove_all) + ' Angstroms is ' + strcompress(string(ew), /remove_all) + ' +/- ' + strcompress(string(ewerr),/remove_all) + ' Angstroms'
	;return, [ew, 0d]
end
;
;
;	continuua is in format [ [[a,b], [c,d]], [[e,f],[g,h]] ]
;	feature is in format [ [a,b], [c,d] ]
