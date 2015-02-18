
function biweight, xi
	
	M = 10^double(median(alog10(xi)))
	MAD = median(abs(xi-M))
	ui = (xi-M)/(MAD*6d)

	r = where(abs(ui) lt 1d, count)
	if count lt 1 then return, M
	Cbi = M + total((xi[r]-M)*(1-ui[r]^2)^2)/total(1-ui[r]^2)^2
	return, cbi
end

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
function measure_simple_ew, $
	spec, $		; The spectrum structure
	bm=bm, $	; Best model spectrum
	continuum, $	; The wavelength range for the continuum [[l0,l1],[l2,l3]]
	feature,$		; The wavelength range for the feature
	quiet=quiet

	usemodel=1
	if keyword_set(bm) eq 0 then usemodel = 0

	fail = [!values.F_NAN, !values.F_NAN]
	roi = where(spec.lambda ge feature[0] and spec.lambda le feature[1], count)
	
	if keyword_set(quiet) then q=1 else q=0
	if count le 2 then begin
		if q eq 0 then $
		message, "Wavelength range of feature is non existent", /info
		return, fail
	endif
	dlf = deriv(spec.lambda[roi])
	ffs = 0
	if usemodel eq 1 then begin
		bms = interpol(bm.spec, bm.lambda, spec.lambda[roi], /spline)
		bms = bms/median(bms)*median(spec.spec[roi])
	endif
	print, usemodel
	roi1 = where(spec.lambda ge continuum[0,0] and spec.lambda le continuum[1,0], count1) 
	roi2 = where(spec.lambda ge continuum[0,1] and spec.lambda le continuum[1,1], count2)

	if count1 le 2 or count2 le 2 then begin	
		if q eq 0 then $
		message, "Wavelength range of continuum is non existent", /info
		return, fail
	endif
	l1 = mean(spec.lambda[roi1])
	l2 = mean(spec.lambda[roi2])

	y1 = biweight(spec.spec[roi1])
	y2 = biweight(spec.spec[roi2])
	if usemodel eq 1 then begin
		bms = interpol(bm.spec, bm.lambda, spec.lambda[roi], /spline)
		bms = bms/median(bms)*median(spec.spec[roi])
	endif
	
	m = (y2-y1)/(l2-l1)
	cnt = m*spec.lambda[roi]+y1-m*l1
	;print, cnt
	;ss = min([[spec.spec[roi]],[spec.spec[roi+1]]],dim=2)
	ew = total((1d - spec.spec[roi]/cnt)*dlf,/nan,/double)

	;if plotew ne 0 then begin
		splot,spec.lambda,spec.spec,xrange=[continuum[0,0],continuum[1,1]],$
			yrange=[0,max(spec.spec[roi])],psym=10
		oplot, [l1,l2],[y1,y2],color=1
		oplot, [continuum[0,0],continuum[0,0]],[0,200],color=4
		oplot, [continuum[1,0],continuum[1,0]],[0,200],color=4
		oplot, [continuum[0,1],continuum[0,1]],[0,200],color=4
		oplot, [continuum[1,1],continuum[1,1]],[0,200],color=4
		oplot, [feature[0],feature[0]],[0,200],color=2
		oplot, [feature[1],feature[1]],[0,200],color=2
	;endif

	return, [ew, 0d]
	return, [(fcm-ffm)/fcm*dlf, $
		sqrt((ffm/fcm^2*dlf*fcs)^2 + (dlf/fcm*ffs)^2)]
end
;
;
;	continuua is in format [ [[a,b], [c,d]], [[e,f],[g,h]] ]
;	feature is in format [ [a,b], [c,d] ]
function measure_simple_ews, spec, continuua, features
	
	n = (size(continuua,/dim))[2]
	nb = (size(features,/dim))[1]
	if n ne nb then begin
		message, 'The number of continuua features is ', n, $
			' the number of spectral features is ', nb, $
			' these numbers must match or else something is wrong'
		return, [!values.F_NAN]
	endif

	v = dblarr(2,n)
	for i = 0, n-1 do begin
		v[*,i] = measure_simple_ew(spec, continuua[*,*,i], features[*,i])
	endfor

	return, v
end

pro measure_simple_ew_test
	spec = {lambda: findgen(10000), spec: fltarr(10000)}
	spec.spec[50:60] = 1
	spec.spec = spec.spec+1

	print, 'The results of the following ratio should be 2'
	print, measure_simple_Ratio(spec, [30,40], [50,60])
	spec.spec = spec.spec*2
	print, measure_simple_Ratio(spec, [30,40], [50,60])
	spec.spec = spec.spec*.5

	print, 'The results of the following should be -11 angstroms (emission)'
	;plot,spec.spec,xrange=[40,70],yrange=[0,3]
	print, measure_simple_ew(spec, [[40,45],[65,67]],[50,60])
	spec.spec=spec.spec*2
	print, measure_simple_ew(spec, [[40,45],[65,70]],[50,60])
	print, measure_simple_ew(spec, [[30,45],[65,80]],[50,60])
	print, measure_simple_ew(spec, [[30,45],[65,66]],[50,60])

	print,'==========='
	spec = {lambda: findgen(10000), spec: fltarr(10000)}
	spec.spec[54:56] = 5
	spec.spec = spec.spec+1
	print, measure_simple_ew(spec, [[40,45],[65,70]],[50,60])
	print,'==========='

	;splot,spec.lambda,spec.spec,xrange=[30,80],yrange=[0,5]

	spec.spec = ((spec.spec-1.)+spec.spec)/2.
	print, 'Adding different spectra together'
	print, measure_simple_ew(spec, [[30,45],[65,66]],[50,60])

	print, 'Need to figure out what the results of the following should be'
	l = findgen(1000.)/10.+3700.
	c = [0., 0., 1., 1., 1.]
	ss = {lambda: l, spec: oii_function(l, c)/total(oii_function(l,c)*deriv(l))+1d}
	
	cs = dindgen(10)/50
	cs = cs+cs[1]
	ews = dblarr(n_elements(cs))
	print,'========================================='
	for i = 0, n_elements(cs)-1 do begin
		c = [0., 0., cs[i], 1., 1.]
		ss = {lambda: l, spec: oii_function(l, c)/total(oii_function(l,c)*deriv(l))+1d}
		ews[i] = (measure_simple_ew_line(ss,'oii'))[0]
	endfor
	plot,cs,ews,/xlog,xrange=[1d-3,1d3],thick=2,yrange=[0,-1]
	print,'========================================='

	ss.spec = ss.spec*2
	print,measure_simple_ews(ss,[[[3700,3710],[3740,3750]], [[3750,3770],[3780,3790]]], $
		[[3710,3740], [3770,3780]])
	print, 'adding some noise'
	ss.spec = ss.spec+randomn(seed,n_elements(ss.spec))
	;plot,ss.lambda,ss.spec,yrange=[0,6]
	print,measure_simple_ews(ss,[[[3700,3710],[3740,3750]], [[3750,3770],[3780,3790]]], $
		[[3710,3740], [3770,3780]])
	
	;splot,ss.lambda,ss.spec,xrange=[3700,3750],yrange=[0,8]
end
