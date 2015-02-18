pro manual_skysub_DEIMOS, sigma
		
	;restore,  file='/Users/blemaux/cvs/spec2d/etc/template.sav'
	mm = mrdfits('/Users/blemaux/cvs/spec2d/etc/uves_sky.fits',1)	
	temp_flux = mm.spec
	temp_wave = mm.lambda

	halfwidth = floor(1.5*sigma)
	kernel=findgen(2*halfwidth+1)-halfwidth
	kernel=exp(-kernel^2/2/sigma^2)
	kernel=kernel/total(kernel)
	stemp_flux =  convol(temp_flux,  kernel,  /center)
	
	marray = findgen(20)
	barray = findgen(40)*0.1 - 2
	initcounts = 10000
	m = 0
	b = 0	

	
	for i=0, n_elements(barray)-1 do begin
		;file = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/16XR1/*/spec1d.16XR1.019.COS_SC2_01188.fits'	
		;spec = fill_gap(file,/tweak, /telluric,/silent,header=header)
		;wave=spec.lambda
        	;airtovac,wave   
        	;spec.lambda=wave
		spec = mrdfits('/Volumes/Data2/orelse/lemaux/deimos/sc1604/cluster_catalog/AGNspectra/spec1d.16XR1.019.COS_SC2_01188.fillgaped.fits',1, /SILENT)
 	       	minlam = min(spec.lambda)
 	       	maxlam = max(spec.lambda)
		;zeros = where(spec.spec eq 0)
        	nonzeros = where(spec.spec ne 0)
		
		range = where(temp_wave ge minlam and temp_wave le maxlam)
        	rangelam = temp_wave[range]+barray[i]
        	rangesky = stemp_flux[range]
        	smoothskyflux = interpol(rangesky, rangelam, spec.lambda)
		
		for j=0, n_elements(marray)-1 do begin
		
		spec.spec = spec.spec -  marray[j]*smoothskyflux
		counts = n_elements(where(spec.spec gt 400))
		if counts lt initcounts then m = marray[j] else m=m
		if counts lt initcounts then b = barray[i] else b=b
		if counts lt initcounts then initcounts = counts else initcounts=initcounts	
		endfor 

	endfor
	
	finalrangelam = temp_wave[range]-b[0]
	finalskyflux = interpol(rangesky*m[0], rangelam, spec.lambda)
	print, m, b

        spec = mrdfits('/Volumes/Data2/orelse/lemaux/deimos/sc1604/cluster_catalog/AGNspectra/spec1d.16XR1.019.COS_SC2_01188.fillgaped.fits',1, /SILENT)	
	spec.spec = spec.spec - finalskyflux
	ss1d = spec.spec
	ss1dlam = spec.lambda	
	ss1divar = spec.ivar

	splot, ss1dlam, ss1d
	;splot, ss1dlam, ss1d
	result = {lambda:ss1dlam,spec:ss1d, ivar:ss1divar}

	mwrfits, result, 'skysubtracted.16XR1.019.fits'
	openw, lun, '16XR1.019.COS_SC2_01188.skysubtracted.dat', /get_lun

                for i=0, n_elements(ss1dlam)-1 do begin
                printf, lun, ss1dlam[i], '   ', ss1d[i], format='(d10.5, 3a, d9.5)'
                endfor

        free_lun, lun		
	;mwrfits, result, 'skysubtracted.16XR1.16.fits'
	
end	
