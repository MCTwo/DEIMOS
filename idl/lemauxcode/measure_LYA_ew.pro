pro measure_LYA_ew, spec, shortlimit, longlimit, mag, filter, z, fluxerr, magerr

	if n_elements(magerr) then magerr=magerr[0] else magerr = 0.
	spec = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/' + spec
	mm =fill_gap(spec,/tweak,/telluric,/silent,header=header)
	wave=mm.lambda
	airtovac,wave
	mm.lambda=wave

	;restspeclam = mm.lambda/(1+z[0])
	
	if filter eq 'r' or filter eq 'i' or filter eq 'z' then begin
	readcol, filter + '.dat', format = 'D,D,D', lam, pointthrough, extthrough, /silent
	through = (pointthrough+extthrough)/2
	endif
	
	if filter ne 'r' or filter ne 'i' or filter ne 'z' then begin
	readcol, filter + '.dat', format = 'D,D', lam, through, /silent
	endif
	;restlam = lam
	minlam = min(lam)       ;min wavelength of the filter transmission curve
	maxlam = max(lam)
	
	print, 'The bluest wavelength this filter covers is ', strcompress(string(minlam), /REMOVE_ALL), ' and the reddest is ', strcompress(string(maxlam), /REMOVE_ALL)
	print, ' ' 

	dummy = where(mm.lambda ge minlam and mm.lambda le maxlam)
	arraylength = n_elements(dummy)
	lambda_in = mm.lambda[dummy]
	spec_in = mm.spec[dummy]
	mocklam = minlam + findgen(arraylength)*(maxlam-minlam)/arraylength
	througharr = interpol(through, lam, mocklam, /quadratic)
	
	;set all values of through array equal to zero blueward of the line
	dummyer = where(mocklam lt shortlimit[1])
	if dummyer[0] eq -1 then begin
		message, 'The Lyman alpha line is not in the filter bandpass, change to a different filter!'
	endif
		
	truncthrougharr = througharr
	truncthrougharr[dummyer] = 0
	
	deltalamband = deriv(lambda_in)	;the lambda spacing in each bin
        meanlam = fltarr(arraylength)

	for j=0, arraylength-2 do begin
                meanlam[j] = (lambda_in[j] + lambda_in[j+1])/2    ;get the effective wavelength in each bin
        endfor


        dummy3 = where(meanlam eq 0)		;clean up last element
        meanlam[dummy3] = max(meanlam)+mean(deltalamband)
	
	fakeavgflam = fltarr(150,n_elements(dummy))
	gridsize = size(fakeavgflam, /dimensions)

	for i=0, gridsize[0]-1 do begin	
	fakeavgflam[i,*] = 10^(-21 + i*0.02)	;generate a grid of possible average flux densities, can change the limits and the grdding size or courseness here
	endfor	

	
	
	bandfluxnum = fltarr(gridsize[0],n_elements(dummy))
        bandfluxdenom = fltarr(n_elements(dummy))
	aquantity = fltarr(gridsize[0])
	aquantitylowerr = fltarr(gridsize[0])
	aquantityhigherr = fltarr(gridsize[0])
	
	for i=0, gridsize[0]-1 do begin		
		for j=0, n_elements(dummy)-1 do begin
                	bandfluxnum[i,j] = througharr[j]*fakeavgflam[i,j]*deltalamband[j]*meanlam[j]/(3.0D+18)    ;this is pedantic see Fukugita 1995 eqn 9, F_nu is corr_spec_in/deltafband
                	bandfluxdenom[j] = througharr[j]*deltalamband[j]/meanlam[j]
		endfor
		;print, total(bandfluxnum[i,*]), total(bandfluxdenom), 10^(-(mag+48.60)/2.5)	
		aquantity[i] = 10^(-(mag+48.60)/2.5)*total(bandfluxdenom)-total(bandfluxnum[i,*])	;solving for the average flux density
		aquantitylowerr[i] = 10^(-(mag+magerr+48.60)/2.5)*total(bandfluxdenom)-total(bandfluxnum[i,*])
		aquantityhigherr[i] = 10^(-(mag-magerr+48.60)/2.5)*total(bandfluxdenom)-total(bandfluxnum[i,*])
	endfor
	
	thats = where(abs(aquantity) eq min(abs(aquantity)))
	specious = where(abs(aquantitylowerr) eq min(abs(aquantitylowerr)))
	reasoningdad = where(abs(aquantityhigherr) eq min(abs(aquantityhigherr)))
	fakemag = -2.5*alog10(total(bandfluxnum[thats,*])/total(bandfluxdenom)) - 48.60
	
	idealflam = fakeavgflam[thats]
	;print, fakeavgflam[specious], fakeavgflam[reasoningdad]
	idealflamlowerr = idealflam - fakeavgflam[specious] 	;the errors are calculated by taking the difference of the new flux density calculated from the 1 sigma upper and lower bounds of the magnitude and and the flux density calculated by the actual magnitude
	idealflamhigherr = fakeavgflam[reasoningdad] - idealflam
	
	strmag = strcompress(string(mag),/REMOVE_ALL)
	strfakemag = strcompress(string(fakemag),/REMOVE_ALL) 
	stridealflam = strcompress(string(idealflam),/REMOVE_ALL)
	print, 'The magnitude of the object was measured as ', strmag, ' and the best fit magnitude was ', strfakemag
	print, 'This best fit magnitude corresponds to an average flux density of ', stridealflam[0], ' ergs/s/cm^2/', string(197B)
	print, ' '
	print, ' '
 
	;calculating background, using calculation from line flux measurement program
	bluebandpass = where(mm.lambda ge shortlimit[0] and mm.lambda le longlimit[0])
	redbandpass = where(mm.lambda ge shortlimit[2] and mm.lambda le longlimit[2])
	featurebandpass = where(mm.lambda ge shortlimit[1] and mm.lambda le longlimit[1])
	blam = mm.lambda[bluebandpass]
	rlam = mm.lambda[redbandpass]
	flam = mm.lambda[featurebandpass]
	ivarB = mm.ivar[bluebandpass]
	ivarR = mm.ivar[redbandpass]	
	bluespec = mm.spec[bluebandpass]
	redspec = mm.spec[redbandpass]

	S = total(ivarB) + total(ivarR)
        Sx = total(blam*ivarB) + total(rlam*ivarR)
        Sxx = total(blam^2*ivarB) + total(rlam^2*ivarR)
        Sy = total(bluespec*ivarB) + total(redspec*ivarR)
        Sxy = total(bluespec*blam*ivarB) + total(redspec*rlam*ivarR)
        delta = Sxx*S - Sx^2

	y0 = (Sxx*Sy - Sx*Sxy)/delta  
        m = (Sxy*S - Sx*Sy)/delta

	;m = (mean(mm.spec[redbandpass])-mean(mm.spec[bluebandpass]))/(mean(rlam)-mean(blam))
	;b = mm.lambda[bluebandpass[0]]*(-m) 
	cnt = m*flam + y0
	plotcnt = m*mm.lambda + y0
	
	
	constfluxcorr = 2e-08/(3600.*!pi*(996./2)^2*0.33)
	bkgdsubspec = mm.spec[featurebandpass]-cnt
	bkgdfluxdens = cnt*constfluxcorr/mm.lambda[featurebandpass]
	bkgdsubfluxdens = bkgdsubspec*constfluxcorr/mm.lambda[featurebandpass]
	;print, mean(bkgdsubfluxdens)
	magcontflam = fltarr(n_elements(featurebandpass)) + idealflam[0]
	dlam = deriv(flam)	
	;print, mean(bkgdfluxdens)

	obsew = total((1-bkgdsubfluxdens/magcontflam)*dlam)
	obsdawsonew = total(bkgdsubfluxdens*0.33)/idealflam
	;print, total(bkgdsubfluxdens*0.33)	
	ewdaw = obsdawsonew/(1+z)
	errdaw = fluxerr/idealflam/(1+z)	;use this error for negative error in cases when the galaxy is not detected in the photometry 	
	lowerrdaw = sqrt((fluxerr/total(bkgdsubfluxdens*0.33))^2 + (idealflamhigherr/idealflam)^2)*ewdaw	;the ideal flux density high error refers to the higher flux density limit, so this is the error on the lower side  
	higherrdaw = sqrt((fluxerr/total(bkgdsubfluxdens*0.33))^2 + (idealflamlowerr/idealflam)^2)*ewdaw
	ew = obsew/(1+z)	
	;print, mean(bkgdsubfluxdens)
	print, 'The standard measurement of Equivalent Width of this line is ', strcompress(string(ew), /REMOVE_ALL), ' ', string(197B)
	print, 'The measurement of Equivalent Width as defined in Dawson et al. is ', strcompress(string(ewdaw[0]), /REMOVE_ALL), ' ', string(197B)
	if magerr ne 0. then print, 'with associated 1 sigma errors + ', strcompress(string(higherrdaw[0]), /REMOVE_ALL), ' - ', strcompress(string(lowerrdaw[0]), /REMOVE_ALL), ' ', string(197B)
	if magerr eq 0. then print, 'with associated 1 sigma error ', strcompress(string(errdaw[0]), /REMOVE_ALL), ' ', string(197B)
	loadcolors
        splot, mm.lambda[bluebandpass[0]-20:redbandpass[n_elements(redbandpass)-1]+20], mm.spec[bluebandpass[0]-20:redbandpass[n_elements(redbandpass)-1]+20]
        oplot, mm.lambda, plotcnt, thick=2, color=5
	oplot, [blam[0],blam[0]], [-2000,2000], color=3
	oplot, [blam[n_elements(blam)-1], blam[n_elements(blam)-1]], [-2000,2000], color=3
	oplot, [rlam[0],rlam[0]], [-2000,2000], color=4
        oplot, [rlam[n_elements(rlam)-1], rlam[n_elements(rlam)-1]], [-2000,2000], color=4
	oplot, [flam[0],flam[0]], [-2000,2000], color=6
        oplot, [flam[n_elements(flam)-1], flam[n_elements(flam)-1]], [-2000,2000], color=6
	
end 
	
		
