pro measure_fitted_ew, spec, lam1, lam2=lam2, sig1, sig2=sig2, amp1, amp2=amp2, order

	if order ne 0 and order ne 1 then begin
                message, 'Continuum can only be fit to 0th or 1st order polynomial, change order'
        endif

	m = mrdfits(spec,1)

	range = where(m.lambda gt lam1-250 and m.lambda lt lam1+250)
	mockarray = findgen(n_elements(range))	

	if n_elements(lam2) and n_elements(sig2) and n_elements(amp2) then begin

	if order eq 0 then expr = 'P[0] + GAUSS1(X, P[1:3]) + GAUSS1(X,P[4:6])'$
		else expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4]) + GAUSS1(X,P[5:7])'

	area1 = sig1*amp1*sqrt(2*!PI)
	area2 = sig2*amp2*sqrt(2*!PI)

	if order eq 0 then start = [5e-05, lam1, sig1, area1, lam2, sig2, area2]$
		else start = [5e-05, 5e-05, lam1, sig1, area1, lam2, sig2, area2]
	
	result = MPFITEXPR(expr,m.lambda[range], m.spec[range], sqrt(1/m.ivar[range]), start, perr=perr, covar=covar)
	
	if order eq 0 then begin

		area = result[3] 
		lam = result[1]
		sig = result[2]
		areaerr = perr[3]
		lamerr = perr[1]
		sigerr = perr[2]
		background_fnu = mean(result[0]*3e18*1e-23/m.lambda[range]^2)
		errback_fnu = sqrt(mean( (perr[0]*3e18*1e-23/m.lambda[range]^2)^2 ) )

	endif else begin
	
		area = result[4]
		lam = result[2]
		sig = result[3]
		areaerr = perr[4]
		lamerr = perr[2]
		sigerr = perr[3]
		background_fnu = mean(result[0]*3e18*1e-23/m.lambda[range]^2+m.lambda[range]*result[1]*3e18*1e-23/m.lambda[range]^2)	;in Jy
		
		errback_fnu = sqrt(mean( (perr[0]*3e18*1e-23/m.lambda[range]^2)^2 + m.lambda[range]*(perr[1]*3e18*1e-23/m.lambda[range]^2)^2 + 2*m.lambda[range]*(covar[0,1]*(3e18*1e-23/m.lambda[range]^2))^2 ) )	;is sigma_F = sqrt(a^2*simga_X^2 + b^2*sigma_Y^2 + 2ab*sigma_XY^2), units are annoying but correct
	endelse

	;bestamp = area/(sqrt(2*!PI)*sig)
	;bestamperr = area/(sqrt(2*!PI)*sig)*sqrt( (sigerr/sig)^2 + (areaerr/area)^2 ) 
	
	; for now just going to take the mean continuum value in the entire range and divide that into the mean
	; Jy*Ang value from the line fit, should instead compute the Gaussian at each point and subtract off
	; the fitted background value over the entire range
	; ACTUALLY MAYBE NOT SINCE IT IS A FLUX DENSITY, JUST TAKING AVERAGE OVER THE ENTIRE RANGE, PROLLY OK

	;linerange = where(m.lambda ge lam - 3*sig and m.lambda le lam + 3*sig)	;aprroximating line bandpass assuming the Gaussian drops to the background level
										;at +/- 3 sigma from the line centroid. This underestimates the background and will
										;overestimate the (magnitude of the) EW

	fnu_area  = area*3e18*1e-23/lam^2		;in Jy*Ang, area output from MPFITEXPR is already background subtracted
	fnu_area_err = areaerr*3e18*1e-23/lam^2	

	ew = -fnu_area/background_fnu
	;print, fnu_area_err/fnu_area, errback_fnu/background_fnu
	;print, fnu_area_err, fnu_area, errback_fnu,background_fnu
	ewerr = ew*sqrt( (fnu_area_err/fnu_area)^2 + (errback_fnu/background_fnu)^2 )
	
	print, 'The Equivalent Width of the Gaussian centered at ' + strcompress(string(lam),/remove_all) + ' is ' + strcompress(string(ew),/remove_all) + ' +/- ' + strcompress(string(ewerr),/remove_all)

	splot, m.lambda[range], m.spec[range], psym=10, yrange=[min(m.spec[range])-mean(m.spec[range]), max(m.spec[range])+mean(m.spec[range])]
	
	 if order eq 0 then oplot, m.lambda, result[0] + gauss1(m.lambda, result[1:3]), color=1, thick=3$
                else  oplot, m.lambda, result[0] + result[1]*m.lambda + gauss1(m.lambda, result[2:4]), color=1, thick=3

	endif
	
	
	if n_elements(lam2) lt 1 then begin

        if order eq 0 then expr = 'P[0] + GAUSS1(X, P[1:3])'$
                else expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4])'

        area1 = sig1*amp1*sqrt(2*!PI)

        if order eq 0 then start = [5e-05, lam1, sig1, area1]$
                else start = [5e-05, 5e-05, lam1, sig1, area1]

        result = MPFITEXPR(expr,m.lambda[range], m.spec[range], sqrt(1/m.ivar[range]), start, perr=perr, covar=covar)

        if order eq 0 then begin
        
                area = result[3]
                lam = result[1]
                sig = result[2]
                areaerr = perr[3]
                lamerr = perr[1]
                sigerr = perr[2]
		background_fnu = mean(result[0]*3e18*1e-23/m.lambda[range]^2)
                errback_fnu = sqrt(mean( (perr[0]*3e18*1e-23/m.lambda[range]^2)^2))

        endif else begin 
         
                area = result[4]
                lam = result[2]
                sig = result[3]
                areaerr = perr[4]
                lamerr = perr[2]
                sigerr = perr[3]
		background_fnu = mean(result[0]*3e18*1e-23/m.lambda[range]^2+m.lambda[range]*result[1]*3e18*1e-23/m.lambda[range]^2)    ;in Jy
		errback_fnu = sqrt(mean( (perr[0]*3e18*1e-23/m.lambda[range]^2)^2 + m.lambda[range]*(perr[1]*3e18*1e-23/m.lambda[range]^2)^2 + 2*m.lambda[range]*(covar[0,1]*(3e18*1e-23/m.lambda[range]^2))^2 ) )
				

        endelse
		
	fnu_area  = area*3e18*1e-23/lam^2          ;in Jy*Ang, area output from MPFITEXPR is already background subtracted
        fnu_area_err = areaerr*3e18*1e-23/lam^2

        ew = -fnu_area/background_fnu
        ewerr = ew*sqrt( (fnu_area_err/fnu_area)^2 + (errback_fnu/background_fnu)^2 )

        print, 'The Equivalent Width of the Gaussian centered at ' + strcompress(string(lam),/remove_all) + ' is ' + strcompress(string(ew),/remove_all) + ' +/- ' + strcompress(string(ewerr),/remove_all)

        splot, m.lambda[range], m.spec[range], psym=10, yrange=[min(m.spec[range])-mean(m.spec[range]), max(m.spec[range])+mean(m.spec[range])]

         if order eq 0 then oplot, m.lambda, result[0] + gauss1(m.lambda, result[1:3]), color=1, thick=3$
                else  oplot, m.lambda, result[0] + result[1]*m.lambda + gauss1(m.lambda, result[2:4]), color=1, thick=3

        endif


end

