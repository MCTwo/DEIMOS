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

function Ez, x

   distance = 1/sqrt(0.3*(1+x)^3+0.7) ;concordance cosmology assumtion omega_M = 0.3, omega_lam = 0.7, to be consistant with DEEP2
   return, distance

end

pro get_lineflux_NIRSPEC, cat, outcat

	readcol, cat, format='A,A,D,D,D,D,D,D,D,I,A', masklist, slitlist, bluelist1, bluelist2, redlist1, redlist2, linelist1, linelist2, zlist, q, file
	openw, lun, outcat, width=200, /get_lun
	tab = '09'XB

	for n=0, n_elements(q)-1 do begin
	mask = masklist[n]
	slit = slitlist[n]
	blueside = [bluelist1[n], bluelist2[n]]
	redside = [redlist1[n], redlist2[n]]
	feature = [linelist1[n],linelist2[n]]
	z = zlist[n]
	spectrum = file[n]	
	spectrum = '/Volumes/Data2/orelse/lemaux/NIRSPEC/06042007/2007jun04/reduced/final_spec/'+spectrum

	spec=mrdfits(spectrum,1)

	print, spectrum	
	; OLD section of the program, I had split up the 1d spectrum into red and blue sides and tried to figure out which CCD the emission line was on, generalized to the entire spectrum through fill_gap instead
	
	;if blueside[1] gt specBhorne.lambda[4095] then spec = specRhorne else spec = specBhorne ;
	;specBhorne = mrdfits(spectrum,3)
        ;specRhorne = mrdfits(spectrum,4)
	;if abs(blueside[1]-specBhorne.lambda[4095]) lt 50. then begin	;if it's close to the CCD gap use the whole spectrum, kinda sloppy...
	;s = dblarr(8192)
	;i = dblarr(8192)
	;l = dblarr(8192)
	;s[0:4095] = specBhorne.spec
	;s[4096:8191] = specRhorne.spec
	;i[0:4095] = specBhorne.ivar
        ;i[4096:8191] = specRhorne.ivar
	;l[0:4095] = specBhorne.lambda
        ;l[4096:8191] = specRhorne.lambda
	;spec = {spec:s, ivar:i, lambda:l}
	;endif 

	bluecontarr = where(spec.lambda ge blueside[0] and spec.lambda le blueside[1])
	redcontarr = where(spec.lambda ge redside[0] and spec.lambda le redside[1])
	linearr = where(spec.lambda ge feature[0] and spec.lambda le feature[1])	
	fullarr = where(spec.lambda ge blueside[0] and spec.lambda le redside[1])
	fullband = spec.lambda[fullarr]
	dlamfull = deriv(spec.lambda[fullarr])
	blueband = spec.lambda[bluecontarr]
	dlamblue = deriv(spec.lambda[bluecontarr])
	redband = spec.lambda[redcontarr]
	dlamred = deriv(spec.lambda[redcontarr])
	lineband = spec.lambda[linearr]
	dlamline = deriv(spec.lambda[linearr])
	bluecontcounts = spec.spec[bluecontarr]*1e-23
	redcontcounts = spec.spec[redcontarr]*1e-23
	linecounts = spec.spec[linearr]*1e-23
	speccounts = spec.spec*1e-23	
	ivarline = spec.ivar[linearr]

	Ny = n_elements(fullband)
	Nx = n_elements(spec.spec)	
	Ni = n_elements(bluecontarr)
	Nj = n_elements(redcontarr)
	Nk = n_elements(linearr)
	dcorrB = dblarr(Ni)
	dcorrR = dblarr(Nj)
	a = dblarr(Nj)
	b = dblarr(Nj)
	dcorrline = dblarr(Nk)
	dcorrlineerr = dblarr(Nk)
	dcorrspec = dblarr(Nx)

	dcorrB = bluecontcounts		
		
	dcorrR = redcontcounts


	dcorrline = linecounts
	dcorrspec = speccounts
	dcorrlineerr = 1/sqrt(ivarline)*1e-23
	;for k=0, Nk-1 do begin
	;dcorrline[k] = deimos_correction(lineband[k])*linecounts[k] ; counts/exptime/aperture/pix in line bandwidth
	;if lineband[k] ge 0 then dcorrlineerr[k]= sqrt(dcorrline[k]) else dcorrlineerr[k] = 0. ; Poisson error counts/exptime/aperture/pix, giving 0 error to negative counts...valid?!?
	;dcorrlineerr[k] = 1/sqrt(ivarline[k])*deimos_correction(lineband[k])
	;dcorrlineerr[k] = sqrt(dcorrlineerr[k]^2+ 1/ivarline[k]) 
	;endfor

	dcorrpoisfluxvar = dblarr(Nk)
	for k=0, Nk-2 do begin
	dcorrpoisfluxvar[k] = (dcorrlineerr[k]*dlamline[k]*3e18/(lineband[k]^2))^2 	;calculating Poisson vairance in flux units erg/s/cm^2 to add in quadrature across the line bandpass
	endfor

	print, mean(dcorrline), mean(dcorrlineerr)
	dcorrpoisfluxerr = sqrt(total(dcorrpoisfluxvar));*(lineband[Nk-1]-lineband[0]) ;Poisson flux error is specific flux, need to multiply by line bandpass ### is old, no longer multiplying by bandapss	

	ivarB = 1/(sqrt(1/spec.ivar[bluecontarr])*1e-23)^2	;1/(dblarr(n_elements(blueband))+errbackB)^2
	ivarR = 1/(sqrt(1/spec.ivar[redcontarr])*1e-23)^2   ;1/(dblarr(n_elements(redband))+errbackR)^2	
	
	S = total(ivarB) + total(ivarR)
	Sx = total(blueband*ivarB) + total(redband*ivarR)
	Sxx = total(blueband^2*ivarB) + total(redband^2*ivarR)
	Sy = total(dcorrB*ivarB) + total(dcorrR*ivarR)
	Sxy = total(dcorrB*blueband*ivarB) + total(dcorrR*redband*ivarR)
	delta = Sxx*S - Sx^2	

	;fakem = (dcorrR[Nj-1] - dcorrB[0])/(redband[Nj-1]-blueband[0])
	;fakey0 = dcorrB[0] - fakem*blueband[0]
	
	y0 = (Sxx*Sy - Sx*Sxy)/delta  
	m = (Sxy*S - Sx*Sy)/delta
	vary0 = Sxx/delta/(n_elements(blueband)+n_elements(redband))*(dlamline*3e18/(lineband^2))^2	;is variance divided by sample size, since it is variance correction is squared
	varm = S/delta/(n_elements(blueband)+n_elements(redband))*(dlamline*3e18/(lineband^2))^2
	;sigb0 = stddev(dcorrB)/sqrt(n_elements(dcorrB)) 	; variance on y0 too large to work correctly, this is also valid since background value in blue bandpass is Gaussian 
	cnt = m*lineband+y0
	fullcnt = m*fullband+y0
	lineback = total(cnt)	
	;fakecnt = fakem*fullband+fakey0
	
	;background subtracting continuum counts from total counts at each pixel
	
	bkgdsubdcorrline = dblarr(Nk)
	for k=0, Nk-1 do begin
	bkgdsubdcorrline[k] = dcorrline[k] - cnt[k]
	endfor 

	print, mask, slit
	;converting to flux units in erg/s/cm^2
	totallinefluxarr = dblarr(Nk)
	for k=0, Nk-1 do begin
	totallinefluxarr[k] = double(bkgdsubdcorrline[k]*dlamline[k]*3e18/(lineband[k]^2))
	endfor
	
	totallineflux = total(totallinefluxarr);*(lineband[Nk-1]-lineband[0]) ;totallinefluxarr is a specific flux in erg/s/cm^2/A, have to multiply by line bandpass ###This is old, calculated a flux by removing 0.33 Ang/pix in constfluxcorr
	
	errlinebackarr = dblarr(Nk)
	for k=0, Nk-1 do begin	
	errlinebackarr[k] = vary0[k] + (lineband[k]*sqrt(varm[k]))^2	; error elements
	endfor	
	
	errlineback = sqrt(total(errlinebackarr));*(lineband[Nk-1]-lineband[0]) 	;add in quadrature
	errtotallineflux = double(sqrt(dcorrpoisfluxerr^2 + errlineback^2))	;adding in quadrature the two random errors: error in background estimation and counting error on the line
	
	;totallineflux = double(total(dcorrline)*(lineband[Nk-1]-lineband[0])-lineback)	;total flux in line in E/s/m^2
	;errline = double(sqrt(errlineback^2+errlinects^2))
	;print, totallineflux, errtotallineflux

	Dcnorm = qromb('Ez', 0.0, z)
	Dc = 1.323e28*Dcnorm	;assuming h = 0.7, Dc is in [cm]
   	lumdist = double((1+z)*Dc)
	linelum = double(totallineflux*4*!pi*lumdist^2)
	errlinelum = double(errtotallineflux*4*!pi*lumdist^2)
	normlinelum = double(linelum)/(10^42D)
	errnormlinelum = double(errlinelum)/(10^42D)
	;errloglinelumneg = double(loglinelum - alog10(linelum-errlinelum))
	printf, lun, mask, '   ', slit, '   ', totallineflux, '   ', errtotallineflux, '   ', linelum, '   ', errlinelum, '   ', normlinelum, '   ', errnormlinelum
; plotting stuff

	loadcolors
	splot, spec.lambda, dcorrspec, xrange=[blueside[0]-20, redside[1]+20], yrange=[-1e-27, 1e-27],psym=10
	oplot, [blueband[0],redband[Nj-1]], [fullcnt[0], fullcnt[Ny-1]], color=2, thick=2
	oplot, [blueband[0],blueband[0]], [-5000,5000], color=4, thick=2
	oplot, [blueband[Ni-1],blueband[Ni-1]], [-5000,5000], color=4, thick=2
	oplot, [redband[0],redband[0]], [-5000,5000], color=4, thick=2
	oplot, [redband[Nj-1],redband[Nj-1]], [-5000,5000], color=4, thick=2
	oplot, [lineband[0],lineband[0]], [-5000,5000], color=1, thick=2
	oplot, [lineband[Nk-1],lineband[Nk-1]], [-5000,5000], color=1, thick=2
	;oplot, [blueband[0],redband[Nj-1]], [fakecnt[0], fakecnt[Ny-1]], color=5
	endfor
	free_lun, lun

	print, median(dcorrspec)
end
	
