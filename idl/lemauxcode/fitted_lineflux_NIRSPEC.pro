;+------------------------------------------------------------
; 
; FITTED_LINEFLUX_NIRSPEC	12/2008
; 
; Program for general flux measurements of NIRSPEC spectra using
; a Gaussian fit to 1 or 2 spectral lines, calculates background 
; using either a constant or linear model. The spectrum is assumed
; to be in units of Jansky's and to be readable by mrdfits on the 
; first extension of the fits table. Will not accept input from
; WCS spectra or multi-extension fits files
;
; INPUTS
; spec - the name of the input spectrum
; lam1 - guess at the central wavelength of the 1st spectral feature
; sig1 - guess at the standard deviation of the 1st spectral feature
; amp1 - guess at the amplitude of the 1st spectral feature
; order - polynomial order of the background fit, can be 0 or 1
;
;
; OPTIONAL KEYWORDS
;
; All optional keywords apply to a 2nd spectral feature, if fitting for
; one feature only ignore these
;
; lam2 - guess at the central wavelength of the 2nd spectral feature 
; sig2 - guess at the standard deviation of the 2nd spectral feature
; amp2 - guess at the amplitude of the 2nd spectral feature
;
;
; OUTPUTS
; As of now, no tangible outputs, only prints out in the IDL command window
; the flux and error of the measured line(s)
;
;
; HISTORY
; written BCL 12/12/08

pro fitted_lineflux_NIRSPEC, spec, lam1, lam2=lam2, sig1, sig2=sig2, amp1, amp2=amp2, order

	if order ne 0 and order ne 1 then begin 
		message, 'Continuum can only be fit to 0th or 1st order polynomial, change order'
	endif

	loadcolors	
	m = mrdfits(spec,1)
	
	range = where(m.lambda gt lam1-250 and m.lambda lt lam1+250)
	
	linespec = m.spec[range]*1e-23                  ;put in ers/s/cm2/Hz, assumes the input spectrum is in Janskys
        linelam = m.lambda[range]
        lineivar = 1/(1/sqrt(m.ivar)*1e-23)^2           ;same

	
	if n_elements(lam2) and n_elements(sig2) and n_elements(amp2) then begin
	
	if order eq 0 then expr = 'P[0] + GAUSS1(X, P[1:3]) + GAUSS1(X,P[4:6])'$
		else expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4]) + GAUSS1(X,P[5:7])'
	
	area1 = sig1*amp1*sqrt(2*!PI)
	area2 = sig2*amp2*sqrt(2*!PI)
	
	if order eq 0 then start = [5e-05, lam1, sig1, area1, lam2, sig2, area2]$
		else start = [5e-05, 5e-05, lam1, sig1, area1, lam2, sig2, area2]
	
	result = MPFITEXPR(expr,m.lambda[range], m.spec[range], sqrt(1/m.ivar[range]), start, perr=perr)
	
	;note area output from MPFITEXPR is already background subtracted so no need to do anything else, just use the outputed error
	
	if order eq 0 then Haflux = result[3]*1e-23*3e18/result[1]^2$	;area is output in units of Jy*Angstrom, convert first to ergs*A/s/cm2/Hz then c/lambda_c^2
		else Haflux = result[4]*1e-23*3e18/result[2]^2 		;this isn't strictly correct, really should do this at every pixel at every wavelength, but close enough

	if order eq 0 then Hafluxerr = perr[3]*3e-05/lam1^2$
		else Hafluxerr = perr[4]*3e-05/lam1^2
	
	if order eq 0 then NIIflux = result[6]*1e-23*3e18/result[4]^2$
		else NIIflux = result[7]*1e-23*3e18/result[5]^2

        if order eq 0 then NIIfluxerr = perr[6]*3e-05/result[4]^2$
                else NIIfluxerr = perr[7]*3e-05/result[5]^2

	splot, m.lambda[range], m.spec[range], psym=10, yrange=[min(m.spec[range])-mean(m.spec[range]), max(m.spec[range])+mean(m.spec[range])]
	
	if order eq 0 then oplot, m.lambda, result[0] + gauss1(m.lambda, result[1:3]) + gauss1(m.lambda, result[4:6]), color=1, thick=5$
		else  oplot, m.lambda, result[0] + result[1]*m.lambda + gauss1(m.lambda, result[2:4]) + gauss1(m.lambda, result[5:7]), color=1, thick=3

	if order eq 0 then reportlam1 = result[1]$
		else reportlam1 = result[2]
		
	if order eq 0 then reportlam2 = result[4]$
		else reportlam2 = result[5]
		
	print, 'The flux of the Gaussian centered at ' + strcompress(string(reportlam1),/REMOVE_ALL) + ' A is ' + strcompress(string(Haflux),/remove_all) + ' +/- ' + strcompress(string(Hafluxerr),/remove_all)
	print, 'The flux of the Gaussian centered at ' + strcompress(string(reportlam2),/REMOVE_ALL) + ' A is ' + strcompress(string(NIIflux),/remove_all) + ' +/- ' + strcompress(string(NIIfluxerr),/remove_all)
	
	endif
	
	if n_elements(lam2) lt 1 then begin

	if order eq 0 then expr = 'P[0] + GAUSS1(X, P[1:3])'$
		else expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4])'

	area1 = sig1*amp1*sqrt(2*!PI)
	
	if order eq 0 then start = [5e-05, lam1, sig1, area1]$
		else start = [5e-05, 5e-05, lam1, sig1, area1]

	result = MPFITEXPR(expr,m.lambda[range], m.spec[range], sqrt(1/m.ivar[range]), start, perr=perr)

	if order eq 0 then Haflux = result[3]*1e-23*3e18/result[1]^2$
		else Haflux = result[4]*1e-23*3e18/result[2]^2

	if order eq 0 then Hafluxerr = perr[3]*3e-05/result[1]^2$
		else Hafluxerr = perr[4]*3e-05/result[2]^2

	splot, m.lambda[range], m.spec[range], psym=10, yrange=[min(m.spec[range])-mean(m.spec[range]), max(m.spec[range])+mean(m.spec[range])]

	if order eq 0 then oplot, m.lambda, result[0] + gauss1(m.lambda, result[1:3]), color=1, thick=3$
		else  oplot, m.lambda, result[0] + result[1]*m.lambda + gauss1(m.lambda, result[2:4]), color=1, thick=3	

	if order eq 0 then reportlam1 = result[1]$
                else reportlam1 = result[2]

	print, 'The flux of the Gaussian centered at ' + strcompress(string(reportlam1),/REMOVE_ALL) + ' A is ' + strcompress(string(Haflux),/remove_all) + ' +/- ' + strcompress(string(Hafluxerr),/remove_all)
	
	endif

end
