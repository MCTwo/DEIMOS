function jjgaussfit_func, p, x = x, y = y, err = err
z = (x - p[1])/p[2]
fit = p[0]*exp(-z^2) + p[3]
if n_elements(err) eq 0 then err = 1.
return, (y - fit)/err
end

function jjgaussfit, x, y, p $
                     , yerr = yerr $
                     , fixsig=fixsig, limitsig=limitsig $
                     , fixcen=fixcen , limitcen=limitcen $
                     , fixamp=fixamp , limitamp=limitamp $
                     , fixbg = fixbg, limitbg=limitbg $
                     , loud=loud
;+
; NAME:
;          JJGAUSSFIT
;
;
; PURPOSE:
;         Fits a gaussian to a 1D distribution of points, This is a wrapper
;         for Craig Markwardt's infinitely useful fitting package MPFIT.
;
;         Unlike IDL's GAUSSFIT, this function A) optionally accounts
;         for measurement errors and B) Works more than 80% of the
;         time without crashing.
;
;         Fits the function:
;
;         f(x) = a0 * exp(-z^2) + a3
;
;         where
;
;         z = (x - a1)/a2
;
; CATEGORY:
;
;         Getting the job done (TCB).
;
; CALLING SEQUENCE:
;
;         fit = jjgaussfit(x, y [, param, yerr=])
;
; INPUTS:
;         x:  independant variable
;         y:  independant variable
;
; OPTIONAL INPUTS:
;
;         yerr: 1 sigma errors associated with each Y
;
; OUTPUTS:
;
;         fit: The best-fit Gaussian evaluated at each X
;         param: 4 element vector containing the fit parameters
;                param = [a0, a1, a2, a3]
;
; MODIFICATION HISTORY:
; 12 Dec 2003 Written by JohnJohn
;-

if keyword_set(yerr) then begin
    fa = {x:double(x), y:double(y), err:double(yerr)} 
endif else begin
    fa = {x:double(x), y:double(y)}
endelse
;;;Set up parameter guesses
mxy = max(y, imx)
mny = min(y)
test = y - mny
w = where(test gt 0.5*(mxy - mny), fwhm)
;p0 = par guesses
p0 = double([mxy, x[imx], fwhm/2.355, mny]) 
npar = n_elements(p0)
parinfo = replicate( { fixed: 0b, $
                       limited: [0b,0b], $
                       limits: dblarr(2) $
                     }, npar)

if n_elements(fixamp) gt 0 then begin
    parinfo[0].fixed = 1b
    p0[0] = fixamp
endif

if n_elements(fixcen) gt 0 then begin
    parinfo[1].fixed = 1b
    p0[1] = fixcen
endif

if keyword_set(fixsig) then begin
    parinfo[2].fixed = 1b
    p0[2] = fixsig
endif

if n_elements(fixbg) gt 0 then begin
    parinfo[3].fixed = 1b
    p0[3] = fixbg
endif

if keyword_set(limitamp) then begin
    parinfo[0].limited = 1b
    parinfo[0].limits = limitamp
endif

if keyword_set(limitcen) then begin
    parinfo[1].limited = 1b
    parinfo[1].limits = limitcen
endif

if keyword_set(limitsig) then begin
    parinfo[2].limited = 1b
    parinfo[2].limits = limitsig
endif

if keyword_set(limitbg) then begin
    parinfo[3].limited = 1b
    parinfo[3].limits = limitbg
endif

p = mpfit('jjgaussfit_func', p0, parinfo=parinfo, funct=fa, maxiter=10 $
          , quiet=1-keyword_set(loud))

z = (x - p[1])/p[2]
fit = p[0]*exp(-z^2) + p[3]
return, fit
end
