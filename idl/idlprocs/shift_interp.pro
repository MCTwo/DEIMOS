function shift_interp,spec,move,spline=spline
;+
; NAME:
;       SHIFT_INTERP
;
; PURPOSE:
;       Shift an array by a non-integer amount 
;       using interpolation.
;
; CALLING SEQUENCE:
;
;       RESULT = SHIFT_INTERP( SPEC, MOVE, /SPLINE )
;
; INPUTS:
;
;       SPEC:   Array to be shifted
;       MOVE:   Amount to shift
;
; OUTPUTS:
;
;       RESULT:  Shifted array
;
; KEYWORD PARAMETERS:
;
;       SPLINE:  Use spline rather than linear interpolation.
;
; SIDE EFFECTS:
;
;       Don't use this with excessively noisy
;       data lest you interpolate noise.
;
; EXAMPLE:
;  
; Create a gaussian and shift it to the right 
; by 0.5
;
;     IDL> x = findgen(81)-40
;     IDL> g = exp(-x^2/10.^2)
;     IDL> plot,x,g,/xs
;     IDL> gs = shift_interp(g,0.5)
;     IDL> oplot,x,gs,lines=3
;     IDL> gsreal = exp(-(x-0.5)^2/10.^2)
;     IDL> print,stdev(gsreal[1:*] - gs[1:*])/max(gs)
;            1.0892190e-06
;
; MODIFICATION HISTORY:
; 02.08.2003   Written by JohnJohn
; 04.28.2003   JJ - Fixed problem where nothing would be returned for
; integer shifts. Now handles integer shifts as if it was SHIFT.PRO
; 05.02.2003   JJ - Now uses faster spline method. SPLINE.PRO is
; slooowwwww...
; 06.04.2003   JJ - Now uses linear interpolation by default and spline
; interpolation as a keyword option.
; 06.09.2003   JJ - Corrected mistake by repacing if bigmove gt 0 with
;                   if abs(bigmove) gt 0
;-
on_error,2 ;If broke, return to sender.

if move eq 0 then return,spec
fracmove = (move mod 1)
if fracmove ne 0 then begin     
    if keyword_set(spline) then begin ;spline onto fractional scale
        x = indgen(n_elements(spec))
        xnew = x + fracmove
        specnew = spl_interp(xnew,spec,spl_init(xnew,spec),x,/double)
    endif else if fracmove gt 0 then $ ;linear interpoltation
          specnew = (1-fracmove)*spec + fracmove*shift(spec,1) else $
          specnew = (1-abs(fracmove))*spec + abs(fracmove)*shift(spec,-1)
endif else specnew = spec
bigmove = fix(move)  ;do integer part of shift 
if abs(bigmove) gt 0 then specbm = shift(specnew, bigmove) else specbm = specnew

return,specbm
end
