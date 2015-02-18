;+
; NAME:
;      VLINE
;     
; PURPOSE:
;      Draw a vertical line on a pre-existing plot window.
;
; CALLING SEQUENCE:
;      VLINE, VAL
;
; INPUTS:
;
;      VAL: The x-value or array of x-values where the vertical
;      line(s) should be drawn
;
; KEYWORD PARAMETERS:
;
;      All keyword parameters are passed to OPLOT.
;
; SIDE EFFECTS:
;
;      Causes a vertical line to appear on your screen.
;
; RESTRICTIONS:
;
;      This program won't do anything else. Sorry, them's the 
;      restrictions.
;
; EXAMPLE:
;
;      Draw a vertical line at x = 0
;      IDL> plot, findgen(10)
;      IDL> vline, 5
;
; MODIFICATION HISTORY:
; Written sometime in 2003 by JohnJohn
;-

pro vline, val,_extra=extra
nv = n_elements(val)
for i = 0,nv-1 do oplot,fltarr(2)+val[i],!y.crange,_extra=extra 
end
