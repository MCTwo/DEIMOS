;+
; NAME:
;       VIEW
;
;
; PURPOSE:
;       Look at an array one row at a time. Useful for viewing
;       echelle spectra arrays one order at a time. Click the
;       screen to move to the next row.
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;       VIEW, ARRAY [,SM=SM]
;
;
; INPUTS:
;       ARRAY - 2D array of stuff you wanna look at
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;       SM - Smooth each row by this amount.
;
;       All other keyword parameters are passed on to PLOT. See
;       ? PLOT for list of plotting keywords.
;
; OUTPUTS:
;       Pretty plot.
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
; Written a long time ago by JohnJohn
;-

pro view,wavin,specin,sm=sm,_extra=e
if n_params() eq 1 then begin
    wav = fan(indgen((size(wavin))[1]), (size(wavin))[2])
    spec = wavin 
endif else begin
    wav = wavin
    spec = specin
endelse
nord = (size(spec))[2]
i = 0
x = !x.crange[0]+1
print,'click on plot for next row'
print,'click to the left of the y-axis to stop'
while i ne nord and x gt !x.crange[0] do begin
    s = spec[*,i]
    w = wav[*,i]
    if n_elements(sm) then s = smooth(s,sm)
    plot,w,s,/xs,title='Row '+strtrim(i,2),_extra=e
    if keyword_set(plot) then psclose
    cursor,x,y,/up
    i = i + 1
endwhile
end

