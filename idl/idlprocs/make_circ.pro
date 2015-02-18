;+
; NAME: 
;          MAKE_CIRC
;
; PURPOSE: 
;          create a circle on a pre-existing array
;
; CATEGORY:
;          
;          array
;
; CALLING SEQUENCE:
;
;          result = make_circ(array, radius [, indices, cent = [xo,yo],
;          val = val])
;
; INPUTS:
;
;          array:    2-D array on which to draw the circle
;          radius:   The radius in pixels of the circle
;
; OPTIONAL INPUTS:
;
;          indices:  array of the index values of ARRAY which fall
;                    within one RADIUS of the center
;
; KEYWORD PARAMETERS:
; 
;          center:   array holding the x and y coordinate of the
;                    centroid, center of the input array chosen by default
;          val:      value to set array elements which fall within one RADIUS
;
; OUTPUTS:
;
;         2-D array
;
; EXAMPLE:
;
;         to draw a circle on a blank array
;         IDL> blank = fltarr(5,5)
;         IDL> result = make_circ(blank, 1)
;         IDL> print, result
;          0.00000      0.00000      0.00000      0.00000      0.00000
;          0.00000      0.00000      1.00000      0.00000      0.00000
;          0.00000      1.00000      1.00000      1.00000      0.00000
;          0.00000      0.00000      1.00000      0.00000      0.00000
;          0.00000      0.00000      0.00000      0.00000      0.00000
;
; MODIFICATION HISTORY:
;         Originally written by JohnJohn sometime around ought-one.
;         Robishaw taught me how to do this. He is wicked cool.
;-

function make_circ, f, r, w,val=val, center=center
on_error,2                      ;Return to sender if broken

  s = size(f)
  arr = f
  if s[0] ne 2 then message,'input must be a 2D array',/ioerror

  ncol = s[1]
  nrow = s[2]

  xgrid = (fltarr(nrow)+1)##indgen(ncol)
  ygrid = indgen(nrow)##(fltarr(ncol)+1)

  if not keyword_set(center) then center = [fix(ncol/2.),fix(nrow/2.)]
  w = where((xgrid-center[0])^2+(ygrid-center[1])^2 le r^2)

  if n_elements(val) eq 0 then val = 1.
  arr[w] = val
return, arr
end

