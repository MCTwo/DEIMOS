;+
; NAME:
;     PICKMIN
;
; PURPOSE:
;     Moves a sampling window through an
;     array and picks out the index of the
;     local minimum (or maximum) in each 
;     window. 
;
; CALLING SEQUENCE:
;
;     RESULT = PICKMIN( ARRAY [, NWIN, WINDOW = WINDOW, /MAX])
;
; INPUTS:
;
;     ARRAY:  Input array
;
; KEYWORD PARAMETERS:
;
;     WINDOW: Window size. Default is 32
;     /MAX:   Find maxima instead of minima
;
; OUTPUTS:
;
;     RESULT: Indeces of local minimum within
;             each window.
;
; OPTIONAL OUTPUTS:
;
;     NWIN:   The number of windows used based
;             on size of WINDOW.
;
; Written by JohnJohn in December 2002.
; 09.11.2003   JJ  Implemented strategic use of double precision
;-

function pickmin,arr,nwin,window=window,max = max

if n_elements(window) eq 0 then window = 32
lo = 0d
hi = double(window-1)
nwin = long(n_elements(arr)/window)
mins = fltarr(nwin)
for i = 0, nwin-1 do begin
    range = fillarr(1,lo,hi)
    segment = arr[range]
    if not keyword_set(max) then dummy = min(segment,ind) else $
      dummy = max(segment,ind)
    mins[i] = ind + double(i)*window
    lo = lo+window
    hi = hi+window
endfor
return,mins
end
