pro make_herm,n,hout,transpose=transpose
;+
; NAME:
;         MAKE_HERM
;
; PURPOSE:
;         Generate the coefficients for the first 
;         N Hermite Polynomials, H_0...H_n
;
; CALLING SEQUENCE:
;
;         MAKE_HERM, N, C_n, [,/TRANSPOSE]
;
; INPUTS:
;         N:   The desired number of terms.
;
; KEYWORD PARAMETERS:
;
;         /TRANSPOSE  Set if you want the
;                     C_n matrix to be
;                     transposed
;
; OUTPUTS:
;         C_n: An (N+1) x (N+1) matrix of 
;              coeffiecients
;
; EXAMPLE:
;
; The first 3 Hermite Polynomials are
; H_0 = 1
; H_1 = 2x
; H_2 = -2 + 4x^2
;
; IDL> make_herm, 2, c
; IDL> print, c
;       1.0000000       0.0000000       0.0000000
;       0.0000000       2.0000000       0.0000000
;      -2.0000000       0.0000000       4.0000000
;
; MODIFICATION HISTORY:
; Created sometime in 2002 by JohnJohn for PSF 
; modeling. If you use this, let me know!
;-

h = dblarr(n+1,n+1)
h[0,0] = 1d
if n gt 0 then h[1,1] = 2d

if not keyword_set(transpose) then begin
    if n gt 1 then begin
        for i = 2,n do begin
            h[*,i] = shift(h[*,i-1],1)*2d - 2d*double(i-1)*h[*,i-2]
        endfor
    endif
endif else begin
    if n gt 1 then begin
        for i = 2,n do begin
            h[i,*] = shift(h[i-1,*],1)*2d - 2d*double(i-1)*h[i-2,*]
        endfor
    endif
endelse

hout = h

end
