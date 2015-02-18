;+
; NAME:
;       RC
;
; PURPOSE:
;       Given a 1D vector of indecies for a 2D array, return the corresponding 
;       column and row indecies. 
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;
;      rc, ind, ncol, r, c
;
; INPUTS:
;
;      IND  - The 1D index
;      NCOL - The number of columns in the 2D array
;
; OUTPUTS:
;
;      R - The row
;      C - The column
;
; EXAMPLE:
;
; IDL> im = findgen(20,200)
; IDL> w = where(im eq 57)
; IDL> rc,w,20,r,c
; IDL> print,r,c
;       2
;      17
; IDL> print,im[c,r]
;      57.0000
;
; MODIFICATION HISTORY:
; Written May 3, 2003 by JohnJohn
; 05.05.2003    Changed fix(ind) to long(ind) to handle looonnnggg arrays.
;-
pro rc,ind,ncol,r,c
ind = long(ind)
ncol = long(ncol)
c = ind mod ncol
r = ind/ncol
end
