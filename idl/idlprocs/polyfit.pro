
;+
; NAME: 
;       POLYFIT
;
; PURPOSE:
;	Fit a polynomial to a function using linear least-squares
; NOTE:
;       2 times faster than POLY_FIT.PRO as tested by Robishaw's 
;       BENCHMARK.PRO
;
; CALLING SEQUENCE:
;       coeffs = polyfit(t, y_of_t, degree [,covariance_matrix, yfit])
;-
function polyfit,t,y,deg,cov,yfit
on_error,2
   n = n_elements(t)
   x = dblarr(n)+1
   pow = indgen(deg+1)
   powarr = fan(pow,n,/trans)
   x =  fan(double(t),deg+1)
   xarr = x^powarr
   xarrt = transpose(xarr)
   alpha = xarr##xarrt
   beta = xarr##double(y)
   cov = invert(alpha)
   a = cov##beta
   if n_params() eq 5 then yfit = poly(t,a)
   return,a
end
