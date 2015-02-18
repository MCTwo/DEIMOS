;+
; NAME:
;   r2z
;   z2r
;   t2z
;   z2t
;   cosmo_init :initiation routine, needed only once/cosmological model
;
; PURPOSE:
;   functions to convert between redshift z and comoving distance r or
;   cosmic time t,  
;   and to setup these relationships. Information maintained in
;   cosmo_params common block .
;
; CALLING SEQUENCE:
;   cosmo_init, Omega,w  
;   z = r2z(r)
;   r = z2r(z)
;   t = z2t(z)
;   z = t2z(t)
; 
; INPUTS:
;   r -- comoving distance (in units of the horizon distance, c/H0)
;   or z -- observed redshift
;   or t -- cosmic time (in units of 1/H0)
;   Omega, w -- parameters for cosmological model (flat)
;
; COMMENTS:
;   to setup a cosmology, execute 'cosmo_init, Omega,w' first. Routines
;   call spline function to evaluate quickly.   Spline has integrated
;   the function 1/E(z) defined by Peebles, 3rd book.
;
; REVISION HISTORY:
;   based on z2r.pro: md 7nov02
;   altered by BFG 8jul04
;-----------------------------------------------------------------
; convert from z to t, given spline and a cosmological model
; computes t(z)
;
function z2t, zz
common cosmo_params, Omegam, w_eq_state, zg, ag, r1, lb, t1, $
  splinerz, splinezr, splinelb, splinezt,  splinetz

;if max(zz) gt 10 then message, 'z2t only works for z<10'
if min(zz) lt 0 then message, 'z2t only works for positive z!'
return, spl_interp(ag,t1,splinetz,1./(1.+zz)) 

end



