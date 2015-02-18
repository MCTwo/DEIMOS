;+
; NAME:
;   r2z
;   z2r
;   t2z
;   z2t
;   cosmo_init :initiation routine, needed only once/cosmological model
;
; PURPOSE:
;   functions to convert between redshift z and comoving distance r
;   and to setup this relationship. Information maintained in
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
; KEYWORDS:
;   lookback -- if set, the routine returns the lookback time to
;               redshift z (in units of 1/H0) instead of the
;               distance. 
;
; COMMENTS:
;   to setup a cosmology, execute 'cosmo_init, Omega,w' first. Routines
;   call spline function to evaluate quickly.   Spline has integrated
;   the function 1/E(z) defined by Peebles, 3rd book.
;
; REVISION HISTORY:
;   md 7nov02
;-----------------------------------------------------------------
; convert from z to r1, given spline and a cosmological model
; computes r1(z)
;
function z2r, zz, lookback = lookback
common cosmo_params, Omegam, w_eq_state, zg, ag, r1, lb, t1, $
  splinerz, splinezr, splinelb, splinezt,  splinetz

if max(zz) gt 10 then message, 'z2r only works for z<10'

if n_elements(lookback) eq 0 then lookback = 0

if lookback then return, spl_interp(zg,lb,splinelb,zz) $
  else return, spl_interp(zg,r1,splinerz,zz)
end








