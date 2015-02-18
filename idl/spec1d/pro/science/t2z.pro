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
; convert from t to z, given spline and a cosmological model
; computes z(t)
;
function t2z, tt
common cosmo_params, Omegam, w_eq_state, zg, ag, r1, lb, t1, $
  splinerz, splinezr, splinelb, splinezt,  splinetz

splineoutput = spl_interp(t1, ag, splinezt, tt) ; scale factor a

if max(splineoutput) gt 10 or min(splineoutput) lt -0.1 then message, 'ERROR: Z<10 is an inaccurate extrapolation in t2z'

return, (1./splineoutput)-1.
end
