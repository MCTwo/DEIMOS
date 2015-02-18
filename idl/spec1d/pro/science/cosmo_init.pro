;+
; NAME:
;   cosmo_init :initiation routine, needed only once/cosmological model
;
; PURPOSE:
;   functions to convert between redshift z and comoving distance r or
;   cosmic time t
;   and to setup this relationship. Information maintained in
;   cosmo_params common block
;
; CALLING SEQUENCE:
;   r = r2z(z, lookback=lookback)
;   z = z2r(r)
;   t = z2t(z)
;   z = t2z(t)
;   cosmo_init, Omega,w  
; 
; INPUTS:
;   r -- comoving distance
;   or z -- observed redshift
;   or t -- cosmic time
;   Omega, w -- parameters for cosmological model (flat)
;
; COMMENTS:
;   to setup a cosmology, execute 'cosmo_init, Omega,w' first. Routines
;   call spline function to evaluate quickly.   Spline has integrated
;   the function 1/E(z) defined by Peebles, 3rd book.  use r2z r2t,
;   t2z and z2r
;   as evaluation of spline
;
; REVISION HISTORY:
;   md 7nov02
;   bfg 8jul04: added cosmic and lookback time routines.
;---------------------------------------------------------------------
; this simple routine is 1/sqrt(E(z) as defined by Peebles, 3rd book.
; it is intended to be called by qsimp within xisp.pro, to evaluate
; comoving distance r1(z)
; specialized for flat cosmologies specified by Omegam and w
; md 11oct00

   function int_r1, z
   common cosmo_params, Omegam, w, zg, ag, r1, lb, t1, $
   splinerz, splinezr, splinelb, splinezt,  splinetz
   return, 1./sqrt(Omegam*(1+z)^3 + (1.-Omegam)*(1+z)^(3*(1+w)))
   end


   function int_lb, z
   common cosmo_params, Omegam, w, zg, ag, r1, lb, t1, $
   splinerz, splinezr, splinelb, splinezt,  splinetz
   return, 1./sqrt(Omegam*(1+z)^3 + (1.-Omegam)*(1+z)^(3*(1+w)))/(1+z)
   end

   function int_tc, a
   common cosmo_params, Omegam, w, zg, ag, r1, lb, t1, $
   splinerz, splinezr, splinelb, splinezt,  splinetz
   int = a*0.
   wn0 = where(a ne 0)
   if wn0[0] ne -1 then $
     int[wn0] =  a[wn0]/sqrt(Omegam*a[wn0]+(1.-Omegam)*a[wn0]^(1-3*w))
   return, int
   end

;_____________________________________________________
; initiation of spline, called only once
pro cosmo_init, Omega, w
common cosmo_params, Omegam, w_eq_state, zg, ag, r1, lb, t1, $
  splinerz, splinezr, splinelb, splinezt,  splinetz

if n_params() eq 0 then begin ;default cosmological model
  Omega = .3
  w = -1.
endif
if n_params() eq 1 then w = -1.

Omegam = Omega ;transfer variable
w_eq_state = w
zg=findgen(501)/50. ; 0<z<10
ag = findgen(501)/500 ; 0<a<1
r1 = qsimp('int_r1',0,zg) ; r1(z) on a grid
lb = qsimp('int_lb',0,zg) ; t_l(z) on a grid
t1 = qromb('int_tc',0,ag) ; t_c(a) on a grid

splinerz=spl_init(zg,r1)  ; initiate spline to compute r1(z)
splinelb=spl_init(zg,lb)  ; initiate spline to compute t_l(z)

splinezr=spl_init(r1,zg)  ; initiate spline for z(r1)

splinetz = spl_init(ag,t1); initiate spline for t1(z)
splinezt = spl_init(t1,ag); initiate spline for z(t1)

;plot, zg,r1*3000.,xtitle='redshift z', ytitle=' comoving coordinate, Mpc'

return
end 













