;+
; NAME:
; dvoldz
;
;
; PURPOSE:
; Computes differential volume element dV/dz as a function of z, given
; certain cosmological parameters.  Relies on cosmo_init.pro and
; z2r.pro.  Assumes a flat universe!
;
;
; CALLING SEQUENCE:
; result = dvoldz, z, OmegaM=omegam, $
;                  wQ=wq,H0=H0, noinit=noinit
;
;
; INPUTS:
; 
;   z         Redshift or array of redshifts
;
; KEYWORD PARAMETERS:
;
;   OmegaM       Matter density parameter
;   wQ           Dark energy equation of state
;   H0           Hubble constant in km/s/Mpc
;   noinit       If set, cosmo_init is not run.  Causes code to run
;                faster if cosmo_init has already been run elsewhere.
;                Causes code to crash if not. 
;
; OUTPUTS:
;
;   result       dV/dz(z)
;
;
; RESTRICTIONS:
;
;   If run with /noinit keyword set, then cosmo_init.pro MUST have
;   been run elsewhere in the code containing the call to dvoldz
;
;
; MODIFICATION HISTORY:
;   BFG 16Mar04
;-

function dvoldz, z, OmegaM=omegam, $
                 wQ=wq,H0=H0,noinit=noinit

if n_elements(omegam) eq 0 then omegam=0.3
if n_elements(wq) eq 0 then wq=-1.
if n_elements(h0) eq 0 then h0=100.
if not keyword_set(noinit) then cosmo_init, omegam, wq

omegaq=1.-omegam
omegar=4.183e-5*(100/h0)^2 ;radiation density, from Peacock p. 664

coh=3e5/h0

rz=coh*z2r(z)

drdz=coh*1./sqrt(omegaq*(1+z)^(3*(1+wq)) + omegam*(1+z)^3 + omegar*(1+z)^4)

;From equn on p. 667 of Peacock (and prior equations)
dvoldz=4*!pi*rz^2*drdz

return, dvoldz

end
