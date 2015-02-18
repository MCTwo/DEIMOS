; -> Converts e- counts to corrected e- counts. <-
;       The inverse of this function returns the throughput of DEIMOS
;	
; 	for use with data with blaze ~ 7000 Ang and GG455 blocking filter
function deimos_correction, $
        lambda ;        Correction factor to DEIMOS e- counts in Angstroms
        return, 1./(-44.248086+0.024340674*lambda-5.0093107e-06*lambda^2+4.6010414e-10*lambda^3-1.5920353e-14*lambda^4)
end

