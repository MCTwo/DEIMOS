; -> Converts e- counts to corrected e- counts. <-
;       The inverse of this function returns the throughput of DEIMOS
;	
; 	for use with data with blaze ~ 6000-7300 Ang and GG455 blocking filter (valid for wavelength coverage 4500->8800
function deimos_correction_masterblue, $
        lambda ;        Correction factor to DEIMOS e- counts in Angstroms
        return, 1./(16.0537-0.0112664*lambda+2.98523e-06*lambda^2-3.68807e-10*lambda^3+2.11040e-14*lambda^4-4.37880e-19*lambda^5)
end

