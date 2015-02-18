; -> Converts e- counts to e- counts. <-
;	The inverse of this function returns the throughput of DEIMOS
;
; npk
; modified BL 5/17/10, this is the program for spectral setups w/ lambda_c ~ 8000,
; should use redular deimos correction for lambda < 7900. Definitely DO NOT use for
; data that goes blueward of 5600 A

function deimos_correction_reallyred, $	
	lambda ;	Correction factor to DEIMOS e- counts in Angstroms
	return, 173.68298-0.12192296*lambda+3.3735416e-05*lambda^2-4.6030713e-09*lambda^3+3.1030358e-13*lambda^4-8.2818479e-18*lambda^5 
end
