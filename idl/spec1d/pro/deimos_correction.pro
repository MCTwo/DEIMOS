; -> Converts e- counts to e- counts. <-
;	The inverse of this function returns the throughput of DEIMOS
;
; npk
function deimos_correction, $	
	lambda ;	Correction factor to DEIMOS e- counts in Angstroms
	return, 1./(-77.9026+0.0395916*lambda-7.49911e-06*lambda^2+6.29692e-10*lambda^3-1.97967e-14*lambda^4)
end
