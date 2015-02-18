; -> Converts e- counts to normalized throughput corrected e- counts. <-
;	The inverse of this function returns the throughput of our LRIS 400 l/mm data (8500 A blaze, tilted to lambda_c = 7500
;	
function lris_correction, $	
	lambda ;	Correction factor to LRIS e- counts in Angstroms
	return, 1./(-11.705724+ 0.0032102149*lambda-1.3619656e-07*lambda^2-8.3788626e-12*lambda^3)
end
