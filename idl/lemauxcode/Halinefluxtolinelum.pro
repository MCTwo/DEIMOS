function Ez, x
    distance = 1/sqrt(0.3*(1+x)^3+0.7) ;concordance cosmology assumtion omega_M = 0.3, omega_lam = 0.7, to be consistant with DEEP2
   return, distance

end

pro Halinefluxtolinelum, cat, outcat

	readcol, cat, format='A,A,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', masklist, slitlist, zlist, EWHa, EWHaerr, EWOII, EWOIIerr, FHa, FHaerr, FNII,FNIIerr, rband, iband, zband, f606w, f814w, RA,Dec, skipline=0
		
	openw, lun, outcat, width=300, /get_lun
        tab = '09'XB

	for n=0, n_elements(masklist)-1 do begin
	
	mask = masklist[n]
	slit = slitlist[n]
	z= zlist[n]
	
	Dcnorm = qromb('Ez', 0.0, z)
        Dc = 1.323e28*Dcnorm    ;assuming h = 0.7, Dc is in [cm]
	lumdist = double((1+z)*Dc)
	Halinelum = double(FHa[n]*1e-18*4*!pi*lumdist^2)
	Halinelumerr = double(FHaerr[n]*1e-18*4*!pi*lumdist^2)
	NIIlinelum = double(FNII[n]*1e-18*4*!pi*lumdist^2)
	NIIlinelumerr = double(FNIIerr[n]*1e-18*4*!pi*lumdist^2)
	Hanormlinelum = double(Halinelum)/(10^42D)
        Hanormlinelumerr = double(Halinelumerr)/(10^42D)
	NIInormlinelum = double(NIIlinelum)/(10^42D)
        NIInormlinelumerr = double(NIIlinelumerr)/(10^42D)

	printf, lun, mask, '   ', slit, '   ', Halinelum, '   ', Halinelumerr, '   ', NIIlinelum, '   ', NIIlinelumerr, '   ', Hanormlinelum, '   ', Hanormlinelumerr, '   ', NIInormlinelum, '   ', NIInormlinelumerr


	endfor
	free_lun, lun

end
