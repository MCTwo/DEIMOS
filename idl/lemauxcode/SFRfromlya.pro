pro SFRfromlya, outcat
	
	readcol, '1604LAE_Flux_and_Lum_newcalc_FINAL_inclslitloss.dat' , format = 'A,A,D,D,D,D,D,D,D', mask,slit,flux,fluxnegerr,fluxposerr, lum, lumnegerr, lumposerr, normlum, normlumnerr, normlumperr
	readcol, '1604LAE_EWs_inclslitloss_watten.dat', format='A,A', slit, mask, z, atten

	SFR = fltarr(n_elements(lum))
	SFRnegerr = fltarr(n_elements(lum))
	SFRposerr = fltarr(n_elements(lum))	
	openw, lun, outcat, /get_lun
	printf, lun, 'Mask  Slit   z  Corrected L(Lya)        SFR[M_solar/yr]   -err   +err'

	for i=0, n_elements(lum)-1 do begin
	SFR[i] = lum[i]*9.52e-43/atten[i]
	SFRnegerr[i] = lumnegerr[i]*9.52e-43/atten[i]
	SFRposerr[i] = lumposerr[i]*9.52e-43/atten[i]	
	printf, lun, mask[i], '  ', slit[i], '  ', z[i], '  ', lum[i]/atten[i], '       ', SFR[i], '  ', SFRnegerr[i], '  ', SFRposerr[i] 
	endfor

	free_lun, lun

end
