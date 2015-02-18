pro individual_coadd_script, cat
	numread = 205
	for i=0, numread-1 do begin
	result = coadd_cluster_spectra(cat, i,normalize=normalize,zero=zero,inorm=inorm,weight=weight)
	endfor
end
