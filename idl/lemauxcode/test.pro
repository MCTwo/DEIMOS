pro test, user
	;read in the cluster members to be coadded from a plain text file
	d2 = getenv('D2_RESULTS')
	
	numread = 100

	readcol, d2 + '/cluster_catalog/coadd/spec1d_coadd/sngl_all_cluster_for_coadd.dat',  format='A,I,D,D,A,A', mask, slit, Iband, z, q, file, skipline=2, numline=numread
	
	
	coaddmembers = dblarr(numread, 20)
	for i=0, numread-1 do begin
	mm = mrdfits(d2 + '/finalzspec/' + 'zspec.' + user + '.' + mask(i) + '.fits',1)
	endfor	
end
