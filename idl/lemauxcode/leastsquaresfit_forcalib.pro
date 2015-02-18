pro leastsquaresfit_forcalib, degree, zmin=zmin, zmax=zmax
		
	if n_elements(zmin) then zmin = zmin[0] else zmin = -50
	if n_elements(zmax) then zmax = zmax[0] else zmax = 50

	;readcol, 'ibandvsspeciband.dat', format = 'A,A,A,D,D,D,D', a,b,c,d,iraw,e,ispecraw, z, q
	readcol, 'uhhhh.cat', format = 'A,A,A,D,D,D,D', a,b,c,d,iraw,e,ispecraw, z, q	

	i = iraw(where(z ge zmin and z le zmax))
	ispec = ispecraw(where(z ge zmin and z le zmax))

	print, n_elements(i)
	
	banddiff = ispec - i
	newiband = i[where(banddiff ge -1.1 and banddiff le 2)]
	newspeciband = ispec[where(banddiff ge -1.1 and banddiff le 2)]
	newbanddiff = banddiff[where(banddiff ge -1.1 and banddiff le 2)]	

	print, n_elements(newiband)

	result = poly_fit(newiband, newspeciband, degree)
	result2 = poly_fit(newiband, newbanddiff, degree)

	meanbanddiff = strcompress(string(mean(newbanddiff)), /remove_all)
	sigbanddiff = strcompress(string(stddev(newbanddiff)), /remove_all)

	x = findgen(500)
	y = result[0] + result[1]*x
	splot, newiband, newspeciband, psym=7
	oplot, x, y, linestyle=2
	y00 = strcompress(string(result[0]), /remove_all)
	m0 = strcompress(string(result[1]), /remove_all)
	y01 = strcompress(string(result2[0]), /remove_all)
	m1 = strcompress(string(result2[1]), /remove_all)

	

	print, 'The best fit for the correlation between i_SDSS and i_spec is i_SDSS =', y00,' +', m0, '*i_spec'
	print, 'The best fit for the correlation between i_SDSS and delta_i is delta_i =', y01,' +', m1, 'i_SDSS'

	print, 'The mean delta_i is ', meanbanddiff, 'and the rms is ', sigbanddiff

end
