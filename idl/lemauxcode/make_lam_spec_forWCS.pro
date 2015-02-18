pro make_lam_spec_forWCS

	readcol, 'NIR_spec_jun07.cat', format = 'A,A,A', n,m, file

	for i=0, n_elements(file)-1 do begin

		m = mrdfits(file[i],1)
		print, file[i]
		mwrfits, m.lambda, 'lambda.' + file[i]
		mwrfits, m.spec, 'spec.' + file[i]
	endfor

end
