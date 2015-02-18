pro continuum_break, spec
	
	spec = mrdfits(spec,1)
	
	bluer = where(spec.lambda ge 1164.7 and spec.lambda lt 1189.7)
	blue = where(spec.lambda ge 1189.7 and spec.lambda lt 1214.7)
	red = where(spec.lambda ge 1218.7 and spec.lambda lt 1243.7)
	redder = where(spec.lambda ge 1243.7 and spec.lambda lt 1268.7)

	bluerflux = mean(spec.spec[bluer])
	blueflux = mean(spec.spec[blue])
	redflux = mean(spec.spec[red])
	redderflux = mean(spec.spec[redder])
	
	bluererr = stddev(spec.spec[bluer])/sqrt(n_elements(spec.spec[bluer]))
	blueerr = stddev(spec.spec[blue])/sqrt(n_elements(spec.spec[blue]))
	rederr = stddev(spec.spec[red])/sqrt(n_elements(spec.spec[red]))
	reddererr = stddev(spec.spec[redder])/sqrt(n_elements(spec.spec[redder]))
	
	bluerlam = (1164.7+1189.7)/2
	bluelam = (1189.7+1214.7)/2
	redlam = (1218.7+1243.7)/2
	redderlam = (1243.7+1268.7)/2
		
	lam = [bluerlam, bluelam, redlam, redderlam]
	flux = [bluerflux, blueflux, redflux, redderflux]
	err = [bluererr, blueerr, rederr, reddererr]
	
	plot, lam, flux
	err_plot, lam, flux-err, flux+err
	print, lam
	print, flux
	print, err

	loadcolors
	set_plot, 'PS'
	device, /color, filename='1604q1.continuumbreak.ps'
	plot, lam, flux, ytitle= textoidl('<F_{\lambda}> (Arbitrary Units)'), xtitle = textoidl('wavelength (\AA)'), title=textoidl('1604 Q=1 Ly\alpha emitting galaxies (4 galaxies)'), /noerase, color=0
	oplot, lam, flux, psym=5
	err_plot, lam, flux-err, flux+err, /noerase
	device, /close_file
end
