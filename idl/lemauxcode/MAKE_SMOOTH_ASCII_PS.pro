pro make_smooth_ascii_ps, sigma
		
	readcol, '16XR3.010.LFC_SC1_01402.new.dat', lambda, spec
	lambda = lambda/1.92706000

	halfwidth = floor(1.5*sigma)
        kernel=findgen(2*halfwidth+1)-halfwidth
        kernel=exp(-kernel^2/2/sigma^2)
        kernel=kernel/total(kernel)
        stemp_flux =  convol(spec,  kernel,  /center)
	
	result = {lambda:lambda, spec: stemp_flux, ivar : 1/abs(stemp_flux)}

	mwrfits, result, '16XR3.010.LFC_SC1_01402.new.smoothed.fits'

	;set_plot, 'PS'
        ;device, /color, xsize = 7, ysize = 8.5, /inches, filename='16XR3.010.spec1d.ps'
        ;loadcolors
        ;xtitle = textoidl('\lambda (Rest Frame)')
        ;ytitle= textoidl('Normalized Flux (Arbitrary Units)')
        ;xtitle = textoidl('\lambda (Rest Frame)')
        ;plot, result.lambda, result.spec, position=[0.12, 0.1, 0.94, 0.6], xrange=[3300,4800], xstyle=1, yrange=[-200,1400], ystyle = 1, $
        ;color=0, title=title, xtitle= xtitle, ytitle = ytitle	
end
