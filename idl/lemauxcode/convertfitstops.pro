pro convertfitstops, name, spec, ps
	set_plot, 'PS'
	spectra = mrdfits(spec,1)
	device, /color, filename= ps
	loadcolors
	xtitle = textoidl('\lambda (Rest Frame)')
	ytitle= textoidl('F_{\nu} (\muJy)')
	plot, spectra.lambda, spectra.spec*1e6, position=[0.12, 0.1, 0.75, 0.94], xrange=[6200,7050], yrange=[-10,50], xstyle = 1, ystyle = 1, color=0, title=title, xtitle= xtitle, ytitle = ytitle, linestyle=0
	;oplot, spectra2.lambda, spectra2.spec*1000, color=1, linestyle=2
	xyouts, 6700, 42, name, charthick=6, charsize = 1.6
	device, /close_file
	
end
