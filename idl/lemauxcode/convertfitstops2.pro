pro convertfitstops, name, spectrum, spectrum2, ps
	set_plot, 'PS'
	spec = mrdfits(spectrum,1)
	spec2 = mrdfits(spectrum2,1)
	device, /color, filename= ps
	loadcolors
	xtitle = textoidl('\lambda (Rest Frame)')
	ytitle= textoidl('Normalized Flux (Arbitrary Units)')
	x= dblarr(2)+1.
	plot, spec.lambda, spec.spec, position=[0.12, 0.1, 0.75, 0.94], /nodata, xrange=[1.1500,1.3050], xstyle=1, yrange=[-50,200], ystyle = 1, linestyle=5, color=0, title=title, xtitle= xtitle, ytitle = ytitle
	oplot, spec.lambda, spec.spec, color=0, linestyle=0
	oplot, spec2.lambda, spec2.spec*median(spec.spec)/median(spec2.spec), color=1, linestyle=2;, linestyle=0
	print, median(spec.spec)/median(spec2.spec)
	y = [100,100]
	y2 = [110,110]
	x = [1.16,1.175]
	oplot, x, y, color=0, linestyle=0
	oplot, x, y2, color=1, linestyle=2	
	xyouts, 1.18, 96, charsize=1.3, charthick=6, 'IDL', color=0
	xyouts, 1.18, 106, charsize=1.3, charthick=6, 'IRAF', color=1
	xyouts, 1.215, 180, charsize=1.5, charthick=7, name, color=0
	device, /close_file
end
