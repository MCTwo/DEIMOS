pro convertfitstops, ps
	set_plot, 'PS'
	device, /color, xoffset=0.2, yoffset=0.2, xsize=8, ysize=6.5, /inches, /portrait, filename= ps
	loadcolors
	title = 'Cl0023 vs. N5281 All Priorities 1st 2 masks'
	xtitle = textoidl('\lambda_{rest} (\AA)')
	ytitle= textoidl('f_{\lambda} (Arbitrary Units)')
	x= dblarr(2)+1.
	spec = mrdfits('cl0023.1st2masks.allpri.fits',1)
	spec2 = mrdfits('nep5281.1st2masks.allpri.fits',1)
	ratio = spec.spec/spec2.spec
	scaling = median(ratio)
	print, scaling
	spec2.spec = spec2.spec+0.08
	plot, spec.lambda, spec.spec, /nodata, position=[0.12, 0.15, 0.93, 0.92], xrange=[3600,4750], xstyle=1, yrange=[0,1.95], ystyle = 1, linestyle=5,$
	color=0, title=title, xtitle= xtitle, ytitle = ytitle, charsize = 1.5, ythick=6, xthick=6, charthick=5
	oplot, spec.lambda, spec.spec, color=0, linestyle=0, thick=2 
	oplot, spec2.lambda, spec2.spec, color=1, linestyle=2, thick=2;, linestyle=0
	linex = dblarr(2) + 3650
	linex[1] = linex[1]+100
	liney = dblarr(2) + 1.75
	oplot, linex, liney, linestyle=0, thick=2, color=0
	xyouts, 3775, 1.73, 'CL0023', color=0, charsize=1.6, charthick=5
	liney = liney - 0.11
	oplot, linex, liney, linestyle=2, thick=2, color=1
	xyouts, 3775, 1.62, 'N5281', color=0, charsize=1.6, charthick=5

	plot, spec.lambda, spec.spec, position=[0.53,0.15,0.73,0.45], xrange=[3715,3740], xstyle=4, yrange=[0.3,1.5], ystyle=4, linestyle=0, color=0, xtitle =' ', ytitle = ' ', charsize=1.5, charthick=5, /noerase
	axis, xaxis=1, xticks=2, xtickv=[3723,3733], charthick=5, charsize=1.5
	axis, yaxis=0, yticks=3, ytickv=[0.5,1, 1.5], charthick=5, charsize=1.5
	;axis, xaxis=1, xticks=2, xtickv=[3725,3733], xtickname=replicate(' ', 2)
        ;axis, yaxis=0, yticks=3, ytickv=[0.5,1, 1.5], ytickname=replicate(' ', 3)
	oplot, spec2.lambda, spec2.spec, color=1, linestyle=2	
	xyouts, 3717, 0.35, '[OII]', color=0, charsize=1.5, charthick=5
		
	plot, spec.lambda, spec.spec, position=[0.73,0.15,0.93,0.45], xrange=[4080,4120], xstyle=4, yrange=[0.7,1.4], ystyle=4, linestyle=0, color=0, xtitle =' ', ytitle = ' ', charsize=1.5, charthick=5, /noerase
	axis, xaxis=1, xticks=2, xtickv=[4092,4109], charthick=5, charsize=1.5
        axis, yaxis=1, yticks=3, ytickv=[0.75,1, 1.3], charthick=5, charsize=1.5
        axis, yaxis=0, yticks=3, ytickv=[0.8,1, 1.2], ytickname=replicate(' ', 3)
	oplot, spec2.lambda, spec2.spec+0.08, color=1, linestyle=2	
	xyouts, 4083, 0.7, 'H' + textoidl('\delta'), color=0, charsize=1.5, charthick=5
	
	device, /close_file
	
end
