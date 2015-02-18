pro lineshapeplot
	set_plot, 'PS'
	spec = mrdfits('smoothed.15sigma.1604q2n3LAE.fits',1)
	spec2 = mrdfits('smoothed.15sigma.1604q1LAE.fits',1)
	spec3 = mrdfits('smoothed.15sigma.1604sngl_nonLYA.fits',1)
	device, xoffset=0.2, yoffset=0.2, xsize=8, ysize=10, /inches, filename= '1604lineshapes.ps'

	plot, spec.lambda, spec.spec, /nodata, position=[0.12, 0.66, 0.9, 0.95], xrange=[1210,1221], xstyle=4, yrange=[-10,490], ystyle = 1,$
	xtitle = textoidl('\lambda_{rest} (\AA)') , ytitle= textoidl('f_{\lambda} (Arbitrary Units)'), $
	thick=7, charthick=4, charsize=1.5, ythick=6, /noerase
	oplot,  spec.lambda, spec.spec, thick=4
	aparams = [13.968, 1215.874, 0.3761, 351.44]
	xyouts, 1210.5, 400, 'Class 2/3 LAE candidates', charthick=5, charsize = 1.25
	xyouts, 1211.5, 350, '(13 Galaxies)', charthick=5, charsize = 1.25
	oplot, spec.lambda, aparams[0] + gauss1(spec.lambda, aparams[1:3]), linestyle=2, thick=4
	axis, xaxis=0, yrange=[1210,1221], xstyle=1, xthick=6, xtickname=replicate(' ', 7)	
	axis, xaxis=1, yrange=[1210,1221], xstyle=1, xthick=6, xtickname=replicate(' ', 7)

	plot, spec2.lambda, spec2.spec, /nodata, position=[0.12, 0.37, 0.9, 0.66], xrange=[1210,1221], xstyle=4, yrange=[-10,490], ystyle = 1,$
	xtitle = textoidl('\lambda_{rest} (\AA)') , ytitle= textoidl('f_{\lambda} (Arbitrary Units)'), $
	thick=7, charthick=4, charsize=1.5, ythick=6, /noerase
	bparams = [21.966,1215.7437, 0.2884,227.05]
	xyouts, 1210.5, 400, 'Class 1 LAE candidates', charthick=5, charsize = 1.25
	xyouts, 1211.5, 350, '(4 Galaxies)', charthick=5, charsize = 1.25
	oplot,  spec2.lambda, spec2.spec, thick=4
	oplot, spec2.lambda, bparams[0] + gauss1(spec2.lambda, bparams[1:3]), linestyle=2, thick=4
	axis, xaxis=0, yrange=[1210,1221], xthick=6, xstyle=1, xtickname=replicate(' ', 7)         
        axis, xaxis=1, yrange=[1210,1221], xthick=6, xstyle=1, xtickname=replicate(' ', 7)
	
	plot, spec3.lambda, spec3.spec, /nodata, position=[0.12, 0.08, 0.9, 0.37], xrange=[1210,1221], xstyle=1, yrange=[-10,499], ystyle = 1,$
        xtitle = textoidl('\lambda_{rest} (\AA)') , ytitle= textoidl('f_{\lambda} (Arbitrary Units)'), $
        thick=7, charthick=4, charsize=1.5, xthick=6, ythick=6,/noerase
	xyouts, 1210.5, 400, 'Single Emission Non-LAE', charthick=5, charsize=1.25
	xyouts, 1211.5, 350, '(22 galaxies)', charthick=5, charsize=1.25
	cparams = [39.51,1215.4159, 0.173, 160.5749]
        oplot,  spec3.lambda, spec3.spec, thick=4
	oplot, spec3.lambda, cparams[0] + gauss1(spec3.lambda, cparams[1:3]), linestyle=2, thick=4
	device, /close_file
end
