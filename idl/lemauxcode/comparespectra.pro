pro comparespectra

a = mrdfits('LFC_SC1_08407.IDL.griddedtostandard.fits',1)
b = mrdfits('standard.fromIDL_microns.fits',1)
c = mrdfits('LFC_SC1_08407.IDL+IRAF.fits',1)
d = mrdfits('standard.fromAPALL.fits',1)

	bnorm = b.spec/mean(b.spec)
	dnorm = d.spec/mean(d.spec)
	IDLstd = fltarr(n_elements(a.spec))   	
	IRAFstd = fltarr(n_elements(c.spec))

	f = where(b.spec eq 0., count)
	b.spec[f] = b.spec[f] + 0.00001	

	for i=0, n_elements(a.spec)-1 do begin
		IDLstd[i] = a.spec[i]*10000/b.spec[i]
	endfor
	
	for i=0, n_elements(c.spec)-1 do begin
		IRAFstd[i] = c.spec[i]*10000/d.spec[i]	
	endfor

	loadcolors
	device, decomposed=0
	entry_device = !d.name
        set_plot, 'PS'
        device, /color, bits_per_pixel=24, filename='standardized.LFC_SC1_08407.ps', xoffset=0.2, yoffset=0.2, xsize=10, ysize=7, /inches, /landscape
	plot, a.lambda, IDLstd, /nodata, color=0, xstyle=1,xrange=[1.16,1.35], yrange = [-50,130], ystyle = 1, xthick=6, xtitle=textoidl('\lambda_{obs} (\mum)'), ytitle='Flux (Arbitrary Units)', title='LFC_SC1_08407', charsize=2, charthick=6
        oplot, a.lambda, IDLstd, color=1, linestyle=0, thick=6
        oplot, c.lambda, IRAFstd, color=6, linestyle=1, thick=6
       
	x = [1.26,1.29]
	y1 = [105,105]
	y2 = [120,120]
	oplot, x, y1, color=1, linestyle=0, thick=6
	oplot, x, y2, color=6, linestyle=1, thick=6
        xyouts, 1.293, 102, 'IDL only', charthick=5, charsize=2, color=1
	xyouts, 1.293, 117, 'IDL+IRAF', charthick=5, charsize=2, color=6
	device, /close_file

	set_plot, entry_device	

end
