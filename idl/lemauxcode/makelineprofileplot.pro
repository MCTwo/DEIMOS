pro makelineprofileplot, bestlam, bestconvol, bestsig, ps

	readcol, 'original_gauss_bestfit.dat', x, bestogauss
	readcol, 'truncated_gauss_bestfit.dat', x, bestgauss
	
	;spectrum = mrdfits('1604q2n3LAE.fits',1)
        ;spectrum = mrdfits('../1604q2n3LAE_unitweight.fits',1)
	spectrum = mrdfits('../smoothed.10sigma.1604q2n3LAE.fits',1)
	range = where(spectrum.lambda gt 1210 and spectrum.lambda lt 1222)      ;limiting the range so things don't get fucked up by edge effects       

        specerr = 1/sqrt(spectrum.ivar[range])
        flux = spectrum.spec[range]
        lambda = spectrum.lambda[range]	
	
	loadcolors 
        entry_device = !d.name
        set_plot, 'PS'
        device, /color, xoffset=0.2, yoffset=0.2, xsize=8, ysize=8, /inches, filename=ps

	legend1x = [1213,1214.5]
	legend1y = [620,620]
	legend2y = [570,570]

        plot, lambda, flux, psym=7, yrange=[-50, 650], xrange=[1212,1220], ystyle=1, color=0, position=[0.14, 0.14, 0.9, 0.9], ytitle= textoidl('f_{\lambda} (Arbitrary Units)'),xstyle=4, xthick=6, ythick=6, thick=4, charsize = 1.5, charthick=5, /noerase

        oplot, x, bestogauss, color=5, linestyle=2, thick=4
        oplot, legend1x, legend1y, color=5, linestyle=2, thick=4
        xyouts, 1212.45, 590, 'Original Gaussian', color=5, charsize=1.5, charthick=3, /data

        oplot, x, bestgauss, color=6, linestyle=0, thick=4
        oplot, legend1x, legend2y, color=6, linestyle=0, thick=4
        xyouts, 1212.25, 540, 'Truncated Smoothed', color=6, /data, charsize=1.5, charthick=3

	xyouts, 1213.15, 515, 'Gaussian', color=6, /data, charsize=1.5, charthick=3
        ;xyouts, 1217.5, 00, textoidl('\lambda_{truncation} =') + strcompress(string(besttrunc), /remove_all), color=0, /data
        ;xyouts, 1217.5, 850, textoidl('\sigma =') + string(bestsig), color=0, /data
        ;xyouts, 1217.5, 800, textoidl('A =') + string(bestamp), color=0, /data 
        ;xyouts, 1217.5, 750, textoidl('C =') + string(bestback), color=0, /data
        sigv = (bestsig/bestlam)*3e+05
	xyouts, 1216.4, 575, textoidl('DEIMOS FWHM = ') + strcompress(string(number_formatter(bestconvol*0.33*2.35, decimal=2)), /remove_all) + textoidl(' \AA'), color=0, charsize=1.4, charthick=4, /data
        xyouts, 1216.4, 540, textoidl('\lambda_{c} = ') + strcompress(string(bestlam), /remove_all) + textoidl(' \AA'), color=0, charsize=1.4, charthick=4, /data
        xyouts, 1216.4, 475, textoidl('\sigma_{v} = ') + strcompress(string(number_formatter(sigv, decimal=2)), /remove_all) + ' (km/s)', color=0, charsize=1.4, charthick=4, /data

        axis, xaxis=0, xstyle=1, xrange=[1212, 1220], xthick=6, charsize=1.5, charthick=5, /save, color=0, xtitle=textoidl('\lambda_{rest} (\AA)')
        axis, xaxis=1, xstyle=1,xrange=[-864,1111], xthick=6, charsize=1.5, charthick=5, /save, color=0, xtitle=textoidl('velocity (km/s)')	;-1367, 1604
        device, /close_file	
        set_plot, entry_device

end
