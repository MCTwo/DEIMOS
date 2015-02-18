pro makeDEEP2throughputps
	
	ps = 'DEEP2throughput.ps'
	;readcol, 'thr_DEIMOS_go1200_8000Ang_og550.dat', format = 'D,D', wredest, thrredest
	;readcol, 'thr_DEIMOS_go1200_7000Ang_og550.dat', format = 'D,D', wred, thrred
	wnick = findgen(10192)*0.33 + 6000
	thrnick = 1/deimos_correction(wnick)
	loadcolors 
	set_plot, 'PS'
	device, /color, bits_per_pixel=24, xsize=8, ysize=6, /inches, filename=ps
	plot, wnick, thrnick, /nodata, xrange=[5700,9500], xstyle=1, yrange=[0,0.4], ystyle=1, xthick=6, xtitle = textoidl('\lambda_{Obs}',font=0), font=0, ytitle= 'DEIMOS Throughput', color=0, thick=4, ythick=6, charthick=6, charsize=1.5
	oplot, wnick, thrnick, color=6, thick=5
		
	;legend
	xleg = findgen(1200)*0.33 + 7500.
	ymasblue = fltarr(1200)+0.08
	ybluest = fltarr(1200)+0.07
	yblue = fltarr(1200) + 0.06
	yred = fltarr(1200) + 0.05
	yredest = fltarr(1200) + 0.04
	ynick = fltarr(1200) + 0.03
	ybrian = fltarr(1200) + 0.02
	;oplot, xleg, ynick, color=1
	;xyouts, 8000., 0.027, 'Nick DEIMOS Response', color=1, charsize=1, /data
	
	; feature ranges for various structures

	x1604 = 7332
	y1604 = 0.18
	x1604err = 776
	yerr = 0
	oploterror, x1604, y1604, x1604err, yerr,  errcolor=5, psym=3, errthick=6
	xyouts, x1604[0]-x1604err[0]+350, y1604[0]-0.024, textoidl('CL1604 Ly\alpha'), color=5, charsize=1.5, charthick=5

        device, /close_file
	

end
