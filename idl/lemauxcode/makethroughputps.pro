pro makethroughputps
	
	ps = 'throughput.ps'
	readcol, 'thr_DEIMOS_go1200_8000Ang_og550.dat', format = 'D,D', wredest, thrredest
	readcol, 'thr_DEIMOS_go1200_7000Ang_og550.dat', format = 'D,D', wred, thrred
	readcol, 'thr_DEIMOS_go1200_7000Ang_gg455.dat', format = 'D,D', wblue, thrblue
	readcol, 'thr_DEIMOS_go1200_6000Ang_gg455.dat', format = 'D,D', wbluest, thrbluest
	readcol, 'DEIMOS_throughput_6000+7000_gg455.dat', format = 'D,D', wmasterblue, thrmasterblue
	wnick = findgen(10192)*0.33 + 6000
	thrnick = 1/deimos_correction(wnick)
	wbrian = findgen(11192)*0.33 + 5200
	thrbrian = 1/deimos_correction_masterblue(wbrian)	
	set_plot, 'PS'
	device, /color, bits_per_pixel=24, filename=ps
	plot, wmasterblue, thrmasterblue, xrange=[4500,9500], xstyle=1, yrange=[0,0.3], ystyle=1, xtitle = textoidl('\lambda',font=0), font=0, ytitle= 'DEIMOS Throughput', color=0, thick=4
	oplot, wbluest, thrbluest, thick=1, color=6
	oplot, wredest, thrredest, color=5
	oplot, wnick, thrnick, color=1
	oplot, wblue, thrblue, color=4, thick=1
	oplot, wred, thrred, color=2
	oplot, wbrian, thrbrian, color=13	
	;legend
	xleg = findgen(1200)*0.33 + 7500.
	ymasblue = fltarr(1200)+0.08
	ybluest = fltarr(1200)+0.07
	yblue = fltarr(1200) + 0.06
	yred = fltarr(1200) + 0.05
	yredest = fltarr(1200) + 0.04
	ynick = fltarr(1200) + 0.03
	ybrian = fltarr(1200) + 0.02
	oplot, xleg, ymasblue, color=0
	xyouts, 8000., 0.077, 'GG455+6000A+7000A', color=0, charsize=1, /data
	oplot, xleg, ybluest, color=6
	xyouts, 8000., 0.067, 'GG455+6000A', color=6, charsize=1, /data	
	oplot, xleg, yblue, color=4
	xyouts, 8000., 0.057, 'GG455+7000A', color=4, charsize=1, /data
	oplot, xleg, yred, color=2
	xyouts, 8000., 0.047, 'OG550+7000A', color=2, charsize=1, /data
	oplot, xleg, yredest, color=5
	xyouts, 8000., 0.037, 'OG550+8000A', color=5, charsize=1, /data
	oplot, xleg, ynick, color=1
	xyouts, 8000., 0.027, 'Nick DEIMOS Response', color=1, charsize=1, /data
	oplot, xleg, ybrian, color=13
	xyouts, 8000., 0.017, 'Brian DEIMOS Response', color=13, charsize=1, /data
	
	; feature ranges for various structures

	y1604 = [0.19,0.19]
	x1604 = [7080., 7792]
	x1604err = [224., 246]
	yerr = [0.,0.]
	oploterror, x1604, y1604, x1604err, yerr,  errcolor=4, psym=3, color=0
	xyouts, x1604[0]-x1604err[0]-50, y1604[0]-0.015, '1604 OII', color=4, /data
	xyouts, x1604[1]-x1604err[1]-50, y1604[1]-0.015, textoidl('1604 H\delta',font=0), color=4, font=0, /data	

	y1322 = [0.15, 0.15]
	x1322 = [6460., 7115.]
	x1322err = [160., 185.]
	oploterror, x1322, y1322, x1322err, yerr,  errcolor=5, psym=3, color=0
	xyouts, x1322[0]-x1322err[0]-50, y1322[0]-0.015, '1324 OII', color=5, /data
	xyouts, x1322[1]-x1322err[1]-50, y1322[1]-0.015, textoidl('1324 H\delta',font=0), color=5, font=0, /data

	y5281 = [0.17, 0.17]
        x5281 = [6760., 7435.]
        x5281err = [50., 55.]
        oploterror, x5281, y5281, x5281err, yerr,  errcolor=0, psym=3, color=0
        xyouts, x5281[0]-x5281err[0]-200, y5281[0]-0.015, '5281 OII', color=0, /data
        xyouts, x5281[1]-x5281err[1]-200, y5281[1]-0.015, textoidl('5281 H\delta',font=0), color=0, font=0, /data

	y2560 = [0.13, 0.13]
        x2560 = [6020., 6620.]
        x2560err = [50., 60.]
        oploterror, x2560, y2560, x2560err, yerr,  errcolor=6, psym=3, color=0
        xyouts, x2560[0]-x2560err[0]-200, y2560[0]-0.015, '2560 OII', color=6, /data
        xyouts, x2560[1]-x2560err[1]-200, y2560[1]-0.015, textoidl('2560 H\delta',font=0), color=6, font=0, /data
        device, /close_file
	

end
