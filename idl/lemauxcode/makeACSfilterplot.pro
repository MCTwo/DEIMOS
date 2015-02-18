pro makeACSfilterplot
	
	loadcolors
	readcol, 'wfc_f606w_t77.dat', lam606, through606
	readcol, 'wfc_f814w_t77.dat', lam814, through814
	set_plot, 'PS'
	device, /color, xoffset=0.2, yoffset=0.2, xsize=8, ysize=10, /inches, filename='ACSfiltercurves.ps'
	plot, lam606, through606, /nodata, position=[0.17, 0.2, 0.9, 0.85], color=0, xrange = [4000,10000],$
	xtitle=textoidl('\lambda (\AA)'), ytitle= textoidl('Total throughtput (HST+ACS)'), charsize=1.6
	oplot, lam606, through606, color=6
	oplot, lam814, through814, color=5
	xyouts, 8300, 0.43, 'F814W', charsize=1.5, color=0, /data
	xyouts, 4400, 0.43, 'F606W', charsize=1.5, color=0, /data
	device, /close_file
end	
