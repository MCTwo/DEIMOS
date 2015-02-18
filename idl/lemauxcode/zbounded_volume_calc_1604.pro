function Ez, z

   distance = 1/sqrt(0.3*(1+z)^3+0.7)
   return, distance

end

pro zbounded_volume_calc_1604, skytol, zmin, zmax, cplot

	if n_elements(zmin) then zmin = zmin[0] else zmin = 0
	if n_elements(zmax) then zmax = zmax[0] else zmax = 1000

	;note skytol is normally set to 60
	readcol, '1604lambdahist.dat', format='D,D', pixlam, pixnum
	
	skytemp = mrdfits('uves_sky.fits',1)
	
	skyarray = where(skytemp.spec gt skytol)
	skylam = skytemp.lambda[skyarray]
	skylammin = skylam - 1.7
	skylammax = skylam + 1.7
	
	lambda = 6100 + findgen(6600)*0.5
	weightarray = fltarr(6600)+1 	;start off with all wavelengths being completely exposed, no sky interference 
	zweightarray = fltarr(6600)+1	;start off with all redshifts being valud
	
	b = where(pixlam lt (zmin+1)*1215.7 or pixlam gt (zmax+1)*1215.7)
	if b[0] ge 0 then zweightarray[b] = 0 	;prevent crashing if there are no redshift limits
	
	for i=0, n_elements(skylam)-1 do begin

	a = where(lambda gt round(skylammin[i]*2)/2 and lambda lt round(skylammax[i]*2)/2) ;rounding so can get it to the nearest 0.5 A bin
	
	if a[0] ge 0 then weightarray[a] = 0.	;prevent crashing when skylines are out of the spectral range
	
	
	endfor
	
	newpixnum = pixnum*weightarray*zweightarray
	loadcolors
	
	Dh = 3000./0.7		;in Mpc can change cosmology here and in Ez, this is h=0.7
	newangle = newpixnum/max(pixnum)*3.08/3600*(2*!PI/360)^2	;in square arcminutes then converted to radians, number exposed for each wavelength
	angle = pixnum*zweightarray/max(pixnum)*3.08/3600*(2*!PI/360)^2
	zpix = (pixlam/1215.7)-1
	
	avglam = fltarr(n_elements(pixlam)-1)	
	avgz = fltarr(n_elements(pixlam)-1)
	Dakinda = fltarr(n_elements(pixlam)-1)
	DLkinda = fltarr(n_elements(pixlam)-1)
	volume = fltarr(n_elements(pixlam)-1)
	newvolume = fltarr(n_elements(pixlam)-1)
	newunitvolume = fltarr(n_elements(pixlam)-1)
	unitvolume = fltarr(n_elements(pixlam)-1)
	DLkindaunitz = fltarr(n_elements(pixlam)-1)
		
	for i=0, n_elements(newpixnum)-2 do begin	;-2 is OK cuz there are a lot of zeros at the end, not truncating any volume
	
		avglam[i] = (pixlam[i] + pixlam[i+1])/2
		avgz[i] = ((zpix[i]+zpix[i+1])/2)
		Dakinda[i] = (qromb('Ez', 0., avgz[i]))^2
		DLkinda[i] = qromb('Ez',zpix[i], zpix[i+1])	
		volume[i] = angle[i]*Dh^3*Dakinda[i]*Dlkinda[i]	;volume for each wavelength
		newvolume[i] = newangle[i]*Dh^3*Dakinda[i]*Dlkinda[i]	;volume for each wavelength incorperating the obstruction from sky lines
		DLkindaunitz[i] = qromb('Ez',zpix[i]-0.5, zpix[i+1]+0.5)
		newunitvolume[i] = newangle[i]*Dh^3*Dakinda[i]*Dlkindaunitz[i]	
		unitvolume[i] = angle[i]*Dh^3*Dakinda[i]*Dlkindaunitz[i]
	endfor	

	cumvol = fltarr(n_elements(volume))
	newcumvol = fltarr(n_elements(volume))
	for i=0, n_elements(volume)-1 do begin
		cumvol[i] = total(volume[0:i])
		newcumvol[i] = total(newvolume[0:i])
	endfor

	lamrange = where(avglam lt 9250)		
	if n_elements(cplot) then cplot = cplot[0] else cplot = 0

	if cplot[0] ne 0 then begin
	entry_device = !d.name
	set_plot, 'PS'	
	device, /color, xoffset=0.2, yoffset=0.2, xsize=8, ysize=10, /inches, filename='higherzLAEcluster.z480to486.ps'
	
	plot, pixlam, unitvolume, /nodata, position=[0.17, 0.55, 0.9, 0.9], color=0, xrange = [7025,7150], $
        ytitle= textoidl('volume / \Delta z (Mpc^{3})'), xstyle=4, ythick=6, charsize=1.5, charthick=6, /noerase
	axis, xaxis=0, xrange=[7025,7150], xstyle=1, /save, xtickname=replicate(' ', 8), color=0, xthick=6;,xrange=[pixlam[0],pixlam[n_elements(pixlam)-1]]
        oplot, pixlam, newunitvolume, color=5, thick=4
        oplot, pixlam, unitvolume, color=6, thick=4
	axis, xaxis=1, xstyle=1, xrange=[4.80,4.86], xthick=6, color=0, charsize=1.5, charthick=6 
	xyouts, 0.51, 0.94, textoidl('z_{Ly\alpha}'), /normal, color=0, charsize=1.5, charthick=6

	plot, avglam[lamrange], cumvol[lamrange], /nodata, position=[0.17, 0.16, 0.9, 0.55], color=0, yrange=[0,6.2e+02], ystyle=1,$
	xrange = [7025,7150], xstyle=1, xtitle=textoidl('\lambda_{Obs} (\AA)'), ytitle= textoidl('Cumulative volume (Mpc^{3})'), $
	ythick=6, xthick=6, charsize=1.5, charthick=5, /noerase
	oplot, avglam[lamrange], cumvol[lamrange], color=6, linestyle=2, thick=4
	oplot, avglam[lamrange], newcumvol[lamrange], color=5, thick=4
	device, /close_file
	set_plot, entry_device
	endif
	
	;splot, pixlam, unitvolume, color=1
	splot, pixlam, newunitvolume, color=5
	print, total(volume), total(newvolume)
end
