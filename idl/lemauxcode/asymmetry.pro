;for calculating the asymmetry of a Lya emission line (or effective Lya line) in the manner of Dawson et al. 2004

pro loadcolors, bottom=bottom, names=names

        if (n_elements(bottom) eq 0) then bottom=0

        ;Load Graphics Colors
        red = [ 0, 255, 0, 255, 0, 255, 0, 255, $
                0, 255, 255, 112, 219, 127, 0, 255]
        grn = [ 0, 0, 255, 255, 255, 0, 0, 255, $
                0, 187, 127, 219, 112, 127, 163, 171]
        blu = [ 0, 255, 255, 0, 0, 0, 255, 255, $
                115, 0, 127, 147, 219, 127, 255, 127]
        tvlct, red, grn, blu, bottom

        ; Set color names

        names = ['Black', 'Magenta', 'Cyan', 'Yellow', $
                 'Green', 'Red', 'Blue', 'White', $
                 'Navy', 'Gold', 'Pink', 'Aquamarine', $
                 'Orchid', 'Gray', 'Sky', 'Beige' ]
end

pro asymmetry, cat, outcat, width

	readcol, cat, format='A,A,D,D,A', mask, slit, z, q, file, skipline=1
	
	openw, lun, outcat, width=400, /get_lun

	for i=0, n_elements(file)-1 do begin
	spec = file[i]
        spec = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/'+spec
	
	corrspec = fill_gap(spec,/tweak,/telluric,/silent,header=header)
	nonreflambda = corrspec.lambda
	wave = corrspec.lambda
        airtovac, wave
        corrspec.lambda = wave
	
	initlam = 1215.7*(z[i]+1)			;can always fake this, the -2 is from airtovac
	maxrange = where(corrspec.lambda ge initlam-2 and corrspec.lambda le initlam+2)			; 4 Angstrom window centered on the initial guess to search for the max valus. Made this really small so that it didn't pick up noise spikes
	range = where(corrspec.lambda ge initlam-width/2 and corrspec.lambda le initlam+width/2)	;set search range for flux analysis, for red and blue 10% flux values
	rangelam = corrspec.lambda[range]
	rangespec = corrspec.spec[range]
	;splot, corrspec.lambda[range], corrspec.spec[range]
	maxrangespec = corrspec.spec[maxrange]
	maxspec = where(rangespec eq max(maxrangespec))	;array value of max flux value in the bandpass
	lamc = rangelam[maxspec]
	maxflux = rangespec[maxspec]
	newz = lamc/1215.7-1
	print, lamc, initlam
	print, file[i]

	;OLD	
	if abs(rangelam[maxspec]-initlam) gt 2 then begin
	loadcolors
	print, abs(lamc-initlam)
	splot, rangelam, rangespec
	oplot, [lamc[0],lamc[0]], [-500, 500], color=4, thick=2
        oplot, [initlam[0], initlam[0]], [-500, 500], color=2, thick=2	
	message, "Initial wavelength guess is off by change redshift input!"
	endif
	;#################	
	bluerange = where(rangelam lt lamc[0])	;blue range between inital central wavlength minus half the width and the new calculated central wavelength...gonna be fucked up if intial guess is way off
	redrange = where(rangelam gt lamc[0])
	
	bluespec = rangespec[bluerange]		;flux array blueward of the line
	redspec = rangespec[redrange]
	bluelam = rangelam[bluerange]
	redlam = rangelam[redrange]	

	blue10 = where(bluespec le 0.1*maxflux[0])	;array values which have flux less than 10% of the central flux on the blue side of the line
	red10 = where(redspec le 0.1*maxflux[0])	;ibid, but red side

	if blue10[0] lt 0 then begin
	print, maxflux
	print, min(bluespec)
	splot, rangelam, rangespec	
	message, "Increase the width, no flux values found with 10% of the central flux on the blue side"
	endif
	
	if red10[0] lt 0 then begin 
	print, maxflux
	print, min(redspec)
	splot, rangelam, rangespec
	message, "Increase the width, no flux values found with 10% of the central flux on the red side"	
	endif
	
	if blue10[0] ge 0 then lambda10blue = bluelam[blue10[n_elements(blue10)-1]] else lambda10blue = 0		;get the last element of the array for which the flux is 10% or less than the central flux blueward of the line. This is the reddest value blueward of the line with flux <= 10%
	if red10[0] ge 0 then lambda10red = redlam[red10[0]] else lambda10red = 0		;get the first element of the array for which the flux is 10% or less than the central flux redward of the line. This is the bluest value redward of the line with flux <= 10%.
	
	asymmetry = (lambda10red - lamc)/(lamc - lambda10blue)
	printf, lun, mask[i], '   ', slit[i], '    ', z[i], '   ', q[i], '   ', asymmetry[0], '   ', file[i]

	loadcolors
	splot, rangelam, rangespec, psym=10
	oplot, [bluelam[blue10[n_elements(blue10)-1]],bluelam[blue10[n_elements(blue10)-1]]], [-5000,5000], color=6, thick=2
	oplot, [redlam[red10[0]],redlam[red10[0]]], [-5000,5000], color=5, thick=2
	oplot, [lamc[0],lamc[0]], [-5000, 5000], color=4, thick=2
	oplot, [initlam[0], initlam[0]], [-5000, 5000], color=2, thick=2
	endfor

	free_lun, lun

end
	
