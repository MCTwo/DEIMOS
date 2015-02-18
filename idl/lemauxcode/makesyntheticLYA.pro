function makesyntheticLYA, z

	result = -1 
	anarray = findgen(19697)*0.33 + 4500	;lambda array from 4500 to 11000
	fluxarray = fltarr(n_elements(anarray))	;associated flux array
	
	;finding various limits
	obslyalam = 1215.7*(1+z)
	bluelam = (1215.7-1.3)*(1+z)
	obslylimitlam = (912+3.5)*(1+z)
	bluerlam = (912-3.5)*(1+z)
	redlam = (1215.7+5.5)*(1+z)
	
	magbreak = 3.8 + 20.3*alog10((1+z)/7)	;magnitude break over the Lya line as defined in Hu et al 2004
	fluxbreak = 10^(-magbreak/2.5)	;corresponding flux break
	lyaheight = 1000	;height of the Lya line above the red continuum
	redcontheight = 1200 	;continuum level redward of the Lya line, this is the basis for all other flux numbers, everything is relative to this
	;populating flux array

	redcontarr = where(anarray ge redlam)
	redlyarr = where(anarray ge obslyalam and anarray lt redlam)	
	bluelyarr = where(anarray ge bluelam and anarray lt obslyalam)
	bluecontarr = where(anarray ge obslylimitlam and anarray lt bluelam)
	lylimitarr = where(anarray ge bluerlam and anarray lt obslylimitlam)
	bluercontarr = where(anarray lt bluerlam)	
	

	if redcontarr[0] ne -1 then fluxarray[redcontarr] = redcontheight 
	if redlyarr[0] ne -1 then fluxarray[redlyarr] = (-lyaheight)/(redlam - obslyalam) * anarray[where(anarray ge obslyalam and anarray lt redlam)] + (lyaheight+redcontheight+redcontheight)/2 + (lyaheight)/(redlam - obslyalam)*(redlam+obslyalam)/2	; transition from Lya line peak to redside continuum; mx+b, m=delta_y/delta_x, b = <v> - m*<x>
	if bluelyarr[0] ne -1 then fluxarray[bluelyarr] = (lyaheight+redcontheight - redcontheight*fluxbreak)/(obslyalam-bluelam) * anarray[where(anarray ge bluelam and anarray lt obslyalam)] + (lyaheight+redcontheight + redcontheight*fluxbreak)/2 - (lyaheight+redcontheight - redcontheight*fluxbreak)/(obslyalam-bluelam)*(obslyalam+bluelam)/2 ;transition from Lya peak to blueside continuum
	if bluecontarr[0] ne -1 then fluxarray[bluecontarr] = 1200*fluxbreak ;blueside continuum
	if lylimitarr[0] ne -1 then fluxarray[lylimitarr] = ((1200*fluxbreak) - 30)/(obslylimitlam - bluerlam) * anarray[where(anarray ge bluerlam and anarray lt obslylimitlam)] + ((1200*fluxbreak) + 30)/2 - ((1200*fluxbreak) - 30)/(obslylimitlam - bluerlam)*(obslylimitlam + bluerlam)/2	+ 2000;transition over the Lyman limit, adding a bunch of noise
	if bluercontarr[0] ne -1 then fluxarray[bluercontarr] = 2000	;continuum blueward of the Lyman limit, only keeping this non-zero to add in noise

	noisyfluxarray = dblarr(n_elements(fluxarray))
	for i=0, n_elements(fluxarray)-1 do begin
		noisyfluxarray[i] = randomu(seed, poisson=fluxarray[i], /double)
	endfor
	
	if lylimitarr[0] ne -1 then noisyfluxarray[lylimitarr] = noisyfluxarray[lylimitarr] - 2000
	if bluercontarr[0] ne -1 then noisyfluxarray[bluercontarr] = noisyfluxarray[bluercontarr] - 2000 ;subtracting off most of the background blueward of the Lyman limit
	
	sigma = 1.4
        spec = {spec:noisyfluxarray, lambda:anarray, ivar:1/abs(noisyfluxarray)}
	smoothnoisyspec = smooth_spec(spec, sigma)
	smoothnoisyspec.spec[n_elements(smoothnoisyspec.spec)-5:n_elements(smoothnoisyspec.spec)-1] = mean(smoothnoisyspec.spec[redcontarr]) 	;cleaning up last few elements
	;splot, smoothnoisyspec.lambda, smoothnoisyspec.spec
	
	result = smoothnoisyspec
	return, result
	
end
