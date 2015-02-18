;This is a program to generate a series of mock observations of the DEEP2 data observed by me 
;in order to determine the distribution of how many LAE are observed when the volume is 10x smaller (or 83/12X).
;Assuming that each mask has 135. slits (since this is a ratio it's probably ok) and drawing 12*135 slits out of
;the possible distribution of 83*135 slits. Also assuming each slit has the same length, this is complimentary 
;to the other assumption since each mask samples the same area of the sky (not including loss from targets which
;grows linearly with the number of slits)

pro sample_variance_sim_forLAEclusters, niter

	;List of all DEEP2 masks zspeced

mask = 	[1100., 1101., 1103., 1112., 1140, 1143, 1146, 1150, 2103, 2105, 2106, 3307, 3346, 4101, 4103, 4104, 4106, 4107, 4109, 4110, 4112, 4115, 4118, 4143, 4144, 4146, 4147, 4149, 4150, 4152, 4153, 4155, 4156, 4158, 4180, 4200, findgen(22)+4201, findgen(20)+4240, findgen(3)+4260, 4280, 4281]

mask = mask*1000			;so that I can seperate out mask/slit names and so every name is unique
slit = (findgen(135)+1)

LAEmasklist = [1112., 1143, 1150, 2103, 3307, 4107, 4110, 4112, 4147, 4180, 4218, 4240, 4243, 4243, 4256, 4257, 4260, 4280, 4280]

LAEmasklist = LAEmasklist*1000

LAEslitlist = [73., 85, 116, 71, 61, 130, 123, 49, 18, 132, 94, 64, 121, 126, 30, 9, 25, 73, 76]

objarray = fltarr(83*135)			;create a array with values XXXXYYY.00 mask=X slit=Y

	for i=0, 83-1 do begin			
		for j=0, 135-1 do begin	
		objarray[j+135*i] = mask[i] + slit[j]
		endfor
	endfor

LAEarray = fltarr(19)				;do the same for the slits that contain LAEs 

LAEarray = LAEmasklist+LAEslitlist

countarr = intarr(niter)
	for i=0, 19-1 do begin
	a = where(objarray eq LAEarray[i])
	endfor	
	
	for i=0, niter-1 do begin		;do niter number of mock observations
	
		dummyarray = objarray		;set up a dummy array so that objarray can be modified
		observedarray = fltarr(round(0.58692*135))
		for j=0, round(0.58692*135)-1 do begin
			randomobs = round(randomu(seed, 1, /uniform)*(83*135-[j+1])) ;generate a random number between 0 and the size of dummy array
			observedarray[j] = objarray[randomobs]	; select an observed mask/slit 
			dummyarray[randomobs] = 0. 		; set the observed array element (mask/slit) to zero
			dummyarray = dummyarray(where(dummyarray ne 0.))	; re-create the list of potential observed mask/slits with the previous observation removed so there are no duplicate observations
		endfor		

		;countarr = intarr(LAEslitlist)
		;for j=0, n_elements(LAEslitlist)-1 do begin
		count = fltarr(19)
		realcount = fltarr(19)
		obsLAEarray = fltarr(19)
		for j=0, 19-1 do begin
				count[j] = where(observedarray eq LAEarray[j]) 	; count how many times a LAE was observed
				if count[j] ne -1 then obsLAEarray[j] = observedarray[count[j]] else obsLAEarray[j] = 0.
				if count[j] ne -1 then realcount[j] = 1 else realcount[j] = 0
		endfor		
		
		countarr[i] = total(realcount)		
	
		if countarr[i] ne 0 then obsLAEarray = obsLAEarray(where(obsLAEarray ne 0.)) else obsLAEarray = obsLAEarray[0]		;if there are some LAEs detectected save these to be printed out, if not then save 0 to be printed out
		;print, obsLAEarray, format='(f11.2)'
		countarr[i] = total(realcount)				
		if i eq 1000 then print, 'Iteration 1000'
		if i eq 2000 then print, 'Iteration 2000'
		if i eq 3000 then print, 'Iteration 3000'
		if i eq 4000 then print, 'Iteration 4000'
		if i eq 5000 then print, 'Iteration 5000'
		if i eq 6000 then print, 'Iteration 6000'
		if i eq 7000 then print, 'Iteration 7000'
		if i eq 8000 then print, 'Iteration 8000'
		if i eq 9000 then print, 'Iteration 9000'
	endfor

openw, lun, 'observedLAEcounts.4.80z4.86.dat', /get_lun
	printf, lun, countarr, format='(f10.3)'
free_lun, lun

loc = fltarr(20)
obsLAEnumberhist = histogram(countarr, binsize=1, location=loc, max=19, min=0)
splot, loc, obsLAEnumberhist


end


