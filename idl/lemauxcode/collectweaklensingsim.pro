pro collectweaklensingsim, nint
	
	dir = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/1604_LYA/weaklensing/realizations'
	
	magfluxcutarray = fltarr(nint,24)
	fluxcutarray = fltarr(nint,24)

	magarray = fltarr(nint,24)
	array = fltarr(nint,24)

	number = sindgen(nint)
	number[0:9] = '000' + number[0:9]
	number[10:99] = '00' + number[10:99]
	number[100:nint-1] = '0' + number[100:nint-1] 	
	number = strcompress(number, /remove_all)

	for i=0, nint-1 do begin

	readcol, dir + '/ndgrid_run_' + number[i] + '.txt', format = 'D,D,D', cutnums, bleep, nums
	readcol, dir + '/magndgrid_run_' + number[i] + '.txt', format = 'D,D,D', cutmagnums, blorp, magnums
	fluxcutarray[i,0:23] = cutnums[0:23]
	magfluxcutarray[i,0:23] = cutmagnums[0:23]
	array[i,0:23] = nums[0:23]
	magarray[i,0:23] = magnums[0:23] 

	endfor

	meanmag = fltarr(24)
	stdmag = fltarr(24) 
	meanreg = fltarr(24)
	stdreg = fltarr(24)
	meancutmag = fltarr(24)
        stdcutmag = fltarr(24)
        meancutreg = fltarr(24)
        stdcutreg = fltarr(24)

	lum = 10^(0.175+findgen(25)*0.075)
	lum2 = 10^(0.250 + findgen(25)*0.075)
	lumavg = (lum+lum2)/2
	openw, lun, 'weaklensing_numbers.cat', width = 200, /get_lun	
	for i=0, 24-1 do begin
	
	meanmag[i] = mean(magarray[0:nint-1,i])
	stdmag[i] = stddev(magarray[0:nint-1,i])
	stdreg[i] = stddev(array[0:nint-1,i])
	meanreg[i] = mean(array[0:nint-1,i])

	meancutmag[i] = mean(magfluxcutarray[0:nint-1,i])
        stdcutmag[i] = stddev(magfluxcutarray[0:nint-1,i])
        stdcutreg[i] = stddev(fluxcutarray[0:nint-1,i])
        meancutreg[i] = mean(fluxcutarray[0:nint-1,i])


	printf, lun, lumavg[i], '   ', meanreg[i], '   ', stdreg[i], '    ', meanmag[i], '   ', stdmag[i], '   ', meancutreg[i], '   ', stdcutreg[i], '    ', meancutmag[i], '   ', stdcutmag[i]
	endfor
	
	free_lun, lun
	
end
