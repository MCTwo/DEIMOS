pro vegafluxcheck, filter
	
readcol, 'vega.full.dat', format='D, D', lambda, flambda 

readcol, filter + '.dat', format = 'D,D,D', lam, pointthrough, extendedthrough, /SILENT
minlam = min(lam)       ;min wavelength of the filter transmission curve
maxlam = max(lam)	

spec_in = flambda
lambda_in = lambda	


dummy = where(lambda_in ge minlam and lambda_in le maxlam)      ;check where the actual spectrum overlaps the filter transmission curve
arraylength = n_elements(dummy)
maxlambda = lambda_in[dummy[arraylength-1]]
minlambda = lambda_in[dummy[0]]
througharray = generate_interpol_filtercurve(filter, arraylength, minlambda, maxlambda)	;interpolate values so that the interpolated grid matches up with the spectrum grid, both in intial wavelength values and spacing

deltalamband = fltarr(n_elements(dummy))
meanlam = fltarr(n_elements(dummy))
           	
		for j=0, n_elements(dummy)-2 do begin
                deltalamband[j] = abs(lambda_in[dummy[j]]-lambda_in[dummy[j+1]])   ;the frequency interval that the filter throughput covers, needed to get a flux density for converting to AB magnitudes 
                meanlam[j] = (lambda_in[dummy[j]] + lambda_in[dummy[j+1]])/2	;get the effective wavelength in each bin
		endfor


dummy2 = where(deltalamband eq 0)		;clean up last element
deltalamband[dummy2] = mean(deltalamband)
dummy3 = where(meanlam eq 0)
meanlam[dummy3] = max(meanlam)
bandfluxnum = fltarr(n_elements(dummy))
bandfluxdenom = fltarr(n_elements(dummy))

	        for j=0, n_elements(dummy)-1 do begin
                bandfluxnum[j] = througharray[j]*spec_in[dummy[j]]*deltalamband[j]*meanlam[j]/(3.0D+18)	;this is pedantic see Fukugita 1995 eqn 9, F_nu is corr_spec_in/deltafband
		bandfluxdenom[j] = througharray[j]*deltalamband[j]/meanlam[j]
                endfor

totbandfluxdens = total(bandfluxnum)/total(bandfluxdenom)


ZP = 48.60	;photometric zero point for AB system
extinction = -0.1	;approximate extinction in magnitudes/airmass at Mauna Kea, see http://www.cfht.hawaii.edu/Instruments/ObservatoryManual/CFHT_ObservatoryManual_(Sec_2).html	
filterABmag = -ZP - 2.5*alog10(totbandfluxdens) ;+ airmass*extinction
print, 'The AB magnitude of Vega is 0.4 in the', filter, 'band, you have calculated', filterABmag
	
end
		
