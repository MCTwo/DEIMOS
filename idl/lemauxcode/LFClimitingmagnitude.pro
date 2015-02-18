;
;
;
;   Name
;	LFClimitingmagnitude.pro
;
;   Pupose
;	To calculate a (N) sigma limiting magnitude for SC1604 objects assuming a circular aperture. This
;	assumes (as is typical with limiting magnitudes) that each pixel in the circular aperture has
;	a value of (N) sigma as calculated in a statistically signifcant region near or on top of the 
;	object.
;
;   Inputs
;  	
;	rsigma = The rms value calculated in the r' image in a statistically significant region near or on
;	top of the object.
;	isigma = i' counterpart to rsigma
;	zsigma = z' counterpart to rsigma
;	r = the effective radius of the galaxy (or an estimate of it) that you are calculating the limiting
;	magnitude for
;	N = the sigma value desired for the calculation (typically 3 sigma limiting magnitudes are quoted)
;	pointing = LFC/Hale telescope pointing either SC1 or SC2 (these are the only two LFC pointings for 
;	SC1604)
;
;   Outputs
;	rmagold, imagold, zmagold = the magnitudes in the old system using the original RG zero points, used to 
;	calculate SDSS calirated magnitudes
;	rsdss, isdss, zsdss = magnitudes calibrated to the SDSS system, these are the actual (N) sigma limiting 
;	magnitudes
	
pro LFClimitingmagnitude, rsigma, isigma, zsigma, r, N, pointing

	if n_elements(pointing) then pointing = pointing[0] else begin
		pointing = 'SC1'
		print, 'No pointing specified, assuming the object is located in the first pointing!'
	endelse
	
	if pointing ne 'SC1' and pointing ne 'SC2' then begin 
		message, 'These transforms are only valid for SC1604 LFC pointings SC1 and SC2!'
	endif

	if n_elements(N) then N=N[0] else begin
	N = 3
	print, 'No threshold set for limiting magnitude, assuming you want 3 sigma'
	endelse

	if n_elements(r) then r = r[0] else begin
	r = 0.5
	print, 'No galaxy size specified, assuming you want r_eff = 0.5"'
	endelse

	A = !PI*(r/0.18)^2	;area subtended in pixels (for LFC 0.18"/pix scale) of a circle of radius r

	rmagold = 26.83014 - 2.5*alog10(A*N*rsigma)
	imagold = 26.52275 - 2.5*alog10(A*N*isigma)
	zmagold = 25.46807 - 2.5*alog10(A*N*zsigma)	

	if pointing eq 'SC1' then begin
	
		rsdss = 1.010709*rmagold + 0.00014*(rmagold - imagold) -0.44597 
		isdss = 0.986003*imagold +0.00139*(rmagold - imagold) +0.07397
		zsdss = 0.982203*zmagold + 0.000004*(imagold-zmagold)+0.005811
	endif

	if pointing eq 'SC2' then begin 
	
		rsdss = 1.0138846*rmagold -0.00114*(rmagold - imagold) -0.43894
		isdss = 0.98598*imagold +0.010645*(rmagold - imagold) -0.040641
		zsdss = 0.989737*zmagold -0.00373*(imagold-zmagold)-0.44112
	endif

	print, '      r_old          i_old          z_old          r_sdss         i_sdss         z_sdss'
	print, rmagold, '  ', imagold, '  ', zmagold, '  ', rsdss, '  ', isdss, '  ', zsdss
	
end	
