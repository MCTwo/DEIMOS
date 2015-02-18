function Ez, z

   distance = 1/sqrt(0.3*(1+z)^3+0.7)
   return, distance

end

function Tz, z
	
   time = 1/((1+z)*sqrt(0.3*(1+z)^3+0.7))
   return, time

end
		
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

pro comoving_volume, h, z1, z2, angle, cplot=cplot

	if n_elements(h) then h = h[0] else begin
		message, 'You must specify a cosmology!'
	endelse
	if n_elements(z1) then z1 = z1[0] else z1 = 0	
	if n_elements(z2) then z2 = z2[0] else z2 = 0
	if n_elements(angle) then angle = angle[0] else angle = 0
	if n_elements(cplot) then cplot = cplot[0] else cplot = 0

	if z1 eq 0 then z1 = findgen(501)/50 
	if z1[0] eq 0 then cplot = 1 		;turn on plotting and show cosmological distances as a function of z
	
	Dh = 3000./h
	a = findgen(2001)
	b = (z2-z1)/2000.*a+z1
	dummy = fltarr(n_elements(b))	
	dummy2 = fltarr(n_elements(b))

	Dcnorm = qromb('Ez',0.,z1)  ;comoving distance to z1 
	Danorm = Dcnorm/(1.+z1)
	DLnorm = Dcnorm*(1+z1)
	
	Dc = Dcnorm*Dh
	Da = Danorm*Dh
	DL = DLnorm*Dh
	scale = 2*!PI*Da*10^3/360/3600	;gives angular scale for Da in kpc/"
	mpctocm = 3.08568e24

	Th = 9.78e+09/h
	TLnorm = qromb('Tz',0,z1)
	TL = TLnorm*Th
	age = Th - TL
	
	if cplot ne 0 then begin
	splot, z1, Dcnorm, xrange=[0,10], yrange=[0,max(DLnorm)], psym=10, color=2
	oplot, z1, Danorm, color=1
	oplot, z1, DLnorm, color=4
	oplot, z1, TLnorm
	endif

	if n_elements(z1) eq 1 then begin
	print, 'The lookback time at z =', z1, ' is', TL, ' years'
	print, 'The comoving distance to an object at z =', z1, ' is', Dc, ' Mpc or', Dc*mpctocm, ' cm'	
	print, 'The angular diameter distance at z =', z1, ' is', Da, ' Mpc or', Da*mpctocm, ' cm'
	print, 'This gives a angular scale at z =', z1, ' of', scale, ' kpc/"'
	print, 'The luminosity distance at z =', z1, ' is', DL, ' Mpc or', DL*mpctocm, ' cm'
	print, 'The lookback time at z =', z1, ' is', TL, ' years'
	print, 'The age of the universe at z =', z1, ' was', age, ' years'
	endif

	if (z2 ne 0) and (angle ne 0) then begin	;if an upper z and angular area are set begin the comoving volume calculation
	
	for i=0, n_elements(b)-2 do begin
		c = (b[i]+b[i+1])/2
		dummy[i] = (qromb('Ez', 0., c))^2
		dummy2[i] = qromb('Ez',b[i], b[i+1])*dummy[i]
	endfor

	integral = total(dummy2)
        Vc = Dh^3*angle*integral
	
	if n_elements(z1) eq 1 then begin
	print, 'The comoving volume between z =', z1, ' and z =', z2, ' for a survey area of', angle, ' steradians is', Vc, ' Mpc^3' 
	endif	
	endif
end
