function Ez, x
    distance = 1/sqrt(0.3*(1+x)^3+0.7) ;concordance cosmology assumtion omega_M = 0.3, omega_lam = 0.7, to be consistant with DEEP2
   return, distance

end

pro MIPS24umfluxtolum, cat, outcat

	readcol, cat, format='D,D,D', m24, m24err, zlist
		
	openw, lun, outcat, width=300, /get_lun
        tab = '09'XB

	for i=0, n_elements(m24)-1 do begin
	
	z= zlist[i]
	
	if m24[i] ne -1 then begin
		fdens24 = (10^(-(m24[i]-4.39799)/2.5))*1e-29
		fdens24negerr = fdens24 - (10^(-(m24[i]+m24err[i]-4.39799)/2.5))*1e-29
		fdens24poserr = (10^(-(m24[i]-m24err[i]-4.39799)/2.5))*1e-29 - fdens24; - fdens24)*1e-29
	endif else begin
		fdens24 = 1e-40
		fdens24negerr = 1e-40
		fdens24poserr = 1e-40
	endelse

	freq24 = 1.25e13 ;freqency of 24 um, doing the dumbest thing, this is assuming the galaxy is at z=0, really it is rest frame 12 um for Cl1604 galaxies

	f24 = fdens24*freq24
	f24negerr = fdens24negerr*freq24
	f24poserr = fdens24poserr*freq24

	Dcnorm = qromb('Ez', 0.0, z)
        Dc = 1.323e28*Dcnorm    ;assuming h = 0.7, Dc is in [cm]
	lumdist = double((1+z)*Dc)
	lum24 = double(f24*4*!pi*lumdist^2)
	lum24poserr = double(f24poserr*4*!pi*lumdist^2)
	lum24negerr = double(f24negerr*4*!pi*lumdist^2)
	loglum24 = alog10(lum24)
	loglum24poserr = alog10(lum24poserr)
	loglum24negerr = alog10(lum24negerr)
	
	printf, lun, lum24, '   ', lum24poserr, '   ', lum24negerr, '   ', loglum24, '   ', loglum24poserr, '   ', loglum24negerr 


	endfor
	free_lun, lun

end
