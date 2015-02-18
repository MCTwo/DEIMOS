
PRO add_serendip, bs, serendip_info
	muy = (size(bs.flux,/dim))[1]
	sigs = randomu(seed, 3)*8
	sigx = sigs[0]+3>4 & sigy = sigs[1]+3>4
	a = sigs[2]*10
	print,'Generated emission: sigx, sigy, a: ', sigx,sigy,a

	nx = 51
	ny = muy

	x = (dindgen(nx)-nx/2)
	y = (dindgen(ny)-ny/2)+muy/4.
	yy = meshgrid(x, y, xx)

	patch = exp(-xx^2/(2*sigx^2)-yy^2/(2*sigy^2))
	patch = patch * (randomu(seed,1))[0]*a
	add_poisson, patch

	c0 = 500-25
	c1 = 500+25
	bs.flux[c0:c1,*] = bs.flux[c0:c1,*] + patch
	bs.ivar[c0:c1,*] = bs.ivar[c0:c1,*] + 1./(patch+50)

	serendip_info = {x: 500, y: muy/4, sigx: sigx, sigy: sigy, $
		a: a}
END

PRO generate_files_for_sex, path, maskname, slitname
	
	mn = strcompress(maskname, /rem)
	sn = string(slitname, format='(i3.3)')
	bs = mrdfits(path+'slit.' + mn + '.' + sn + 'B.fits.gz', 1, /silent)
	rs = mrdfits(path+'slit.' + mn + '.' + sn + 'R.fits.gz', 1, /silent)

	add_serendip, bs, serendip_info


	bl = bs.dlambda
	l0 = bs.lambda0
	for i = 0, (size(bl, /dim))[1] - 1 do begin
		bl[*,i] = bl[*,i] + l0
	endfor

	rl = rs.dlambda
	l0 = rs.lambda0
	for i = 0, (size(rl, /dim))[1] - 1 do begin
		rl[*,i] = rl[*,i] + l0
	endfor
	id = 'slit.' + mn + '.' + sn
	mwrfits, bs.flux, '/Volumes/Data/lemaux/deimos/1dfits/' + id + 'B.flux.fits', /create
	mwrfits, rs.flux, '/Volumes/Data/lemaux/deimos/1dfits/' + id + 'R.flux.fits', /create
	mwrfits, bl, '/Volumes/Data/lemaux/deimos/1dfits/' + id + 'B.lambda.fits', /create
	mwrfits, rl, '/Volumes/Data/lemaux/deimos/1dfits/' + id + 'R.lambda.fits', /create
	mwrfits, bs.ivar, '/Volumes/Data/lemaux/deimos/1dfits/' + id + 'B.ivar.fits', /create
	mwrfits, rs.ivar, '/Volumes/Data/lemaux/deimos/1dfits/' + id + 'R.ivar.fits', /create
	mwrfits, serendip_info, '/Volumes/Data/lemaux/deimos/1dfits/' + id + '.info.fits', /create
END
