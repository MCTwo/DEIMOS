pro make2dcutout
	
	readcol, '1604LYA_all_for2d.dat', format='A,A,D,D,A', mask, slit, z, ypos, Q

	for i=0, n_elements(mask)-1 do begin 
	
	slitB = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/' + mask[i] + '/*/slit.' + mask[i] + '.' + slit[i] + 'B.fits.gz'
	slitR = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/' + mask[i] + '/*/slit.' + mask[i] + '.' + slit[i] + 'R.fits.gz'

	print, mask[i] + '.' + slit[i] + '.' + Q[i]
	print, slitB

	specB = mrdfits(slitB,1, /silent)
	specR = mrdfits(slitR,1, /silent)
	
	lambdac = (1+z[i])*1215.7 

	if specB.lambda0[4095] + specB.dlambda[4095,ypos[i]] gt lambdac then begin ;check if the emission is on the red slit or the blue

		spec = specB
		print, 'Using blue slit file!'
	endif
	
	if specB.lambda0[4095] + specB.dlambda[4095,ypos[i]] lt lambdac then begin
		spec=specR
		print, 'Using red slit file!'
	endif
	
	pixc =	where(abs(spec.lambda0 + mean(spec.dlambda[*,ypos[i]]) - lambdac) eq min(abs(spec.lambda0 + mean(spec.dlambda[*,ypos[i]]) - lambdac))) ;find the pixel which closely corresponds to the observed wavlenegth

	specsize =  size(spec.flux, /dimensions)
	
	if pixc-40 lt 0. or pixc+40 gt specsize[0] then begin
		message, 'Unable to cutout spectrum, change central wavelength or change cutout size in dispersion dimnension!'
	endif

	if ypos[i]-7 lt 0. or ypos[i]+7 gt specsize[1] then begin
		message, 'Unable to cutout spectrum, change y-position or change cutout size in spatial dimnension!'
	endif
	
	cutout = spec.flux[pixc[0]-40:pixc[0]+40, ypos[i]-14:ypos[i]+14]
	lambda0 = spec.lambda0[pixc[0]-40:pixc[0]+40]
	dlambda = spec.dlambda[pixc[0]-40:pixc[0]+40, ypos[i]-14:ypos[i]+14]
	ivar = spec.ivar[pixc[0]-40:pixc[0]+40, ypos[i]-14:ypos[i]+14]
	result = {flux: cutout, ivar:ivar, lambda0:lambda0, dlambda:dlambda}

	cutoutinv = cutout*(-1)

	if mask[i] eq 'GHF2' then begin
		mwrfits, result, 'cutout.' + mask[i] + '.' + slit[i] + '.' + Q[i] + '.fits'
	endif 

	if mask[i] ne 'GHF2' then begin
		mwrfits, result, 'cutout.' + mask[i] + '.' + slit[i] + '.fits' 
	endif
	endfor

	;set_plot, 'PS'	
	;device, /color, bits_per_pixel=24, filename = 'narf.ps', xoffset=0.1, yoffset= 0.1, xsize=7.75, ysize=10.25, /inches
	;tv, bytscl(cutoutinv, min=-100, max=50), 0.5, 0.5, xsize=0.4, ysize=0.4, /normal
	;endfor
	;device, /close_file
end	
