pro HSTpostagestamps

        readcol, 'NIRSPEC_forHSTpostagestamps.dat', format='A,A,D,D,A,A', mask, slit, ra, dec, f606wfile, f814wfile

	for i=0, n_elements(mask)-1 do begin

		print, mask[i] + '.' + slit[i]
		
		;opening HST images
		f606w = '/Volumes/Data2/orelse/lemaux/HST/ACS/individual_pointings/' + f606wfile[i]
	        f814w = '/Volumes/Data2/orelse/lemaux/HST/ACS/individual_pointings/' + f814wfile[i]

		f606whdr = headfits(f606w)
		f814whdr = headfits(f814w)
	
	        f606wCRVAL1 = double(sxpar(f606whdr,'CRVAL1'))    ;initial RA for any movement along the image
	        f606wCRVAL2 = double(sxpar(f606whdr, 'CRVAL2'))   ;initial dec for any movement along the image
	
		f814wCRVAL1 = double(sxpar(f814whdr,'CRVAL1'))
		f814wCRVAL2 = double(sxpar(f814whdr, 'CRVAL2'))
	
		
       	 	f606waa = double(sxpar(f606whdr, 'CD1_1'))        ;get the transformation matrix between pixels and global coordinates 
       	 	f606wab = double(sxpar(f606whdr, 'CD1_2'))
       	 	f606wba = double(sxpar(f606whdr, 'CD2_1'))
       	 	f606wbb = double(sxpar(f606whdr, 'CD2_2'))
	
		f814waa = double(sxpar(f814whdr, 'CD1_1'))       
	        f814wab = double(sxpar(f814whdr, 'CD1_2'))
       	 	f814wba = double(sxpar(f814whdr, 'CD2_1'))
       	 	f814wbb = double(sxpar(f814whdr, 'CD2_2'))
	
	        f606wrefxpix = sxpar(f606whdr, 'CRPIX1')          ;initial x pix corresponding to CRVAL1
       	 	f606wrefypix = sxpar(f606whdr, 'CRPIX2')  
	
		f814wrefxpix = sxpar(f814whdr, 'CRPIX1')
		f814wrefypix = sxpar(f814whdr, 'CRPIX2')
	

		; make the 606W image 1st	
		cdmatrix = dblarr(2,2)
        	cdmatrix[0,0] = f606waa
        	cdmatrix[0,1] = f606wab
        	cdmatrix[1,0] = f606wba
        	cdmatrix[1,1] = f606wbb
        	cdinvmatrix = invert(cdmatrix)

        	;full spherical calculation of PA, d, d_RA, d_dec from intial point from fits header to center of slit
		p = atan( (-cos(dec[i]*2*!PI/360)*sin( (f606wCRVAL1-ra[i])*2*!PI/360 ))/(cos(f606wCRVAL2*2*!PI/360)*sin(dec[i]*2*!PI/360)-sin(f606wCRVAL2*2*!PI/360)*cos(dec[i]*2*!PI/360)*cos( 2*!PI/360*(f606wCRVAL1 - ra[i]))))
		if p lt 0 and f606wCRVAL1 gt ra[i] then p = 2*!PI + p else if p lt 0 and f606wCRVAL1 lt ra[i] then p = !PI + p else if p gt 0 and f606wCRVAL1 gt ra[i] then p = !PI + p else if p gt 0 and f606wCRVAL1 lt ra[i] then p = p

        	deltadis = 360/(2*!PI)*acos(sin(dec[i]*2*!PI/360)*sin(f606wCRVAL2*2*!PI/360)+cos(f606wCRVAL2*2*!PI/360)*cos(dec[i]*2*!PI/360)*cos(2*!PI/360*(f606wCRVAL1 - ra[i])))
        	deltaRA = deltadis*sin(p)
        	deltadec = deltadis*cos(p)

        	coormatrix =dblarr(2)
        	coormatrix[0] = deltaRA
        	coormatrix[1] = deltadec

        	;transform between difference in Ra/Dec to x/y pixels
        	newpixarray = cdinvmatrix # coormatrix
		newxcenter = floor(f606wrefxpix + newpixarray[0])
       	 	newycenter = floor(f606wrefypix + newpixarray[1])
		
		imagesize = 5.        ;distance from center of image to vertical or horizonal edge in "
        	imagepixsize = imagesize/0.03 ;imagesize rescaled by ACS for multi-drizzled binned data
        	EWedge = [newxcenter-imagepixsize, newxcenter+imagepixsize]
        	NSedge = [newycenter-imagepixsize, newycenter+imagepixsize]	
		
		;print, EWedge, NSedge	

		;remember ORIENTAT when using PA and making slits
		sxaddpar, f606whdr, 'CRVAL1', ra[i], ' R.A. (degrees) of reference pixel'
		sxaddpar, f606whdr, 'CRVAL2', dec[i], 'Declination of reference pixel'
		sxaddpar, f606whdr, 'CRPIX1', imagepixsize, 'Reference Pixel in X'
		sxaddpar, f606whdr, 'CRPIX2', imagepixsize, 'Reference Pixel in Y'
		;sxaddpar, f606whdr, 'slitPA', PA, 'Position angle of DEIMOS slit'
		;sxaddpar, f606whdr, 'slitlen', slitlength, 'Length (arcsec) of DEIMOS slit'
		sxaddpar, f606whdr, 'object', mask[i]+slit[i], 'cutout DEIMOS object'
		
		mm = mrdfits(f606w,0, /silent)
		if EWedge[0] lt 0 then begin
			finder = fltarr(EWedge[1]-EWedge[0]+2, NSedge[1]-NSedge[0]+2)
			finderpseudo = mm(0:EWedge[1], NSedge[0]:NSedge[1])
			offset = abs(EWedge[0])
			;print, size(finderpseudo, /dimensions)
			finder[floor(offset):floor(EWedge[1]+offset), 0:floor(NSedge[1]-NSedge[0])+1] = finderpseudo
		endif			 

		if EWedge[0] ge 0 then begin
		finder = mm(EWedge[0]:EWedge[1], NSedge[0]:NSedge[1])
		endif

		;if mask[i] eq 'GHF2' and LAEra[i] gt 240.996 then begin
                ;        mwrfits, finder, 'ACScutout.f606w.' + mask[i] + '.' + slit[i] +'s1.fits', f606whdr, /silent, /create
                ;endif

                ;if mask[i] eq 'GHF2' and LAEra[i] lt 240.996 then begin
                        ;mwrfits, finder, 'ACScutout.f606w.' + mask[i] + '.' + slit[i] +'s2.fits', f606whdr, /silent, /create
                ;endif

                ;if mask[i] ne 'GHF2' then begin
                        mwrfits, finder, 'ACScutout.f606w.' + mask[i] + '.' + slit[i] + '.fits', f606whdr, /silent, /create

                ;endif

	
		;814W image next

	 	cdmatrix = dblarr(2,2)
                cdmatrix[0,0] = f814waa
                cdmatrix[0,1] = f814wab
                cdmatrix[1,0] = f814wba
                cdmatrix[1,1] = f814wbb
                cdinvmatrix = invert(cdmatrix)

                ;full spherical calculation of PA, d, d_RA, d_dec from intial point from fits header to center of slit
                p = atan( (-cos(dec[i]*2*!PI/360)*sin( (f814wCRVAL1-ra[i])*2*!PI/360 ))/(cos(f814wCRVAL2*2*!PI/360)*sin(dec[i]*2*!PI/360)-sin(f814wCRVAL2*2*!PI/360)*cos(dec[i]*2*!PI/360)*cos( 2*!PI/360*(f814wCRVAL1 - ra[i]))))
		if p lt 0 and f814wCRVAL1 gt ra[i] then p = 2*!PI + p else if p lt 0 and f814wCRVAL1 lt ra[i] then p = !PI + p else if p gt 0 and f814wCRVAL1 gt ra[i] then p = !PI + p else if p gt 0 and f814wCRVAL1 lt ra[i] then p = p
	
                deltadis = 360/(2*!PI)*acos(sin(dec[i]*2*!PI/360)*sin(f814wCRVAL2*2*!PI/360)+cos(f814wCRVAL2*2*!PI/360)*cos(dec[i]*2*!PI/360)*cos(2*!PI/360*(f814wCRVAL1 - ra[i])))
                deltaRA = deltadis*sin(p)
                deltadec = deltadis*cos(p)

                coormatrix =dblarr(2)
                coormatrix[0] = deltaRA
                coormatrix[1] = deltadec

                ;transform between difference in Ra/Dec to x/y pixels
                newpixarray = cdinvmatrix # coormatrix
                newxcenter = floor(f814wrefxpix + newpixarray[0])
                newycenter = floor(f814wrefypix + newpixarray[1])
	
                EWedge = [newxcenter-imagepixsize, newxcenter+imagepixsize]
                NSedge = [newycenter-imagepixsize, newycenter+imagepixsize]

		print, newxcenter, newycenter

		mm = mrdfits(f814w,0, /silent)
	
		if EWedge[0] lt 0 then begin
                        finder = fltarr(EWedge[1]-EWedge[0]+2, NSedge[1]-NSedge[0]+2)
                        finderpseudo = mm(0:EWedge[1], NSedge[0]:NSedge[1])
                        offset = abs(EWedge[0])
                	finder[floor(offset):floor(EWedge[1]+offset), 0:floor(NSedge[1]-NSedge[0])+1] = finderpseudo
		endif

                if EWedge[0] ge 0 then begin
                finder = mm(EWedge[0]:EWedge[1], NSedge[0]:NSedge[1])
                endif
		
		sxaddpar, f814whdr, 'CRVAL1', ra[i], ' R.A. (degrees) of reference pixel'
                sxaddpar, f814whdr, 'CRVAL2', dec[i], 'Declination of reference pixel'
                sxaddpar, f814whdr, 'CRPIX1', imagepixsize, 'Reference Pixel in X'
                sxaddpar, f814whdr, 'CRPIX2', imagepixsize, 'Reference Pixel in Y'
                ;sxaddpar, f814whdr, 'slitPA', PA, 'Position angle of DEIMOS slit'
                ;sxaddpar, f814whdr, 'slitlen', slitlength, 'Length (arcsec) of DEIMOS slit'
                sxaddpar, f814whdr, 'object', mask[i]+slit[i], 'cutout DEIMOS object'                

		;if mask[i] eq 'GHF2' and LAEra[i] gt 240.996 then begin 
		;	mwrfits, finder, 'ACScutout.f814w.' + mask[i] + '.' + slit[i] +'s1.fits', f814whdr, /silent, /create
		;endif

		;if mask[i] eq 'GHF2' and LAEra[i] lt 240.996 then begin
                ;        mwrfits, finder, 'ACScutout.f814w.' + mask[i] + '.' + slit[i] +'s2.fits', f814whdr, /silent, /create
                ;endif

		;if mask[i] ne 'GHF2' then begin
			mwrfits, finder, 'ACScutout.f814w.' + mask[i] + '.' + slit[i] + '.fits', f814whdr, /silent, /create

		;endif

	endfor

end
