pro LFCpostagestamps_forNIRSPEC
	
	readcol, 'NIRSPEC_forLFCpostagestamps.dat', format='A,A,D,D,A,A,A', mask, slit, ra, dec, rfile, ifile, zfile

	for i=0, n_elements(mask)-1 do begin

		print, mask[i] + '.' + slit[i]
		; getting slit PA and length from the bintabs file
		
		
	
		;opening LFC photometry files	
		rprime = '/Volumes/Data2/orelse/lemaux/LFC/sc1604/' + rfile[i]
		iprime = '/Volumes/Data2/orelse/lemaux/LFC/sc1604/' + ifile[i]
		zprime = '/Volumes/Data2/orelse/lemaux/LFC/sc1604/' + zfile[i]

		rhdr = headfits(rprime)
		ihdr = headfits(iprime)
		zhdr = headfits(zprime)
	
	        rCRVAL1 = double(sxpar(rhdr,'CRVAL1'))    ;initial RA for any movement along the image
	        rCRVAL2 = double(sxpar(rhdr, 'CRVAL2'))   ;initial dec for any movement along the image
	
		iCRVAL1 = double(sxpar(ihdr,'CRVAL1'))
		iCRVAL2 = double(sxpar(ihdr, 'CRVAL2'))
	
		zCRVAL1 = double(sxpar(zhdr,'CRVAL1'))
		zCRVAL2 = double(sxpar(zhdr,'CRVAL2'))
		
       	 	raa = double(sxpar(rhdr, 'CD1_1'))        ;get the transformation matrix between pixels and global coordinates 
       	 	rab = double(sxpar(rhdr, 'CD1_2'))
       	 	rba = double(sxpar(rhdr, 'CD2_1'))
       	 	rbb = double(sxpar(rhdr, 'CD2_2'))
	
		iaa = double(sxpar(ihdr, 'CD1_1'))       
	        iab = double(sxpar(ihdr, 'CD1_2'))
       	 	iba = double(sxpar(ihdr, 'CD2_1'))
       	 	ibb = double(sxpar(ihdr, 'CD2_2'))
	
		zaa = double(sxpar(zhdr, 'CD1_1'))        
	        zab = double(sxpar(zhdr, 'CD1_2'))
	        zba = double(sxpar(zhdr, 'CD2_1'))
	        zbb = double(sxpar(zhdr, 'CD2_2'))
		
	        rrefxpix = sxpar(rhdr, 'CRPIX1')          ;initial x pix corresponding to CRVAL1
       	 	rrefypix = sxpar(rhdr, 'CRPIX2')  
	
		irefxpix = sxpar(ihdr, 'CRPIX1')
		irefypix = sxpar(ihdr, 'CRPIX2')
	
		zrefxpix = sxpar(zhdr, 'CRPIX1')
       	 	zrefypix = sxpar(zhdr, 'CRPIX2')

		; make the r band image 1st	
		cdmatrix = dblarr(2,2)
        	cdmatrix[0,0] = raa
        	cdmatrix[0,1] = rab
        	cdmatrix[1,0] = rba
        	cdmatrix[1,1] = rbb
        	cdinvmatrix = invert(cdmatrix)

        	;full spherical calculation of PA, d, d_RA, d_dec from intial point from fits header to center of slit
		p = atan( (-cos(dec[i]*2*!PI/360)*sin( (rCRVAL1-ra[i])*2*!PI/360 ))/(cos(rCRVAL2*2*!PI/360)*sin(dec[i]*2*!PI/360)-sin(rCRVAL2*2*!PI/360)*cos(dec[i]*2*!PI/360)*cos( 2*!PI/360*(rCRVAL1 - ra[i]))))
        	if p lt 0 and rCRVAL1 gt ra[i] then p = 2*!PI + p else if p lt 0 and rCRVAL1 lt ra[i] then p = !PI + p else if p gt 0 and rCRVAL1 gt ra[i] then p = !PI + p else if p gt 0 and rCRVAL1 lt ra[i] then p = p

        	deltadis = 360/(2*!PI)*acos(sin(dec[i]*2*!PI/360)*sin(rCRVAL2*2*!PI/360)+cos(rCRVAL2*2*!PI/360)*cos(dec[i]*2*!PI/360)*cos(2*!PI/360*(rCRVAL1 - ra[i])))
        	deltaRA = deltadis*sin(p)
        	deltadec = deltadis*cos(p)

        	coormatrix =dblarr(2)
        	coormatrix[0] = deltaRA
        	coormatrix[1] = deltadec

        	;transform between difference in Ra/Dec to x/y pixels
        	newpixarray = cdinvmatrix # coormatrix
		newxcenter = floor(rrefxpix + newpixarray[0])
       	 	newycenter = floor(rrefypix + newpixarray[1])
		
		imagesize = 10.        ;distance from center of image to vertical or horizonal edge in "
        	imagepixsize = imagesize/0.18 ;imagesize rescaled by LFC pixel size for 1x1 binned data
        	EWedge = [newxcenter-imagepixsize, newxcenter+imagepixsize]
        	NSedge = [newycenter-imagepixsize, newycenter+imagepixsize]	

		;remember ORIENTA when using PA and making slits
		sxaddpar, rhdr, 'CRVAL1', ra[i], ' R.A. (degrees) of reference pixel'
                sxaddpar, rhdr, 'CRVAL2', dec[i], 'Declination of reference pixel'
                sxaddpar, rhdr, 'CRPIX1', imagepixsize, 'Reference Pixel in X'
                sxaddpar, rhdr, 'CRPIX2', imagepixsize, 'Reference Pixel in Y'
                sxaddpar, rhdr, 'object', mask[i]+slit[i], 'cutout NIRSPEC object'	

		mm = mrdfits(rprime,0, /silent)
		finder = mm(EWedge[0]:EWedge[1], NSedge[0]:NSedge[1])
        	mwrfits, finder, 'LFCcutout.r.' + mask[i] + '.' + slit[i] +'.fits', rhdr, /silent, /create
	
		;iband image next
		
		cdmatrix = dblarr(2,2)
                cdmatrix[0,0] = iaa
                cdmatrix[0,1] = iab
                cdmatrix[1,0] = iba
                cdmatrix[1,1] = ibb
                cdinvmatrix = invert(cdmatrix)

                p = atan( (-cos(dec[i]*2*!PI/360)*sin( (iCRVAL1-ra[i])*2*!PI/360 ))/(cos(iCRVAL2*2*!PI/360)*sin(dec[i]*2*!PI/360)-sin(iCRVAL2*2*!PI/360)*cos(dec[i]*2*!PI/360)*cos( 2*!PI/360*(iCRVAL1 - ra[i]))))
                if p lt 0 and iCRVAL1 gt ra[i] then p = 2*!PI + p else if p lt 0 and iCRVAL1 lt ra[i] then p = !PI + p else if p gt 0 and iCRVAL1 gt ra[i] then p = !PI + p else if p gt 0 and iCRVAL1 lt ra[i] then p = p

                deltadis = 360/(2*!PI)*acos(sin(dec[i]*2*!PI/360)*sin(iCRVAL2*2*!PI/360)+cos(iCRVAL2*2*!PI/360)*cos(dec[i]*2*!PI/360)*cos(2*!PI/360*(iCRVAL1 - ra[i])))
                deltaRA = deltadis*sin(p)
                deltadec = deltadis*cos(p)

                coormatrix =dblarr(2)
                coormatrix[0] = deltaRA
                coormatrix[1] = deltadec

                newpixarray = cdinvmatrix # coormatrix
                newxcenter = floor(irefxpix + newpixarray[0])
                newycenter = floor(irefypix + newpixarray[1])
		
                EWedge = [newxcenter-imagepixsize, newxcenter+imagepixsize]
                NSedge = [newycenter-imagepixsize, newycenter+imagepixsize]

		;remember ORIENTA when using PA and making slits
                sxaddpar, ihdr, 'CRVAL1', ra[i], ' R.A. (degrees) of reference pixel'
                sxaddpar, ihdr, 'CRVAL2', dec[i], 'Declination of reference pixel'
                sxaddpar, ihdr, 'CRPIX1', imagepixsize, 'Reference Pixel in X'
                sxaddpar, ihdr, 'CRPIX2', imagepixsize, 'Reference Pixel in Y'
                sxaddpar, ihdr, 'object', mask[i]+slit[i], 'cutout NIRSPEC object'	
	
		mm = mrdfits(iprime,0, /silent)
                finder = mm(EWedge[0]:EWedge[1], NSedge[0]:NSedge[1])
                mwrfits, finder, 'LFCcutout.i.' + mask[i] + '.' + slit[i] +'.fits', ihdr, /silent, /create

		;finally z band image 
		cdmatrix = dblarr(2,2)
                cdmatrix[0,0] = zaa
                cdmatrix[0,1] = zab
                cdmatrix[1,0] = zba
                cdmatrix[1,1] = zbb
                cdinvmatrix = invert(cdmatrix)

                p = atan( (-cos(dec[i]*2*!PI/360)*sin( (zCRVAL1-ra[i])*2*!PI/360 ))/(cos(zCRVAL2*2*!PI/360)*sin(dec[i]*2*!PI/360)-sin(zCRVAL2*2*!PI/360)*cos(dec[i]*2*!PI/360)*cos( 2*!PI/360*(zCRVAL1 - ra[i]))))
                if p lt 0 and zCRVAL1 gt ra[i] then p = 2*!PI + p else if p lt 0 and zCRVAL1 lt ra[i] then p = !PI + p else if p gt 0 and zCRVAL1 gt ra[i] then p = !PI + p else if p gt 0 and zCRVAL1 lt ra[i] then p = p

                deltadis = 360/(2*!PI)*acos(sin(dec[i]*2*!PI/360)*sin(zCRVAL2*2*!PI/360)+cos(zCRVAL2*2*!PI/360)*cos(dec[i]*2*!PI/360)*cos(2*!PI/360*(zCRVAL1 - ra[i])))
                deltaRA = deltadis*sin(p)
                deltadec = deltadis*cos(p)

                coormatrix =dblarr(2)
                coormatrix[0] = deltaRA
                coormatrix[1] = deltadec

                newpixarray = cdinvmatrix # coormatrix
                newxcenter = floor(zrefxpix + newpixarray[0])
                newycenter = floor(zrefypix + newpixarray[1])

                EWedge = [newxcenter-imagepixsize, newxcenter+imagepixsize]
                NSedge = [newycenter-imagepixsize, newycenter+imagepixsize]

		sxaddpar, zhdr, 'CRVAL1', ra[i], ' R.A. (degrees) of reference pixel'
                sxaddpar, zhdr, 'CRVAL2', dec[i], 'Declination of reference pixel'
                sxaddpar, zhdr, 'CRPIX1', imagepixsize, 'Reference Pixel in X'
                sxaddpar, zhdr, 'CRPIX2', imagepixsize, 'Reference Pixel in Y'
                sxaddpar, zhdr, 'object', mask[i]+slit[i], 'cutout NIRSPEC object'

		mm = mrdfits(zprime,0, /silent)
                finder = mm(EWedge[0]:EWedge[1], NSedge[0]:NSedge[1])
		mwrfits, finder, 'LFCcutout.z.' + mask[i] + '.' + slit[i] +'.fits', zhdr, /silent, /create
	endfor

end
