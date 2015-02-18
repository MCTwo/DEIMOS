pro serendip_finder_subaru, field, mask, slit, offset, multiband=multiband, band=band, fieldnum=fieldnum

; Adapted from original serendip_finder.pro, BCL 12/4/10
;
; Add more detailed notes later...
; Slit must be a string of the form XXX (e.g., 002, 032, 132)
; offset should be positive or negative depending on if serendip is above or below the target. Is in DEIMOS pixels. 
; band is the band for the atv image, the other bands come up as static plots (if multiband is set). Default is r', can change if desired
; if multiband is set then will output static plot windows of all three bands in addition to the atv window 

	; first get the target info, slit PA, and slit length from the bintabs file in the mask directory

	bintab = mrdfits('/Volumes/Data2/orelse/lemaux/deimos/ORELSEmasks/'+mask+'/*/'+mask+'.bintabs.fits',1, /silent)
	slittab = mrdfits('/Volumes/Data2/orelse/lemaux/deimos/ORELSEmasks/'+mask+'/*/'+mask+'.bintabs.fits',3, /silent)
	maskinfo = mrdfits('/Volumes/Data2/orelse/lemaux/deimos/ORELSEmasks/'+mask+'/*/'+mask+'.bintabs.fits',2,/silent)

	index = where(slit eq strcompress(slittab.slitname, /remove_all) )		;find the bintab entry that corresponds to the target slit


	if index[0] lt 0 then begin
		message, 'Slit was not found in bintabs file, check that the slit exists and that it is written correctly'
	endif

	IDtar = bintab[index].object		;get information about the target, im assuming that the indexing is the same for the 1st and 3rd 
	RAtar = bintab[index].ra_obj		;extensions of the bintab file, maybe not true always in which case have to use ID
	Dectar = bintab[index].dec_obj
	ibandtar = bintab[index].mag

	print, ' '	
	print, 'The target object should be ' + strcompress(string(IDtar), /remove_all) + ' at RA = ' + strcompress(string(RAtar), /remove_all) + $
		' and Dec = ' + strcompress(string(dectar), /remove_all) + ' and i magnitude = ' + strcompress(string(ibandtar), /remove_all)

	slitlen = slittab[index].slitlen
	slitpa = slittab[index].slitlpa
	slitwidth = slittab[index].slitwid		;always gonna be 1" for our data, but generalizing	
	maskpa = maskinfo.PA_PNT

	print, 'The PA (E of N) of this slit is ' + strcompress(string(slitpa), /remove_all)
	;print, 'The PA (E of N) of this mask is ' + strcompress(string(maskpa), /remove_all)

	checkPA = abs(abs(slitpa) - abs(maskpa))         ;total slit+mask PA 
		
	if checkPA le 10. and slitPA lt 0 then begin
                print, ''
                print, '------'
                print, 'WARNING - Slit may have flipped orientation!! The slit and mask PAs are similar and the slit PA is negative'
                print, ''
                print, 'If no serendip is found try changing the sign of the offset'
                print, '------'
                print, ''
        endif

	; next read in the obj_info file to get the position of the trace from the bottom of the slit
	objinfo = mrdfits('/Volumes/Data2/orelse/lemaux/deimos/ORELSEmasks/'+mask+'/*/obj_info.'+mask+'.fits',1,/silent)

	indexobj = where(strcompress(objinfo.objno, /remove_all) eq strcompress(IDtar, /remove_all))	;not a one to one mapping between the bintabs and obj_info cuz of 
													;serendips and 2 entries for each slit in obj_info (B & R), find a match	
	objposas = objinfo[indexobj[0]].objpos * 0.1185	;gets the position of the target trace relative to the bottom of the slit in arcsec (objpos is in DEIMOS pixels). 
							;take the 0th element of indexobj to only get the blue side info, red side is essentially redundant

	objpos = objposas/0.20 -  slitlen/2./0.20		;target trace position from middle of slit 

	print, 'The slit length is ' + strcompress(string(slitlen), /remove_all) + '" and the target object is ' +  strcompress(string(objpos*0.20), /remove_all) + '" from the middle of the slit'

	; now read in the LFC image corresponding to the field of interest

	if field ne 'XL005' and field ne 'XLSS005' then begin
		filecount = findfile('/Volumes/MegaMegaData/orelse/lemaux/Subaru/*' + field + '*/*' + field + '*r.fits', count=numfields)
	
		if numfields eq 1 then begin
	
			fileR = '/Volumes/MegaMegaData/orelse/lemaux/Subaru/*' + field + '*/*' + field + '*r.fits'
			fileI = '/Volumes/MegaMegaData/orelse/lemaux/Subaru/*' + field + '*/*' + field + '*i.fits' 
			fileZ = '/Volumes/MegaMegaData/orelse/lemaux/Subaru/*' + field + '*/*' + field + '*z.fits'
			
			if n_elements(band) le 0 then begin
			
				hdr = headfits(fileR)
				image = mrdfits(fileR,0,header, /silent)
				imageR = image
	
				if n_elements(multiband) gt 0 then begin
	
		                	imageI = mrdfits(fileI,0, /silent) 
        		        	imageZ = mrdfits(fileZ,0, /silent)
				endif			

			endif else begin
			
				if band eq 'r' then hdr = headfits(fileR)
		                if band eq 'i' then hdr = headfits(fileI)
		                if band eq 'z' then hdr = headfits(fileZ)

				if band eq 'r' then image = mrdfits(fileR,0,header, /silent)	
				if band eq 'i' then image = mrdfits(fileI,0,header, /silent)
				if band eq 'z' then image = mrdfits(fileI,0,header, /silent)
			
				if n_elements(multiband) gt 0 then begin
					imageR = mrdfits(fileR,0, /silent)
       	        	       	 	imageI = mrdfits(fileI,0, /silent)
       		               	 	imageZ = mrdfits(fileZ,0, /silent)
				endif

			endelse
	

			if n_elements(image) le 0 then begin
				message, 'You must pick either r, i, or z for image band'
			endif

		endif else begin
		
			if n_elements(fieldnum) eq 0 then begin
                	        message, 'There is more than one pointing, you must specify which number field your target is in'
                	endif      

			fileR = '/Volumes/MegaMegaData/orelse/lemaux/Subaru/*' + field + '*/*' + field + '_' + fieldnum + '_r.fits'
                	fileI = '/Volumes/MegaMegaData/orelse/lemaux/Subaru/*' + field + '*/*' + field + '_' + fieldnum + '_i.fits'
                	fileZ = '/Volumes/MegaMegaData/orelse/lemaux/Subaru/*' + field + '*/*' + field + '_' + fieldnum + '_z.fits'
                
			if n_elements(band) le 0 then begin
                		hdr = headfits(fileR)
			        image = mrdfits(fileR,0,header, /silent)
                	        imageR = image
                	        if n_elements(multiband) gt 0 then begin 
			
					imageI = mrdfits(fileI,0, /silent)
                	        	imageZ = mrdfits(fileZ,0, /silent)
				endif

                	endif else begin

				if band eq 'r' then hdr = headfits(fileR)
	                	if band eq 'i' then hdr = headfits(fileI)
       		        	if band eq 'z' then hdr = headfits(fileZ)
	
	                        if band eq 'r' then image = mrdfits(fileR,0,header, /silent)
	                        if band eq 'i' then image = mrdfits(fileI,0,header, /silent)
	                        if band eq 'z' then image = mrdfits(fileI,0,header, /silent)
	
		               if n_elements(multiband) gt 0 then begin
					imageR = mrdfits(fileR,0, /silent)
	                	      	imageI = mrdfits(fileI,0, /silent)
		                        imageZ = mrdfits(fileZ,0, /silent)
				endif
	
                	endelse


	                if n_elements(image) le 0 then begin
        	                message, 'You must pick either r, i, or z for image band'
                	endif

		endelse

	endif else begin

		if field ne 'XL005' and field ne 'XLSS005' then begin
			message, 'Field not recognized, check to make sure there are associated optical images'
		endif

		fileU = '/Volumes/Data2/orelse/lemaux/CFHT/*' + field + '*/D1U.2008B.CFHTLS.cutout.fits'
                fileG = '/Volumes/Data2/orelse/lemaux/CFHT/*' + field + '*/D1G.2008B.CFHTLS.cutout.fits'
                fileR = '/Volumes/Data2/orelse/lemaux/CFHT/*' + field + '*/D1R.2008B.CFHTLS.cutout.fits'
		fileI = '/Volumes/Data2/orelse/lemaux/CFHT/*' + field + '*/D1I.2008B.CFHTLS.cutout.fits'

                if n_elements(band) le 0 then begin

                        hdr = headfits(fileR)
                        image = mrdfits(fileR,0,header, /silent)
                        imageR = image

                        if n_elements(multiband) gt 0 then begin

                                imageG = mrdfits(fileG,0, /silent)
                                imageI = mrdfits(fileI,0, /silent)
                        endif

                endif else begin

                        if band eq 'u' then hdr = headfits(fileU)
                        if band eq 'g' then hdr = headfits(fileG)
                        if band eq 'r' then hdr = headfits(fileR)
			if band eq 'i' then hdr = headfits(fileI)

			if band eq 'u' then image = mrdfits(fileU,0,header, /silent)
                        if band eq 'g' then image = mrdfits(fileG,0,header, /silent)
                        if band eq 'r' then image = mrdfits(fileR,0,header, /silent)
                        if band eq 'i' then image = mrdfits(fileI,0,header, /silent)

                        if n_elements(multiband) gt 0 then begin
                                imageU = mrdfits(fileU,0, /silent)
				imageG = mrdfits(fileG,0, /silent)
                                imageR = mrdfits(fileR,0, /silent)
                                imageI = mrdfits(fileI,0, /silent)
                        endif

                endelse


                if n_elements(image) le 0 then begin
                        message, 'This is a CFHT image must pick either u, g, r, or i for image band'
                endif
	endelse
		
	CRVAL1 = double(sxpar(hdr,'CRVAL1'))	;initial RA for any movement along the image
	CRVAL2 = double(sxpar(hdr, 'CRVAL2'))	;initial dec for any movement along the image
	
	aa = double(sxpar(hdr, 'CD1_1'))	;get the transformation matrix between pixels and global coordinates 
	ab = double(sxpar(hdr, 'CD1_2'))
	ba = double(sxpar(hdr, 'CD2_1'))
	bb = double(sxpar(hdr, 'CD2_2'))

	refxpix = sxpar(hdr, 'CRPIX1')		;initial x pix corresponding to CRVAL1
	refypix = sxpar(hdr, 'CRPIX2')		;initial y pix coresponding to CRVAL2


	; next read in full phot. catalog with ID/Ra/Dec/mags 
	if field ne 'XL005' then readcol, '/Volumes/MegaMegaData/orelse/lemaux/Subaru/phot_cats/subaru.*' + field + '*.cat', format='A,D,D,D,D,D', LFCID, LFCRA, LFCDec, LFCrband, LFCiband, LFCzband, /silent
	if field eq 'XL005' then readcol, '/Volumes/Data2/orelse/lemaux/LFC/phot_cats/lfc.*XLSS005*',format='A,D,D,D,D,D', LFCID, LFCRA, LFCDec, LFCrband, LFCiband, LFCzband, /silent	
	
	;making range of RAs and decs so don't have to calculate distances for every object, range is the length of the slit in every direction, allows for maximal possible offset, adding 3" to each side so I don't miss any serendip IDs
	
	closeobj = where(LFCRA gt RAtar - slitlen/3600. - 3./3600. and LFCRA lt RAtar + slitlen/3600. + 3./3600. and LFCDec gt Dectar-slitlen/3600. - 3./3600. and LFCDec lt Dectar+slitlen/3600. + 3./3600. and LFCID ne 'F' + strcompress(IDtar,/remove_all) and LFCID ne strcompress(IDtar,/remove_all))

	if closeobj[0] ne -1 then begin
		print, 'The number of potential serendips detected in the LFC image is ' + strcompress(string(n_elements(closeobj)), /remove_all)
		dtarLFC = fltarr(n_elements(closeobj))
	        ptarLFC = fltarr(n_elements(closeobj)) 

	        closeLFCID = LFCID[closeobj]
	        closeLFCRA = LFCRA[closeobj]
	        closeLFCDec = LFCDec[closeobj]
		
		 for i=0, n_elements(closeLFCID)-1 do begin

                ; calculate distance and PA from target to each close LFC object in arcsec & radians
        	        dtarLFC[i] = 3600.*360./(2.*!PI)*acos(sin(dectar[0]*2.*!PI/360.)*sin(closeLFCdec[i]*2.*!PI/360.)+cos(dectar[0]*2.*!PI/360.)*cos(closeLFCDec[i]*2.*!PI/360.)*cos(2.*!PI/360.*(RAtar[0]-closeLFCRA[i])))
        	        ptarLFC[i] = atan( (-cos(dectar*2*!PI/360)*sin( (RAtar-closeLFCRA[i])*2*!PI/360 ))/(cos(dectar*2*!PI/360)*sin(closeLFCDec[i]*2*!PI/360)-sin(dectar*2*!PI/360)*cos(closeLFCDec[i]*2*!PI/360)*cos( 2*!PI/360*(RAtar-closeLFCRA[i]))))

                	if ptarLFC[i] lt 0 and RAtar gt closeLFCRA[i] then ptarLFC[i] = 2.*!PI + ptarLFC[i] else $
                	if ptarLFC[i] lt 0 and RAtar lt closeLFCRA[i] then ptarLFC[i] = !PI + ptarLFC[i] else $
                	if ptarLFC[i] gt 0 and RAtar gt closeLFCRA[i] then ptarLFC[i] = !PI + ptarLFC[i] else $ 
                	if ptarLFC[i] gt 0 and RAtar lt closeLFCRA[i] then ptarLFC[i] = ptarLFC[i] 

	        endfor

        	sindex = where(dtarLFC lt 20.)          ;make a cut at 10", potentially missing some really far away serendips, but can lengthen

        ; get information on potential serendips

        	sID = LFCID[sindex]
        	sRA = LFCRA[sindex]
        	sdec = LFCDec[sindex]
        	sdtarser = dtarLFC[sindex]              ;distance from target to potential serendips
	
	endif else begin
		print, 'There are no potential serendips detected in the LFC image'
	endelse

	offsetas = offset*0.1185 		;offset from target to actual serendip in arcsec, 0.1185"/pix is the plate scale of DEIMOS
	
	;getting the new image center

	p = atan( (-cos(dectar*2*!PI/360)*sin( (CRVAL1-RAtar)*2*!PI/360 ))/(cos(CRVAL2*2*!PI/360)*sin(dectar*2*!PI/360)-sin(CRVAL2*2*!PI/360)*cos(dectar*2*!PI/360)*cos( 2*!PI/360*(CRVAL1 - RAtar))))
        if p lt 0 and CRVAL1 gt RAtar then p = 2*!PI + p else if p lt 0 and CRVAL1 lt RAtar then p = !PI + p else if p gt 0 and CRVAL1 gt RAtar then p = !PI + p else if p gt 0 and CRVAL1 lt RAtar then p = p

        deltadis = 360/(2*!PI)*acos(sin(dectar*2*!PI/360)*sin(CRVAL2*2*!PI/360)+cos(CRVAL2*2*!PI/360)*cos(dectar*2*!PI/360)*cos(2*!PI/360*(CRVAL1 - RAtar)))


        deltaRA = deltadis*sin(p)
        deltadec = deltadis*cos(p)

        coormatrix =dblarr(2)
        coormatrix[0] = deltaRA
        coormatrix[1] = deltadec

	cdmatrix = dblarr(2,2)
        cdmatrix[0,0] = aa
        cdmatrix[0,1] = ab
        cdmatrix[1,0] = ba
        cdmatrix[1,1] = bb
        cdinvmatrix = invert(cdmatrix)

        ;transform between difference in Ra/Dec to x/y pixels
        newpixarray = cdinvmatrix # coormatrix

        newxcenter = floor(refxpix + newpixarray[0])
        newycenter = floor(refypix + newpixarray[1])

        ;112x112 guide image can change value if needed, corresponds to 20"x20" 

        imagesize = 15.        ;distance from center of image to vertical or horizonal edge in "
        imagepixsize = imagesize/0.20 ;imagesize rescaled by Suprime-Cam pixel size for 1x1 binned data
        EWedge = [newxcenter-imagepixsize, newxcenter+imagepixsize]
        NSedge = [newycenter-imagepixsize, newycenter+imagepixsize]

	sxaddpar, header, 'NAXIS1', round(2*imagepixsize)
	sxaddpar, header, 'NAXIS2', round(2*imagepixsize)
        sxaddpar, header, 'CRVAL1', RAtar
        sxaddpar, header, 'CRVAL2', Dectar
        sxaddpar, header, 'CRPIX1', imagepixsize
        sxaddpar, header, 'CRPIX2', imagepixsize
	 

	atv, image(EWedge[0]:EWedge[1], NSedge[0]:NSedge[1]), header=header
	;atv_setcompass
	;atv_setscalebar
	
	;finderI = mm(EWedge[0]:EWedge[1], NSedge[0]:NSedge[1])
	;finderZ = mm(EWedge[0]:EWedge[1], NSedge[0]:NSedge[1])

	slitsize = dblarr(2)    ;2 array with slit width and slit length taken from the info in the bintabs file 
        slitsize[0] = slitwidth/0.20
        slitsize[1] = slitlen/0.20

	loadcolors
	
	; crap for creating box

	if strpos(field, '0224') lt 0 then begin		;checking if using the RCS0224 cuz that image is orinted North right, East up...

		if slitPA gt 0 then begin
			leftbottomx = imagepixsize + objpos*sin(slitPA*2.*!pi/360.) + slitsize[1]/2*sin(slitPA*2.*!pi/360.) - slitsize[0]/2*cos(slitPA*2.*!pi/360.)     ;is slitsize/2 because starting from center 
	                lefttopx = imagepixsize + objpos*sin(slitPA*2.*!pi/360.) - slitsize[1]/2*sin(slitPA*2.*!pi/360.) - slitsize[0]/2*cos(slitPA*2.*!pi/360.)
	                rightbottomx = imagepixsize + objpos*sin(slitPA*2.*!pi/360.) + slitsize[1]/2*sin(slitPA*2.*!pi/360.) + slitsize[0]/2*cos(slitPA*2.*!pi/360.)
	                righttopx = imagepixsize + objpos*sin(slitPA*2.*!pi/360.) - slitsize[1]/2*sin(slitPA*2.*!pi/360.) + slitsize[0]/2*cos(slitPA*2.*!pi/360.)
	        endif else begin
	                leftbottomx = imagepixsize + objpos*sin(slitPA*2.*!pi/360.) + slitsize[1]/2*sin(slitPA*2.*!pi/360.) - slitsize[0]/2*cos(slitPA*2.*!pi/360.)
	                lefttopx = imagepixsize + objpos*sin(slitPA*2.*!pi/360.) - slitsize[1]/2*sin(slitPA*2.*!pi/360.) - slitsize[0]/2*cos(slitPA*2.*!pi/360.)
	                rightbottomx = imagepixsize + objpos*sin(slitPA*2.*!pi/360.) + slitsize[1]/2*sin(slitPA*2.*!pi/360.) + slitsize[0]/2*cos(slitPA*2.*!pi/360.)
	                righttopx = imagepixsize + objpos*sin(slitPA*2.*!pi/360.) - slitsize[1]/2*sin(slitPA*2.*!pi/360.) + slitsize[0]/2*cos(slitPA*2.*!pi/360.)
	        endelse
	
	        if slitPA gt 0 then begin
	                leftbottomy = imagepixsize - objpos*cos(slitPA*2.*!pi/360.) - slitsize[1]/2*cos(slitPA*2.*!pi/360.) - slitsize[0]/2*sin(slitPA*2.*!pi/360.)
	                lefttopy = imagepixsize - objpos*cos(slitPA*2.*!pi/360.) + slitsize[1]/2*cos(slitPA*2.*!pi/360.) - slitsize[0]/2*sin(slitPA*2.*!pi/360.)
	                rightbottomy = imagepixsize - objpos*cos(slitPA*2.*!pi/360.) - slitsize[1]/2*cos(slitPA*2.*!pi/360.) + slitsize[0]/2*sin(slitPA*2.*!pi/360.)
	                righttopy = imagepixsize - objpos*cos(slitPA*2.*!pi/360.) + slitsize[1]/2*cos(slitPA*2.*!pi/360.) + slitsize[0]/2*sin(slitPA*2.*!pi/360.)
	        endif else begin
	                leftbottomy = imagepixsize - objpos*cos(slitPA*2.*!pi/360.) - slitsize[1]/2*cos(slitPA*2.*!pi/360.) - slitsize[0]/2*sin(slitPA*2.*!pi/360.)
	                lefttopy = imagepixsize - objpos*cos(slitPA*2.*!pi/360.) + slitsize[1]/2*cos(slitPA*2.*!pi/360.) - slitsize[0]/2*sin(slitPA*2.*!pi/360.)
	                rightbottomy = imagepixsize - objpos*cos(slitPA*2.*!pi/360.) - slitsize[1]/2*cos(slitPA*2.*!pi/360.) + slitsize[0]/2*sin(slitPA*2.*!pi/360.)
	                righttopy = imagepixsize - objpos*cos(slitPA*2.*!pi/360.) + slitsize[1]/2*cos(slitPA*2.*!pi/360.) + slitsize[0]/2*sin(slitPA*2.*!pi/360.)
        endelse


	endif else begin		;change some of the signs of the box endpoints and the x and y geometrical functions to accomidate the RCS0224 image orientation, y -> -x & x -> y
		
		if slitPA gt 0 then begin
                        leftbottomy = imagepixsize - objpos*sin(slitPA*2.*!pi/360.) - slitsize[1]/2*sin(slitPA*2.*!pi/360.) + slitsize[0]/2*cos(slitPA*2.*!pi/360.)     ;is slitsize/2 because starting from center 
                        lefttopy = imagepixsize - objpos*sin(slitPA*2.*!pi/360.) + slitsize[1]/2*sin(slitPA*2.*!pi/360.) + slitsize[0]/2*cos(slitPA*2.*!pi/360.)
                        rightbottomy = imagepixsize - objpos*sin(slitPA*2.*!pi/360.) - slitsize[1]/2*sin(slitPA*2.*!pi/360.) - slitsize[0]/2*cos(slitPA*2.*!pi/360.)
                        righttopy = imagepixsize - objpos*sin(slitPA*2.*!pi/360.) + slitsize[1]/2*sin(slitPA*2.*!pi/360.) - slitsize[0]/2*cos(slitPA*2.*!pi/360.)
                endif else begin
                        leftbottomy = imagepixsize - objpos*sin(slitPA*2.*!pi/360.) - slitsize[1]/2*sin(slitPA*2.*!pi/360.) + slitsize[0]/2*cos(slitPA*2.*!pi/360.)
                        lefttopy = imagepixsize - objpos*sin(slitPA*2.*!pi/360.) + slitsize[1]/2*sin(slitPA*2.*!pi/360.) + slitsize[0]/2*cos(slitPA*2.*!pi/360.)
                        rightbottomy = imagepixsize - objpos*sin(slitPA*2.*!pi/360.) - slitsize[1]/2*sin(slitPA*2.*!pi/360.) - slitsize[0]/2*cos(slitPA*2.*!pi/360.)
                        righttopy = imagepixsize - objpos*sin(slitPA*2.*!pi/360.) + slitsize[1]/2*sin(slitPA*2.*!pi/360.) - slitsize[0]/2*cos(slitPA*2.*!pi/360.)
                endelse

                if slitPA gt 0 then begin
                        leftbottomx = imagepixsize - objpos*cos(slitPA*2.*!pi/360.) - slitsize[1]/2*cos(slitPA*2.*!pi/360.) - slitsize[0]/2*sin(slitPA*2.*!pi/360.)
                        lefttopx = imagepixsize - objpos*cos(slitPA*2.*!pi/360.) + slitsize[1]/2*cos(slitPA*2.*!pi/360.) - slitsize[0]/2*sin(slitPA*2.*!pi/360.)
                        rightbottomx = imagepixsize - objpos*cos(slitPA*2.*!pi/360.) - slitsize[1]/2*cos(slitPA*2.*!pi/360.) + slitsize[0]/2*sin(slitPA*2.*!pi/360.)
                        righttopx = imagepixsize - objpos*cos(slitPA*2.*!pi/360.) + slitsize[1]/2*cos(slitPA*2.*!pi/360.) + slitsize[0]/2*sin(slitPA*2.*!pi/360.)
                endif else begin
                        leftbottomx = imagepixsize - objpos*cos(slitPA*2.*!pi/360.) - slitsize[1]/2*cos(slitPA*2.*!pi/360.) - slitsize[0]/2*sin(slitPA*2.*!pi/360.)
                        lefttopx = imagepixsize - objpos*cos(slitPA*2.*!pi/360.) + slitsize[1]/2*cos(slitPA*2.*!pi/360.) - slitsize[0]/2*sin(slitPA*2.*!pi/360.)
                        rightbottomx = imagepixsize - objpos*cos(slitPA*2.*!pi/360.) - slitsize[1]/2*cos(slitPA*2.*!pi/360.) + slitsize[0]/2*sin(slitPA*2.*!pi/360.)
                        righttopx = imagepixsize - objpos*cos(slitPA*2.*!pi/360.) + slitsize[1]/2*cos(slitPA*2.*!pi/360.) + slitsize[0]/2*sin(slitPA*2.*!pi/360.)
                endelse
	endelse
		

	xleft = [leftbottomx, lefttopx]
	xright = [rightbottomx, righttopx]
	xtop = [lefttopx, righttopx]
        xbottom = [leftbottomx, rightbottomx]
	
	yleft = [leftbottomy, lefttopy]
        yright = [rightbottomy, righttopy]
        ytop = [lefttopy, righttopy]
        ybottom = [leftbottomy, rightbottomy]

	; crap for creating target circle
	
	xcirc = findgen(500)-250+imagepixsize
	ycircpos = sqrt(64-(xcirc-imagepixsize)^2) + imagepixsize	;circle of radius 8 pix
	ycircneg = -1.*sqrt(64-(xcirc-imagepixsize)^2) + imagepixsize

	atvplot, xcirc, ycircpos, color=4, THICK=2
	atvplot, xcirc, ycircneg, color=4, THICK=2

	atvplot, xleft, yleft, color=6, thick=2
	atvplot, xright, yright, color=6, thick=2
	atvplot, xtop, ytop, color=6, thick=2
	atvplot, xbottom, ybottom, color=6, thick=2


	if strpos(field, '0224') lt 0 then begin	; have to check cuz 0224 image is annoyingly rotated 90 degrees
		sercenter = [imagepixsize - offsetas/0.20*sin(slitPA*2.*!pi/360.), imagepixsize + offsetas/0.20*cos(slitPA*2.*!pi/360.)]
	endif else begin
		sercenter = [imagepixsize + offsetas/0.20*cos(slitPA*2.*!pi/360.), imagepixsize + offsetas/0.20*sin(slitPA*2.*!pi/360.)]	;switching signs and coordinates to match 0224 orientation
	endelse 

	;making an X of 10 pixels on the serendip location
	Xser_xtopleftbtmright = [sercenter[0]-5., sercenter[0]+5.]
	Xser_ytopleftbtmright = [sercenter[1]+5., sercenter[1]-5.]
	Xser_xbtmlefttopright = [sercenter[0]-5., sercenter[0]+5.]
	Xser_ybtmlefttopright = [sercenter[1]-5., sercenter[1]+5.]

	atvplot, Xser_xtopleftbtmright, Xser_ytopleftbtmright, thick=4, color=1
        atvplot, Xser_xbtmlefttopright, Xser_ybtmlefttopright, thick=4, color=1

	; now adding the LFC IDs of the nearby objects (~within 7" of the target) if there are potential serendips

	if closeobj[0] ne -1 then begin	
		dtarLFCRA = -1.*dtarLFC*sin(ptarLFC)/0.20		;get the offset from the target to the closeby objects from the previous calculation, both RA & dec output in pixels, RA is negative because of of the orientation 
	dtarLFCdec = dtarLFC*cos(ptarLFC)/0.20	

		if strpos(field, '0224') lt 0 then begin	;RCS0224 image is rotated 90 degrees from sky orintation, compensating if using this field
		
		for i=0, n_elements(dtarLFCRA)-1 do begin
			atvxyouts, imagepixsize+dtarLFCRA[i], imagepixsize+dtarLFCdec[i], closeLFCID[i], color=5, charsize=1.9, charthick=22, thick=15	;plot LFC IDs of potential serenedips in atv window
		endfor	
		
		endif else begin
		
		for i=0, n_elements(dtarLFCRA)-1 do begin
                        atvxyouts, imagepixsize+dtarLFCdec[i], imagepixsize-dtarLFCRA[i], closeLFCID[i], color=5, charsize=1.9, charthick=22, thick=15  ;North is right and East is up in this image, so flipped directions
                endfor

		endelse
	endif	

	; put title on the top 
	
	atvxyouts, -20, 2*imagepixsize+5, 'Finder for ' + strcompress(mask,/remove_all) + '.' + strcompress(slit, /remove_all) + ', target is ' + strcompress(IDtar,/remove_all), color=7, charsize=1.8, charthick=10
	
	;done with ATV plotting, now make 3 plot windows with the other bands displayed (if multiband is set)

	if n_elements(multiband) gt 0 then begin 
	
		sersize = [1., 10.]
		window, xsize=(imagepixsize*2), ysize=(imagepixsize*2), title = 'r band'
		device, decomposed=1
		tvscl, -1.*imageR(EWedge[0]:EWedge[1], NSedge[0]:NSedge[1])
		tvcircle, 8., imagepixsize, imagepixsize, color='FFFF00'XL, THICK=3
		
		if slitpa gt 0 then begin
                
                        tvbox, slitsize, imagepixsize + objpos*sin(slitPA*2.*!pi/360.), imagepixsize - objpos*cos(slitPA*2.*!pi/360.), color='FF0000'XL, angle=-slitPA, thick=2	;angle is -PA because angle is defined as clockwise and PA is counterclockwise
                endif else begin

                        tvbox, slitsize, imagepixsize + objpos*sin(slitPA*2.*!pi/360.), imagepixsize - objpos*cos(slitPA*2.*!pi/360.), color='FF0000'XL, angle=-slitPA, thick=2

                endelse
		
		tvbox, sersize, sercenter[0], sercenter[1], color='0000FF'XL, angle=45, thick=2
		tvbox, sersize, sercenter[0], sercenter[1], color='0000FF'XL, angle=-45, thick=2
	
		window, 2, xsize=(imagepixsize*2), ysize=(imagepixsize*2), Title='i band'
	        device, decomposed=1
		tvscl, -1.*imageI(EWedge[0]:EWedge[1], NSedge[0]:NSedge[1])
	        tvcircle, 8., imagepixsize, imagepixsize, color='FFFF00'XL, THICK=3
	
		if slitpa gt 0 then begin
	        
			tvbox, slitsize, imagepixsize + objpos*sin(slitPA*2.*!pi/360.), imagepixsize - objpos*cos(slitPA*2.*!pi/360.), color='FF0000'XL, angle=-slitPA, thick=2
	       	
		endif else begin
		
			tvbox, slitsize, imagepixsize + objpos*sin(slitPA*2.*!pi/360.), imagepixsize - objpos*cos(slitPA*2.*!pi/360.), color='FF0000'XL, angle=-slitPA, thick=2 
	
		endelse
	
		tvbox, sersize, sercenter[0], sercenter[1], color='0000FF'XL, angle=45, thick=2
	        tvbox, sersize, sercenter[0], sercenter[1], color='0000FF'XL, angle=-45, thick=2

		window, 1, xsize=(imagepixsize*2), ysize=(imagepixsize*2), Title='z band'
	        device, decomposed=1
		if field ne 'XL005' then tvscl, -1.*imageZ(EWedge[0]:EWedge[1], NSedge[0]:NSedge[1])
		if field eq 'XL005' or field eq 'XLSS005' then tvscl, -1.*imageG(EWedge[0]:EWedge[1], NSedge[0]:NSedge[1])
	        tvcircle, 8., imagepixsize, imagepixsize, color='FFFF00'XL, THICK=3
		
		if slitpa gt 0 then begin
                
                        tvbox, slitsize, imagepixsize + objpos*sin(slitPA*2.*!pi/360.), imagepixsize - objpos*cos(slitPA*2.*!pi/360.), color='FF0000'XL, angle=-slitPA, thick=2

                endif else begin

                        tvbox, slitsize, imagepixsize + objpos*sin(slitPA*2.*!pi/360.), imagepixsize - objpos*cos(slitPA*2.*!pi/360.), color='FF0000'XL, angle=-slitPA, thick=2

                endelse        
	
		tvbox, sersize, sercenter[0], sercenter[1], color='0000FF'XL, angle=45, thick=2
	        tvbox, sersize, sercenter[0], sercenter[1], color='0000FF'XL, angle=-45, thick=2
	
		device, decomposed=0
	endif

end
