pro LYAmosaic
	
	set_plot, 'PS'
        loadcolors
	device, /color, xoffset=0.1, yoffset=0.2, xsize=10.7, ysize=14, /inches, filename='LYAspectralmosaic1.ps'
	
	readcol, '1604LYA_all_formosaic.dat', format='A,A,D,D,A', mask, slit, z, q, spec, /silent
	readcol, 'LAE_wtargets_all_RAdec_forHST.dat', format='A,A,D,D,A,A,D,D,D,D', maskagain, slitagain, HSTra, HSTdec, f606wfile, f814wfile, HSTLAEra, HSTLAEdec, tarypos, yslit, /silent
	readcol, 'LAE_wtargets_all_RAdec_forLFC.dat', format='A,A,D,D,A,A,A,D,D', maskagain, slitagain, LFCra, LFCdec, rfile, ifile, zfile, LFCLAEra, LFCLAEdec, LFCtarypos, LFCyslit, /silent
	
	; doing the first 6, 2nd 6, and then the last 5 on a seperate plot

	twodxsize = 0.11
	twodysize = 0.07
	postagexsize = 0.155
	postageysize = 0.11625
	startHST606 = 0.59
	startHST814 = 0.78
	startr = 0.52
        starti = 0.685
        startz = 0.85
	start2d = 0.40
	start1d = 0.18
	end1d = 0.38	

	; first the 1d spectra 
	
	for i=0, 6-1 do begin 
		
	spec[i] = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/' + spec[i]
	corrspec = fill_gap(spec[i],/tweak,/telluric,/silent,header=header)
	wave = corrspec.lambda
        airtovac, wave
        corrspec.lambda = wave
	
	range = where(corrspec.lambda ge (1+z[i])*1215.7-20 and corrspec.lambda le (1+z[i])*1215.7+20)

	constfluxcorr = double(2e-08/(3600.*!pi*(449.)^2*0.33*3e18*1e-29)) ;3e+18 is the speed of light in Ang/s and 10e-29 ergs/s/cm^2/Hz = 1 uJy

	lambda = corrspec.lambda[range]
	
	uJyspec = dblarr(n_elements(corrspec.spec[range]))	
	
		for j=0, n_elements(range)-1 do begin	
		uJyspec[j] = constfluxcorr*corrspec.spec[range[j]]*lambda[[j]];spectrum in microJanskys, converted from ergs/s/cm^2/Angstom
		endfor
		
	ypos = (5-i)/6.2 + 0.1
	
	if slit[i] eq '61s2' then begin
	xyouts, 0.02, ypos, mask[i] + '.' + '61s3', charsize=1.5, charthick=6, /normal
	endif
	
	if slit[i] eq '61s1' then begin 
	xyouts, 0.02, ypos, mask[i] + '.' + '61s2', charsize=1.5, charthick=6, /normal
	endif
	
	if slit[i] ne '61s1' and slit[i] ne '61s2' then begin
	xyouts, 0.02, ypos, mask[i] + '.' + slit[i], charsize=1.5, charthick=6, /normal
	endif

	if i eq 5 then begin 

	plot, corrspec.lambda[range], uJyspec, yrange=[-1, max(uJyspec)], ystyle=1, position=[start1d, ypos-0.057, end1d, ypos+0.057],$
        color=0, xrange = [min(lambda), max(lambda)], xstyle=1, xtitle=textoidl('\lambda (\AA)') , ytickformat = '(f5.1)', $
        thick=6, charthick=4, charsize=1.2, xthick=6, ythick=6, xticks=2, yticks=2, xminor=5, yminor=5, xticklen=0.11, yticklen=0.11, /noerase
	endif

	;xyouts, 0.18, 0.5, textoidl('f_{\nu} \muJy'), charthick=4, charsize=1.2, /normal
	
	if i ne 5 then begin
	plot, corrspec.lambda[range], uJyspec, yrange=[-1, max(uJyspec)], ystyle=1, position=[start1d, ypos-0.057, end1d, ypos+0.057],$
        color=0, xrange = [min(lambda), max(lambda)], xstyle=1, ytickformat = '(f5.1)', $ 
        thick=6, charthick=4, charsize=1.2, xthick=6, ythick=6, xticks=2, yticks=2, xminor=5, yminor=5, xticklen=0.11, yticklen=0.11, /noerase
	endif


	twodcutout = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/1604_LYA/mosaic/2dcutouts/cutout.' + mask[i] + '.0' + slit[i] + '.fits'
	twodspec = mrdfits(twodcutout,1, /silent)
	cutoutinv = twodspec.flux*(-1)
	tv, bytscl(cutoutinv, min=-13.51, max=20.5), start2d, ypos - 0.04, xsize=twodxsize, ysize=twodysize, /normal 
	
	f606w = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/1604_LYA/mosaic/postage_stamps/ACScutout.f606w.' + mask[i] + '.0' + slit[i] + '.fits'
	f814w = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/1604_LYA/mosaic/postage_stamps/ACScutout.f814w.' + mask[i] + '.0' + slit[i] + '.fits'

	image606 = mrdfits(f606w,0, /silent)*(-1)
	image814 = mrdfits(f814w,0, /silent)*(-1)

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

	cdmatrix606 = dblarr(2,2)
        cdmatrix606[0,0] = f606waa
        cdmatrix606[0,1] = f606wab
        cdmatrix606[1,0] = f606wba
        cdmatrix606[1,1] = f606wbb
        cdinvmatrix606 = invert(cdmatrix606)
	
	cdmatix814 = dblarr(2,2)
        cdmatix814[0,0] = f814waa
        cdmatix814[0,1] = f814wab
        cdmatix814[1,0] = f814wba
        cdmatix814[1,1] = f814wbb
        cdinvmatrix = invert(cdmatix814)

	slitPA = double(sxpar(f606whdr,'slitPA'))
	slitlen = double(sxpar(f606whdr,'SLITLEN'))

	orientat606 = double(sxpar(f606whdr,'ORIENTAT'))
	orientat814 = double(sxpar(f814whdr,'ORIENTAT'))

	; LAE position in 606 image	
	p = atan( (-cos(HSTLAEdec[i]*2*!PI/360)*sin( (f606wCRVAL1-HSTLAEra[i])*2*!PI/360 ))/(cos(f606wCRVAL2*2*!PI/360)*sin(HSTLAEdec[i]*2*!PI/360)-sin(f606wCRVAL2*2*!PI/360)*cos(HSTLAEdec[i]*2*!PI/360)*cos( 2*!PI/360*(f606wCRVAL1 - HSTLAEra[i]))))
        if p lt 0 and f606wCRVAL1 gt HSTLAEra[i] then p = 2*!PI + p else if p lt 0 and f606wCRVAL1 lt HSTLAEra[i] then p = !PI + p else if p gt 0 and f606wCRVAL1 gt HSTLAEra[i] then p = !PI + p else if p gt 0 and f606wCRVAL1 lt HSTLAEra[i] then p = p

        deltadis = 360/(2*!PI)*acos(sin(HSTLAEdec[i]*2*!PI/360)*sin(f606wCRVAL2*2*!PI/360)+cos(f606wCRVAL2*2*!PI/360)*cos(HSTLAEdec[i]*2*!PI/360)*cos(2*!PI/360*(f606wCRVAL1 - HSTLAEra[i])))
        deltaRA = deltadis*sin(p)
        deltadec = deltadis*cos(p)

        coormatrix =dblarr(2)
        coormatrix[0] = deltaRA
        coormatrix[1] = deltadec

        ;transform between difference in Ra/Dec to x/y pixels
        newpixarray = cdinvmatrix # coormatrix
        LAExcenter606 = floor(f606wrefxpix + newpixarray[0])
        LAEycenter606 = floor(f606wrefypix + newpixarray[1])


	; LAE position in 814 image     
        p = atan( (-cos(HSTLAEdec[i]*2*!PI/360)*sin( (f814wCRVAL1-HSTLAEra[i])*2*!PI/360 ))/(cos(f814wCRVAL2*2*!PI/360)*sin(HSTLAEdec[i]*2*!PI/360)-sin(f814wCRVAL2*2*!PI/360)*cos(HSTLAEdec[i]*2*!PI/360)*cos( 2*!PI/360*(f814wCRVAL1 - HSTLAEra[i]))))
        if p lt 0 and f814wCRVAL1 gt HSTLAEra[i] then p = 2*!PI + p else if p lt 0 and f814wCRVAL1 lt HSTLAEra[i] then p = !PI + p else if p gt 0 and f814wCRVAL1 gt HSTLAEra[i] then p = !PI + p else if p gt 0 and f814wCRVAL1 lt HSTLAEra[i] then p = p

        deltadis = 360/(2*!PI)*acos(sin(HSTLAEdec[i]*2*!PI/360)*sin(f814wCRVAL2*2*!PI/360)+cos(f814wCRVAL2*2*!PI/360)*cos(HSTLAEdec[i]*2*!PI/360)*cos(2*!PI/360*(f814wCRVAL1 - HSTLAEra[i])))
        deltaRA = deltadis*sin(p)
        deltadec = deltadis*cos(p)

        coormatrix =dblarr(2)
        coormatrix[0] = deltaRA
        coormatrix[1] = deltadec

        ;transform between difference in Ra/Dec to x/y pixels
        newpixarray = cdinvmatrix # coormatrix
        LAExcenter814 = floor(f814wrefxpix + newpixarray[0])
        LAEycenter814 = floor(f814wrefypix + newpixarray[1])		

	;tweaking some cuz there is a discrepancy between the measured slit lengths in zspec and the lengths from bintab files
	if mask[i] + '.' + slit[i] eq 'GHF2.61s1' then begin
	LAEycenter606 = LAEycenter606-20
	LAEycenter814 = LAEycenter814-20
	endif

	if mask[i] eq 'FG1' then begin
        LAEycenter606 = LAEycenter606+30
        LAEycenter814 = LAEycenter814+30
	LAExcenter606 = LAExcenter606+10
	LAExcenter814 = LAExcenter814+10
        endif

	if mask[i] + '.' + slit[i] eq 'SC1NM1.85' then begin
	LAExcenter606 = LAExcenter606-15
        LAExcenter814 = LAExcenter814-15
	endif

	print, mask[i] + '.' + slit[i]	
	;print, newpixarray[0], newpixarray[1]
	;print, HSTLAEra[i], HSTLAEdec[i]
	;print, LAExcenter606, LAEycenter606
	;print, orientat606
	;print, f606wCRVAL1, f606wCRVAL2

	windowXSize = !d.x_size
	windowYSize = !d.y_size
	x0 = (startHST606) * windowXSize
	x00 = (startHST814) * windowXSize	
	y0 = (ypos - postageysize/2 - 0.014 + 0.2/14) * windowYSize
	xsize = postagexsize * windowXSize
	ysize = postageysize * windowYSize		
	slitsize = [1/0.03*xsize/668, slitlen/0.03*xsize/668]

	taroffsetx =  (yslit[i]/2 -tarypos[i])*0.118*sin(slitPA*2*!PI/360)/0.03
	taroffsety = (yslit[i]/2 -tarypos[i])*0.118*cos(slitPA*2*!PI/360)/0.03
	
	if mask[i] ne 'FG2' then begin
		tv, bytscl(image606, min=-0.029, max=0.045), startHST606, ypos - postageysize/2, xsize=postagexsize, ysize=postageysize, /normal
		tvbox, slitsize, (334+taroffsetx)*xsize/668 + x0, (334+taroffsety)*ysize/668 + y0, color=2, angle = -slitPA, thick=4
		tvcircle, 1./0.03*xsize/668, LAExcenter606*xsize/668 + x0, LAEycenter606*ysize/668 + y0, color=6, thick=7
		tv, bytscl(image814, min=-0.0418, max=0.0551), startHST814, ypos - postageysize/2, xsize=postagexsize, ysize=postageysize, /normal
		tvbox, slitsize, (334+taroffsetx)*xsize/668 + x00, (334+taroffsety)*ysize/668 + y0, color=2, angle = -slitPA, thick=4
        	tvcircle, 1./0.03*xsize/668, LAExcenter814*xsize/668 + x00, LAEycenter814*ysize/668 + y0, color=6, thick=7
	endif
	
	if mask[i] eq 'FG2' then begin
                tv, bytscl(image606, min=-0.0185, max=0.035), startHST606, ypos - postageysize/2, xsize=postagexsize, ysize=postageysize, /normal
                tvbox, slitsize, (334+taroffsetx)*xsize/668 + x0, (334+taroffsety)*ysize/668 + y0, color=2, angle = -slitPA, thick=4
                tvcircle, 1./0.03*xsize/668, LAExcenter606*xsize/668 + x0, LAEycenter606*ysize/668 + y0, color=6, thick=7
		tv, bytscl(image814, min=-0.0314, max=0.0451), startHST814, ypos - postageysize/2, xsize=postagexsize, ysize=postageysize, /normal
                tvbox, slitsize, (334+taroffsetx)*xsize/668 + x00, (334+taroffsety)*ysize/668 + y0, color=2, angle = -slitPA, thick=4
                tvcircle, 1./0.03*xsize/668, LAExcenter814*xsize/668 + x00, LAEycenter814*ysize/668 + y0, color=6, thick=7
	endif
	endfor
	device, /close_file




	set_plot, 'PS'
        device, /color, xoffset=0.1, yoffset=0.2, xsize=10.7, ysize=14, /inches, filename='LYAspectralmosaic2.ps'
	
	for i=6, 12-1 do begin
	
	spec[i] = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/' + spec[i]
        corrspec = fill_gap(spec[i],/tweak,/telluric,/silent,header=header)
        wave = corrspec.lambda
        airtovac, wave
        corrspec.lambda = wave

        range = where(corrspec.lambda ge (1+z[i])*1215.7-20 and corrspec.lambda le (1+z[i])*1215.7+20)

        constfluxcorr = double(2e-08/(3600.*!pi*(449.)^2*0.33*3e18*1e-29)) ;3e+18 is the speed of light in Ang/s and 10e-29 ergs/s/cm^2/Hz = 1 uJy

        lambda = corrspec.lambda[range]

        uJyspec = dblarr(n_elements(corrspec.spec[range]))	

	        for j=0, n_elements(range)-1 do begin
                uJyspec[j] = constfluxcorr*corrspec.spec[range[j]]*lambda[[j]];spectrum in microJanskys, converted from ergs/s/cm^2/Angstom
                endfor

        ypos = (5-(i-6))/6.2 + 0.1
	xyouts, 0.02, ypos, mask[i] + '.' + slit[i], charsize=1.5, charthick=6, /normal
	
	if i eq 11 then begin

        plot, corrspec.lambda[range], uJyspec, yrange=[-1.0, max(uJyspec)], ystyle=1, position=[start1d, ypos-0.057, end1d, ypos+0.057],$
        color=0, xrange = [min(lambda), max(lambda)], xstyle=1, xtitle=textoidl('\lambda (\AA)') , ytickformat='(f5.1)', $
        thick=6, charthick=4, charsize=1.2, xthick=6, ythick=6, xticks=2, yticks=2, xminor=5, yminor=5, xticklen=0.11, yticklen=0.11, /noerase  
        endif

	;xyouts, 0.18, 0.5, textoidl('f_{\nu} \muJy'), charthick=4, charsize=1.2       
 
        if i ne 11 then begin
        plot, corrspec.lambda[range], uJyspec, yrange=[-1.0, max(uJyspec)], ystyle=1, position=[start1d, ypos-0.057, end1d, ypos+0.057],$
        color=0, xrange = [min(lambda), max(lambda)], xstyle=1, ytickformat='(f5.1)', $ 
        thick=6, charthick=4, charsize=1.2, xthick=6, ythick=6, xticks=2, yticks=2, xminor=5, yminor=5, xticklen=0.11, yticklen=0.11, /noerase 
        endif

	twodcutout = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/1604_LYA/mosaic/2dcutouts/cutout.' + mask[i] + '.0' + slit[i] + '.fits'
        twodspec = mrdfits(twodcutout,1, /silent)
        cutoutinv = twodspec.flux*(-1)
        tv, bytscl(cutoutinv, min=-13.51, max=20.5), start2d, ypos - 0.04, xsize=twodxsize, ysize=twodysize, /normal

	 if mask[i] + '.' + slit[i] ne 'SC2NM1.45' and mask[i] + '.' + slit[i] ne 'SC2NM2.61' then begin
                f606w = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/1604_LYA/mosaic/postage_stamps/ACScutout.f606w.' + mask[i] + '.0' + slit[i] + '.fits'
                f814w = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/1604_LYA/mosaic/postage_stamps/ACScutout.f814w.' + mask[i] + '.0' + slit[i] + '.fits'

                image606 = mrdfits(f606w,0, /silent)*(-1)
                image814 = mrdfits(f814w,0, /silent)*(-1)

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

                cdmatrix606 = dblarr(2,2)
                cdmatrix606[0,0] = f606waa
                cdmatrix606[0,1] = f606wab
                cdmatrix606[1,0] = f606wba
                cdmatrix606[1,1] = f606wbb
                cdinvmatrix606 = invert(cdmatrix606)

                cdmatix814 = dblarr(2,2)
                cdmatix814[0,0] = f814waa
                cdmatix814[0,1] = f814wab
                cdmatix814[1,0] = f814wba
                cdmatix814[1,1] = f814wbb
                cdinvmatrix = invert(cdmatix814)

                slitPA = double(sxpar(f606whdr,'slitPA'))
                slitlen = double(sxpar(f606whdr,'SLITLEN'))

                orientat606 = double(sxpar(f606whdr,'ORIENTAT'))
                orientat814 = double(sxpar(f814whdr,'ORIENTAT'))

                ; LAE position in 606 image     

		 p = atan( (-cos(HSTLAEdec[i]*2*!PI/360)*sin( (f606wCRVAL1-HSTLAEra[i])*2*!PI/360 ))/(cos(f606wCRVAL2*2*!PI/360)*sin(HSTLAEdec[i]*2*!PI/360)-sin(f606wCRVAL2*2*!PI/360)*cos(HSTLAEdec[i]*2*!PI/360)*cos( 2*!PI/360*(f606wCRVAL1 - HSTLAEra[i]))))
                if p lt 0 and f606wCRVAL1 gt HSTLAEra[i] then p = 2*!PI + p else if p lt 0 and f606wCRVAL1 lt HSTLAEra[i] then p = !PI + p else if p gt 0 and f606wCRVAL1 gt HSTLAEra[i] then p = !PI + p else if p gt 0 and f606wCRVAL1 lt HSTLAEra[i] then p = p

                deltadis = 360/(2*!PI)*acos(sin(HSTLAEdec[i]*2*!PI/360)*sin(f606wCRVAL2*2*!PI/360)+cos(f606wCRVAL2*2*!PI/360)*cos(HSTLAEdec[i]*2*!PI/360)*cos(2*!PI/360*(f606wCRVAL1 - HSTLAEra[i])))
                deltaRA = deltadis*sin(p)
                deltadec = deltadis*cos(p)

                coormatrix =dblarr(2)
                coormatrix[0] = deltaRA
                coormatrix[1] = deltadec

                ;transform between difference in Ra/Dec to x/y pixels
                newpixarray = cdinvmatrix # coormatrix
                LAExcenter606 = floor(f606wrefxpix + newpixarray[0])
                LAEycenter606 = floor(f606wrefypix + newpixarray[1])


                ; LAE position in 814 image     
                p = atan( (-cos(HSTLAEdec[i]*2*!PI/360)*sin( (f814wCRVAL1-HSTLAEra[i])*2*!PI/360 ))/(cos(f814wCRVAL2*2*!PI/360)*sin(HSTLAEdec[i]*2*!PI/360)-sin(f814wCRVAL2*2*!PI/360)*cos(HSTLAEdec[i]*2*!PI/360)*cos( 2*!PI/360*(f814wCRVAL1 - HSTLAEra[i]))))
                if p lt 0 and f814wCRVAL1 gt HSTLAEra[i] then p = 2*!PI + p else if p lt 0 and f814wCRVAL1 lt HSTLAEra[i] then p = !PI + p else if p gt 0 and f814wCRVAL1 gt HSTLAEra[i] then p = !PI + p else if p gt 0 and f814wCRVAL1 lt HSTLAEra[i] then p = p

                deltadis = 360/(2*!PI)*acos(sin(HSTLAEdec[i]*2*!PI/360)*sin(f814wCRVAL2*2*!PI/360)+cos(f814wCRVAL2*2*!PI/360)*cos(HSTLAEdec[i]*2*!PI/360)*cos(2*!PI/360*(f814wCRVAL1 - HSTLAEra[i])))
                deltaRA = deltadis*sin(p)
                deltadec = deltadis*cos(p)

                coormatrix =dblarr(2)
                coormatrix[0] = deltaRA
                coormatrix[1] = deltadec

                ;transform between difference in Ra/Dec to x/y pixels
                newpixarray = cdinvmatrix # coormatrix
                LAExcenter814 = floor(f814wrefxpix + newpixarray[0])
                LAEycenter814 = floor(f814wrefypix + newpixarray[1])
		
	endif

 	print, mask[i] + '.' + slit[i]

	windowXSize = !d.x_size
        windowYSize = !d.y_size
        x0 = (startHST606) * windowXSize
        x00 = (startHST814) * windowXSize
	y0 = (ypos +  postageysize/2 - 0.134 + 0.2/14) * windowYSize
	xsize = postagexsize * windowXSize
        ysize = postageysize * windowYSize       

	slitsize = [1/0.03*xsize/668, slitlen/0.03*xsize/668]

	if mask[i] ne '16XR2' then begin
        taroffsetx =  (yslit[i]/2 -tarypos[i])*0.118*sin(slitPA*2*!PI/360)/0.03
        taroffsety = (yslit[i]/2 -tarypos[i])*0.118*cos(slitPA*2*!PI/360)/0.03
	endif
	
	if mask[i] eq '16XR2' then begin
	taroffsetx =  (yslit[i]/2 -tarypos[i])*0.118*sin(-slitPA*2*!PI/360)/0.03
        taroffsety = (yslit[i]/2 -tarypos[i])*0.118*cos(-slitPA*2*!PI/360)/0.03
	endif

 
	if mask[i] + '.' + slit[i] ne 'SC2NM1.45' and mask[i] + '.' + slit[i] ne 'SC2NM2.61' then begin
                f606w = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/1604_LYA/mosaic/postage_stamps/ACScutout.f606w.' + mask[i] + '.0' + slit[i] + '.fits'
                f814w = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/1604_LYA/mosaic/postage_stamps/ACScutout.f814w.' + mask[i] + '.0' + slit[i] + '.fits'

                image606 = mrdfits(f606w,0, /silent)*(-1)
                image814 = mrdfits(f814w,0, /silent)*(-1)

                tv, bytscl(image606, min=-0.029, max=0.045), startHST606, ypos - postageysize/2, xsize=postagexsize, ysize=postageysize, /normal
                tvbox, slitsize, (334+taroffsetx)*xsize/668 + x0, (334+taroffsety)*ysize/668 + y0, color=2, angle = -slitPA, thick=4

                if mask[i] + '.' + slit[i] ne '16XR1.72' and mask[i] + '.' + slit[i] ne 'SC2NM1.78' then begin
		tvcircle, 1./0.03*xsize/668, LAExcenter606*xsize/668 + x0, LAEycenter606*ysize/668 + y0, color=6, thick=7
		endif

		if mask[i] + '.' + slit[i] eq '16XR1.72' then begin
                tvcircle, 1./0.03*xsize/668, 334*xsize/668 + x0, 334*ysize/668 + y0, color=6, thick=7                             
                endif

		if mask[i] + '.' + slit[i] eq 'SC2NM1.78' then begin
                tvcircle, 1./0.03*xsize/668, (334+10)*xsize/668 + x0, (334-75)*ysize/668 + y0, color=6, thick=7
                endif

		tv, bytscl(image814, min=-0.0418, max=0.0551), startHST814, ypos - postageysize/2, xsize=postagexsize, ysize=postageysize, /normal
		tvbox, slitsize, (334+taroffsetx)*xsize/668 + x00, (334+taroffsety)*ysize/668 + y0, color=2, angle = -slitPA, thick=4
                
		if mask[i] + '.' + slit[i] ne '16XR1.72' and mask[i] + '.' + slit[i] ne 'SC2NM1.78' then begin
		tvcircle, 1./0.03*xsize/668, LAExcenter814*xsize/668 + x00, LAEycenter814*ysize/668 + y0, color=6, thick=7	
		endif

		if mask[i] + '.' + slit[i] eq '16XR1.72' then begin
                tvcircle, 1./0.03*xsize/668, 334*xsize/668 + x00, 334*ysize/668 + y0, color=6, thick=7       
                endif
		
		if mask[i] + '.' + slit[i] eq 'SC2NM1.78' then begin
		tvcircle, 1./0.03*xsize/668, (334+10)*xsize/668 + x00, (334-75)*ysize/668 + y0, color=6, thick=7
		endif

        endif

	if mask[i] + '.' + slit[i] eq 'SC2NM1.45' or mask[i] + '.' + slit[i] eq 'SC2NM2.61' then begin

		rprime = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/1604_LYA/mosaic/postage_stamps/LFCcutout.r.' + mask[i] + '.0' + slit[i] + '.fits'
                iprime = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/1604_LYA/mosaic/postage_stamps/LFCcutout.i.' + mask[i] + '.0' + slit[i] + '.fits'
                zprime = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/1604_LYA/mosaic/postage_stamps/LFCcutout.z.' + mask[i] + '.0' + slit[i] + '.fits'

                rhdr = headfits(rprime)

                rCRVAL1 = double(sxpar(rhdr,'CRVAL1'))    ;initial RA for any movement along the image
                rCRVAL2 = double(sxpar(rhdr, 'CRVAL2'))   ;initial dec for any movement along the image

                raa = double(sxpar(rhdr, 'CD1_1'))        ;get the transformation matrix between pixels and global coordinates 
                rab = double(sxpar(rhdr, 'CD1_2'))
                rba = double(sxpar(rhdr, 'CD2_1'))
                rbb = double(sxpar(rhdr, 'CD2_2'))

                rrefxpix = sxpar(rhdr, 'CRPIX1')          ;initial x pix corresponding to CRVAL1
                rrefypix = sxpar(rhdr, 'CRPIX2')
                
		cdmatrixr = dblarr(2,2)
                cdmatrixr[0,0] = raa
                cdmatrixr[0,1] = rab
                cdmatrixr[1,0] = rba
                cdmatrixr[1,1] = rbb
                cdinvmatrixr = invert(cdmatrixr)

                slitPA = double(sxpar(rhdr,'slitPA'))
                slitlen = double(sxpar(rhdr,'SLITLEN'))

		windowXSize = !d.x_size
 	       	windowYSize = !d.y_size
 	       	x0 = startr * windowXSize
       		x00 = starti * windowXSize
        	x000 = startz * windowXSize
		y0 = (ypos + postageysize/2 - 0.134 + 0.2/14 ) * windowYSize
        	xsize = postagexsize * windowXSize
        	ysize = postageysize * windowYSize

        	slitsize = [1/0.18*xsize/112, slitlen/0.18*xsize/112]
	
			
		if mask[i] + '.' + slit[i] eq 'SC2NM1.45' then begin 			

			 p = atan( (-cos(LFCLAEdec[0]*2*!PI/360)*sin( (rCRVAL1-LFCLAEra[0])*2*!PI/360 ))/(cos(rCRVAL2*2*!PI/360)*sin(LFCLAEdec[0]*2*!PI/360)-sin(rCRVAL2*2*!PI/360)*cos(LFCLAEdec[0]*2*!PI/360)*cos( 2*!PI/360*(rCRVAL1 - LFCLAEra[0]))))
	                if p lt 0 and rCRVAL1 gt LFCLAEra[0] then p = 2*!PI + p else if p lt 0 and rCRVAL1 lt LFCLAEra[0] then p = !PI + p else if p gt 0 and rCRVAL1 gt LFCLAEra[0] then p = !PI + p else if p gt 0 and rCRVAL1 lt LFCLAEra[0] then p = p

        	        deltadis = 360/(2*!PI)*acos(sin(LFCLAEdec[0]*2*!PI/360)*sin(rCRVAL2*2*!PI/360)+cos(rCRVAL2*2*!PI/360)*cos(LFCLAEdec[0]*2*!PI/360)*cos(2*!PI/360*(rCRVAL1 - LFCLAEra[0])))
        	        deltaRA = deltadis*sin(p)
        	        deltadec = deltadis*cos(p)

        	        coormatrix =dblarr(2)
        	        coormatrix[0] = deltaRA
        	        coormatrix[1] = deltadec

        	        ;transform between difference in Ra/Dec to x/y pixels
        	        newpixarray = cdinvmatrixr # coormatrix
			LAExcenterr = floor(rrefxpix + newpixarray[0])
        	        LAEycenterr = floor(rrefypix + newpixarray[1])                	


			taroffsetx =  (LFCyslit[0]/2 -LFCtarypos[0])*0.118*sin(-slitPA*2*!PI/360)/0.18
 			taroffsety = (LFCyslit[0]/2 -LFCtarypos[0])*0.118*cos(-slitPA*2*!PI/360)/0.18

			imager = mrdfits(rprime,/silent)*(-1)
                	imagei = mrdfits(iprime, /silent)*(-1)
                	imagez = mrdfits(zprime, /silent)*(-1)

                	tv, bytscl(imager, min=-4.411, max=-4.117), startr, ypos - postageysize/2, xsize=postagexsize, ysize=postageysize, /normal
                	tvbox, slitsize, (56+taroffsetx)*xsize/112 + x0, (56+taroffsety)*ysize/112 + y0, color=2, angle = -slitPA, thick=4
                	tvcircle, 1./0.18*xsize/112, LAExcenterr*xsize/112 + x0, LAEycenterr*ysize/112 + y0, color=6, thick=7
                	tv, bytscl(imagei, min=-13.19, max=-12.67), starti, ypos - postageysize/2, xsize=postagexsize, ysize=postageysize, /normal
                	tvbox, slitsize, (56+taroffsetx)*xsize/112 + x00, (56+taroffsety)*ysize/112 + y0, color=2, angle = -slitPA, thick=4
                        tvcircle, 1./0.18*xsize/112, LAExcenterr*xsize/112 + x00, LAEycenterr*ysize/112 + y0, color=6, thick=7
			tv, bytscl(imagez, min=-10.473, max=-9.89), startz, ypos - postageysize/2, xsize=postagexsize, ysize=postageysize, /normal
			tvbox, slitsize, (56+taroffsetx)*xsize/112 + x000, (56+taroffsety)*ysize/112 + y0, color=2, angle = -slitPA, thick=4
                        tvcircle, 1./0.18*xsize/112, LAExcenterr*xsize/112 + x000, LAEycenterr*ysize/112 + y0, color=6, thick=7
		endif

        	if mask[i] + '.' + slit[i] eq 'SC2NM2.61' then begin

                	p = atan( (-cos(LFCLAEdec[1]*2*!PI/360)*sin( (rCRVAL1-LFCLAEra[1])*2*!PI/360 ))/(cos(rCRVAL2*2*!PI/360)*sin(LFCLAEdec[1]*2*!PI/360)-sin(rCRVAL2*2*!PI/360)*cos(LFCLAEdec[1]*2*!PI/360)*cos( 2*!PI/360*(rCRVAL1 - LFCLAEra[1]))))
                        if p lt 0 and rCRVAL1 gt LFCLAEra[1] then p = 2*!PI + p else if p lt 0 and rCRVAL1 lt LFCLAEra[1] then p = !PI + p else if p gt 0 and rCRVAL1 gt LFCLAEra[1] then p = !PI + p else if p gt 0 and rCRVAL1 lt LFCLAEra[1] then p = p

                        deltadis = 360/(2*!PI)*acos(sin(LFCLAEdec[1]*2*!PI/360)*sin(rCRVAL2*2*!PI/360)+cos(rCRVAL2*2*!PI/360)*cos(LFCLAEdec[1]*2*!PI/360)*cos(2*!PI/360*(rCRVAL1 - LFCLAEra[1])))
                        deltaRA = deltadis*sin(p)
                        deltadec = deltadis*cos(p)

                        coormatrix =dblarr(2)
                        coormatrix[0] = deltaRA
                        coormatrix[1] = deltadec

                        ;transform between difference in Ra/Dec to x/y pixels
                        newpixarray = cdinvmatrixr # coormatrix
                        LAExcenterr = floor(rrefxpix + newpixarray[0])-1	;have to trick it circle is not showing up right
                        LAEycenterr = floor(rrefypix + newpixarray[1])+12	

                        taroffsetx =  (LFCyslit[1]/2 -LFCtarypos[1])*0.118*sin(slitPA*2*!PI/360)/0.18
                        taroffsety = (LFCyslit[1]/2 -LFCtarypos[1])*0.118*cos(slitPA*2*!PI/360)/0.18

			imager = mrdfits(rprime, /silent)*(-1)
                	imagei = mrdfits(iprime, /silent)*(-1)
                	imagez = mrdfits(zprime, /silent)*(-1)

                	tv, bytscl(imager, min=-10.72, max=-9.48), startr, ypos - 0.05, xsize=postagexsize, ysize=postageysize, /normal
                	tvbox, slitsize, (56+taroffsetx)*xsize/112 + x0, (56+taroffsety)*ysize/112 + y0, color=2, angle = -slitPA, thick=4
                        tvcircle, 1./0.18*xsize/112, LAExcenterr*xsize/111 + x0, LAEycenterr*ysize/111 + y0, color=6, thick=7
			tv, bytscl(imagei, min=-9.25, max=-7.656), starti, ypos - 0.05, xsize=postagexsize, ysize=postageysize, /normal
			tvbox, slitsize, (56+taroffsetx)*xsize/112 + x00, (56+taroffsety)*ysize/111 + y0, color=2, angle = -slitPA, thick=4
                        tvcircle, 1./0.18*xsize/112, LAExcenterr*xsize/111 + x00, LAEycenterr*ysize/111 + y0, color=6, thick=7
                	tv, bytscl(imagez, min=-10.754, max=-9.43), startz, ypos - 0.05, xsize=postagexsize, ysize=postageysize, /normal
			tvbox, slitsize, (56+taroffsetx)*xsize/112 + x000, (56+taroffsety)*ysize/111 + y0, color=2, angle = -slitPA, thick=4
                        tvcircle, 1./0.18*xsize/112, LAExcenterr*xsize/111 + x000, LAEycenterr*ysize/111 + y0, color=6, thick=7
		endif	
		endif
        endfor

	device, /close_file




	set_plot, 'PS'
        device, /color, xoffset=0.1, yoffset=0.2, xsize=10.7, ysize=14, /inches, filename='LYAspectralmosaic3.ps'

        for i=12, n_elements(mask)-1 do begin

        spec[i] = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/' + spec[i]
        corrspec = fill_gap(spec[i],/tweak,/telluric,/silent,header=header)
        wave = corrspec.lambda
        airtovac, wave
        corrspec.lambda = wave

        range = where(corrspec.lambda ge (1+z[i])*1215.7-20 and corrspec.lambda le (1+z[i])*1215.7+20)

        constfluxcorr = double(2e-08/(3600.*!pi*(449.)^2*0.33*3e18*1e-29)) ;3e+18 is the speed of light in Ang/s and 10e-29 ergs/s/cm^2/Hz = 1 uJy

        lambda = corrspec.lambda[range]

        uJyspec = dblarr(n_elements(corrspec.spec[range]))

                for j=0, n_elements(range)-1 do begin
                uJyspec[j] = constfluxcorr*corrspec.spec[range[j]]*lambda[[j]];spectrum in microJanskys, converted from ergs/s/cm^2/Angstom
                endfor

        ypos = (4-(i-12))/5.2 + 0.1
        xyouts, 0.02, ypos, mask[i] + '.' + slit[i], charsize=1.5, charthick=6, /normal

        if i eq 16 then begin

        plot, corrspec.lambda[range], uJyspec, yrange=[-0.5, max(uJyspec)], ystyle=1, position=[start1d, ypos-0.057, end1d, ypos+0.057],$
        color=0, xrange = [min(lambda), max(lambda)], xstyle=1, xtitle=textoidl('\lambda (\AA)') , ytickformat = '(f5.1)', $
        thick=6, charthick=4, charsize=1.2, xthick=6, ythick=6, xticks=2, yticks=2, xminor=5, yminor=5, xticklen=0.11, yticklen=0.11, /noerase
        endif

        if i ne 16 then begin
        plot, corrspec.lambda[range], uJyspec, yrange=[-0.5, max(uJyspec)], ystyle=1, position=[start1d, ypos-0.057, end1d, ypos+0.057],$
        color=0, xrange = [min(lambda), max(lambda)], xstyle=1, $
        thick=6, charthick=4, charsize=1.2, xthick=6, ythick=6, xticks=2, yticks=2, xminor=5, yminor=5, xticklen=0.11, yticklen=0.11, /noerase
        endif

        twodcutout = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/1604_LYA/mosaic/2dcutouts/cutout.' + mask[i] + '.0' + slit[i] + '.fits'
        twodspec = mrdfits(twodcutout,1, /silent)
        cutoutinv = twodspec.flux*(-1)
        tv, bytscl(cutoutinv, min=-13.51, max=20.5), start2d, ypos - 0.04, xsize=twodxsize, ysize=twodysize, /normal

        f606w = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/1604_LYA/mosaic/postage_stamps/ACScutout.f606w.' + mask[i] + '.0' + slit[i] + '.fits'
        f814w = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/1604_LYA/mosaic/postage_stamps/ACScutout.f814w.' + mask[i] + '.0' + slit[i] + '.fits'
	
	image606 = mrdfits(f606w,0, /silent)*(-1)
        image814 = mrdfits(f814w,0, /silent)*(-1)

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
        
	cdmatrix606 = dblarr(2,2)
        cdmatrix606[0,0] = f606waa
        cdmatrix606[0,1] = f606wab
        cdmatrix606[1,0] = f606wba
        cdmatrix606[1,1] = f606wbb
        cdinvmatrix606 = invert(cdmatrix606)
        
	cdmatix814 = dblarr(2,2)
        cdmatix814[0,0] = f814waa
        cdmatix814[0,1] = f814wab
        cdmatix814[1,0] = f814wba
        cdmatix814[1,1] = f814wbb
        cdinvmatrix = invert(cdmatix814)
        
	slitPA = double(sxpar(f606whdr,'slitPA'))
        slitlen = double(sxpar(f606whdr,'SLITLEN'))
        
	orientat606 = double(sxpar(f606whdr,'ORIENTAT'))
        orientat814 = double(sxpar(f814whdr,'ORIENTAT'))

        ; LAE position in 606 image     

        p = atan( (-cos(HSTLAEdec[i]*2*!PI/360)*sin( (f606wCRVAL1-HSTLAEra[i])*2*!PI/360 ))/(cos(f606wCRVAL2*2*!PI/360)*sin(HSTLAEdec[i]*2*!PI/360)-sin(f606wCRVAL2*2*!PI/360)*cos(HSTLAEdec[i]*2*!PI/360)*cos( 2*!PI/360*(f606wCRVAL1 - HSTLAEra[i]))))
        if p lt 0 and f606wCRVAL1 gt HSTLAEra[i] then p = 2*!PI + p else if p lt 0 and f606wCRVAL1 lt HSTLAEra[i] then p = !PI + p else if p gt 0 and f606wCRVAL1 gt HSTLAEra[i] then p = !PI + p else if p gt 0 and f606wCRVAL1 lt HSTLAEra[i] then p = p
	

	deltadis = 360/(2*!PI)*acos(sin(HSTLAEdec[i]*2*!PI/360)*sin(f606wCRVAL2*2*!PI/360)+cos(f606wCRVAL2*2*!PI/360)*cos(HSTLAEdec[i]*2*!PI/360)*cos(2*!PI/360*(f606wCRVAL1 - HSTLAEra[i])))
        deltaRA = deltadis*sin(p)
        deltadec = deltadis*cos(p)

        coormatrix =dblarr(2)
        coormatrix[0] = deltaRA
        coormatrix[1] = deltadec

        ;transform between difference in Ra/Dec to x/y pixels
        newpixarray = cdinvmatrix # coormatrix
        LAExcenter606 = floor(f606wrefxpix + newpixarray[0])
        LAEycenter606 = floor(f606wrefypix + newpixarray[1])


        ; LAE position in 814 image     
        p = atan( (-cos(HSTLAEdec[i]*2*!PI/360)*sin( (f814wCRVAL1-HSTLAEra[i])*2*!PI/360 ))/(cos(f814wCRVAL2*2*!PI/360)*sin(HSTLAEdec[i]*2*!PI/360)-sin(f814wCRVAL2*2*!PI/360)*cos(HSTLAEdec[i]*2*!PI/360)*cos( 2*!PI/360*(f814wCRVAL1 - HSTLAEra[i]))))
        if p lt 0 and f814wCRVAL1 gt HSTLAEra[i] then p = 2*!PI + p else if p lt 0 and f814wCRVAL1 lt HSTLAEra[i] then p = !PI + p else if p gt 0 and f814wCRVAL1 gt HSTLAEra[i] then p = !PI + p else if p gt 0 and f814wCRVAL1 lt HSTLAEra[i] then p = p

        deltadis = 360/(2*!PI)*acos(sin(HSTLAEdec[i]*2*!PI/360)*sin(f814wCRVAL2*2*!PI/360)+cos(f814wCRVAL2*2*!PI/360)*cos(HSTLAEdec[i]*2*!PI/360)*cos(2*!PI/360*(f814wCRVAL1 - HSTLAEra[i])))
        deltaRA = deltadis*sin(p)
        deltadec = deltadis*cos(p)

        coormatrix =dblarr(2)
        coormatrix[0] = deltaRA
        coormatrix[1] = deltadec

        ;transform between difference in Ra/Dec to x/y pixels
        newpixarray = cdinvmatrix # coormatrix
        LAExcenter814 = floor(f814wrefxpix + newpixarray[0])
        LAEycenter814 = floor(f814wrefypix + newpixarray[1])
              
        print, mask[i] + '.' + slit[i]

        windowXSize = !d.x_size
        windowYSize = !d.y_size
        x0 = (startHST606) * windowXSize
        x00 = (startHST814) * windowXSize
        y0 = (ypos - 0.074  + 0.2/14) * windowYSize
        xsize = postagexsize * windowXSize
        ysize = postageysize * windowYSize


        slitsize = [1/0.03*xsize/668, slitlen/0.03*xsize/668]

	if mask[i] ne '16XR2' and mask[i] + '.' + slit[i] ne '16XR1.97' then begin
        taroffsetx =  (yslit[i]/2 -tarypos[i])*0.118*sin(slitPA*2*!PI/360)/0.03
        taroffsety = (yslit[i]/2 -tarypos[i])*0.118*cos(slitPA*2*!PI/360)/0.03
	endif
		
	if mask[i] eq '16XR2' then begin
	taroffsetx = (yslit[i]/2 -tarypos[i])*0.118*sin(-slitPA*2*!PI/360)/0.03
        taroffsety = (yslit[i]/2 -tarypos[i])*0.118*(-cos(slitPA*2*!PI/360)/0.03) + 40	;not coming out right have to trick it
	endif

	if mask[i] + '.' + slit[i] eq '16XR1.97' then begin
	taroffsetx = (yslit[i]/2 -tarypos[i])*0.118*sin(slitPA*2*!PI/360)/0.03
        taroffsety = (yslit[i]/2 -tarypos[i])*0.118*(-cos(slitPA*2*!PI/360)/0.03) + 32  ;not coming out right have to trick it
        endif


        f606w = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/1604_LYA/mosaic/postage_stamps/ACScutout.f606w.' + mask[i] + '.0' + slit[i] + '.fits'
        f814w = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/1604_LYA/mosaic/postage_stamps/ACScutout.f814w.' + mask[i] + '.0' + slit[i] + '.fits'

        image606 = mrdfits(f606w,0, /silent)*(-1)
        image814 = mrdfits(f814w,0, /silent)*(-1)

        tv, bytscl(image606, min=-0.029, max=0.045), startHST606, ypos - 0.05, xsize=postagexsize, ysize=postageysize, /normal
        tvbox, slitsize, (334+taroffsetx)*xsize/668 + x0, (334+taroffsety)*ysize/668 + y0, color=2, angle = -slitPA, thick=4

        if mask[i] + '.' + slit[i] ne '16XR1.72' and mask[i] ne '16XR2' then begin
                tvcircle, 1./0.03*xsize/668, (LAExcenter606+20)*xsize/668 + x0, LAEycenter606*ysize/668 + y0, color=6, thick=7	;+20 on 606 and 814 x cuz the circles are off again
        endif

        if mask[i] + '.' + slit[i] eq '16XR1.72' then begin
                tvcircle, 1./0.03*xsize/668, (334+taroffsetx)*xsize/668 + x0, (355+taroffsety)*ysize/668 + y0, color=6, thick=7
        endif

	if mask[i] eq '16XR2' then begin
		tvcircle, 1./0.03*xsize/668, LAExcenter606*xsize/668 + x0, (LAEycenter606+40)*ysize/668 + y0, color=6, thick=7
	endif

        tv, bytscl(image814, min=-0.0418, max=0.0551), startHST814, ypos - 0.05, xsize=postagexsize, ysize=postageysize, /normal
        tvbox, slitsize, (334+taroffsetx)*xsize/668 + x00, (334+taroffsety)*ysize/668 + y0, color=2, angle = -slitPA, thick=4
        
	if mask[i] + '.' + slit[i] ne '16XR1.72' and mask[i] ne '16XR2' then begin
                tvcircle, 1./0.03*xsize/668, (LAExcenter814+20)*xsize/668 + x00, LAEycenter814*ysize/668 + y0, color=6, thick=7
        endif

        if mask[i] + '.' + slit[i] eq '16XR1.72' then begin
        tvcircle, 1./0.03*xsize/668, (334+taroffsetx)*xsize/668 + x00, (355+taroffsety)*ysize/668 + y0, color=6, thick=7
	endif


	if mask[i] eq '16XR2' then begin
                tvcircle, 1./0.03*xsize/668, LAExcenter814*xsize/668 + x00, (LAEycenter814+40)*ysize/668 + y0, color=6, thick=7
        endif	
	
        endfor

	device, /close_file

end	
