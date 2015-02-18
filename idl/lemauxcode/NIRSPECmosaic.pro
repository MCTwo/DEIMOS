pro NIRSPECmosaic
	
	set_plot, 'PS'
        loadcolors
	device, /color, xoffset=0.1, yoffset=0.2, xsize=7.5, ysize=10, /inches, filename='NIRSPECmosaic1.ps'
	
	readcol, 'linefluxnEW_mag_all.cat', format='A,A,D,D,A', mask, slit, z, EWHa, EWHaerr, EWOII, EWOIIerr, FHa, FHaerr, FNII, FNIIerr, rmag, imag, zmag, f606mag, f814mag, RA, dec, /silent
	
	; doing the first 6, 2nd 6, 3rd 6, and then the last 2 on a seperate plot

	postagexsize = 0.225
	postageysize = 0.17625
	startHST606 = 0.18
	startHST814 = 0.425

	lineratio = fltarr(n_elements(mask))

	for i=0, n_elements(mask)-1 do begin
		if FHaerr[i] gt 0 then lineratio[i] = alog10(FNII[i]/FHa[i]) else lineratio[i] = 99.	;If using only 3sigma line fluxes then set the NII/Ha line ratio to ambiguous
	endfor

	EWdelim = -5.*EWHa/(1.+z) - 8	;if -1*EWOII is greater than this then high OII/Ha else low OII/Ha
	posOII = -1*EWOII

	for i=0, 6-1 do begin

	ypos = (5-i)/5.5 + 0.07
        xyouts, 0.01, ypos, mask[i] + '.' + slit[i], charsize=1.3, charthick=5, /normal
	xyouts, 0.008, ypos-0.03, [i+1], charsize=1.5, charthick=6, color=1, /normal
	if lineratio[i] ge -0.22 and lineratio[i] lt 98. then begin
	xyouts, 0.655, ypos, 'AGN NII/Ha', charsize=1.5, charthick=6, color=5, /normal
	endif
	
	if lineratio[i] lt -0.22 then begin
	xyouts, 0.655, ypos, 'SF NII/Ha', charsize=1.5, charthick=6, color=6, /normal
        endif

	if lineratio[i] gt 98. then begin
	xyouts, 0.655, ypos, 'Ambiguous NII/Ha', charsize=1.5, charthick=6, color=4, /normal
        endif

	if posOII[i] ge EWdelim[i] then begin
	xyouts, 0.91, ypos, 'High OII/Ha', charsize=1.5, charthick=6, color=5, /normal
	endif

	if posOII[i] lt EWdelim[i] then begin
	xyouts, 0.91, ypos, 'Low OII/Ha', charsize=1.5, charthick=6, color=6, /normal
        endif

	f606w = '/Volumes/Data2/orelse/lemaux/NIRSPEC/06042007/2007jun04/reduced/final_spec/final_data/mosaic/postage_stamps/ACScutout.f606w.' + mask[i] + '.' + slit[i] + '.fits'
        f814w = '/Volumes/Data2/orelse/lemaux/NIRSPEC/06042007/2007jun04/reduced/final_spec/final_data/mosaic/postage_stamps/ACScutout.f814w.' + mask[i] + '.' + slit[i] + '.fits'

	image606 = mrdfits(f606w,0, /silent)*(-1)
	image814 = mrdfits(f814w,0, /silent)*(-1)

	f606whdr = headfits(f606w)
	f814whdr = headfits(f814w)
	
	f606wCRVAL1 = double(sxpar(f606whdr,'CRVAL1'))    ;initial RA for any movement along the image
        f606wCRVAL2 = double(sxpar(f606whdr, 'CRVAL2'))   ;initial dec for any movement along the image

        f814wCRVAL1 = double(sxpar(f814whdr,'CRVAL1'))
        f814wCRVAL2 = double(sxpar(f814whdr, 'CRVAL2'))


	windowXSize = !d.x_size
	windowYSize = !d.y_size
	x0 = (startHST606) * windowXSize
	x00 = (startHST814) * windowXSize	
	y0 = (ypos - postageysize/2 - 0.014 + 0.2/14) * windowYSize
	xsize = postagexsize * windowXSize
	ysize = postageysize * windowYSize		

	
		tv, bytscl(image606, min=-0.029, max=0.045), startHST606, ypos - postageysize/2, xsize=postagexsize, ysize=postageysize, /normal
		;tvbox, slitsize, (334+taroffsetx)*xsize/668 + x0, (334+taroffsety)*ysize/668 + y0, color=2, angle = -slitPA, thick=4
		;tvcircle, 1./0.03*xsize/668, LAExcenter606*xsize/668 + x0, LAEycenter606*ysize/668 + y0, color=6, thick=7
		tv, bytscl(image814, min=-0.0418, max=0.0551), startHST814, ypos - postageysize/2, xsize=postagexsize, ysize=postageysize, /normal
		;tvbox, slitsize, (334+taroffsetx)*xsize/668 + x00, (334+taroffsety)*ysize/668 + y0, color=2, angle = -slitPA, thick=4
        	;tvcircle, 1./0.03*xsize/668, LAExcenter814*xsize/668 + x00, LAEycenter814*ysize/668 + y0, color=6, thick=7
	;endif
	
	endfor
	device, /close_file




	set_plot, 'PS'
        device, /color, xoffset=0.1, yoffset=0.2, xsize=7.5, ysize=10, /inches, filename='NIRSPECmosaic2.ps'
	
	for i=6, 12-1 do begin
	
        ypos = (5-(i-6))/5.5 + 0.075

	xyouts, 0.01, ypos, mask[i] + '.' + slit[i], charsize=1.3, charthick=5, /normal
	xyouts, 0.008, ypos-0.03, [i+1], charsize=1.5, charthick=6, color=1, /normal

        if lineratio[i] ge -0.22 and lineratio[i] lt 98. then begin
        xyouts, 0.655, ypos, 'AGN NII/Ha', charsize=1.5, charthick=6, color=5, /normal
        endif

        if lineratio[i] lt -0.22 then begin
        xyouts, 0.655, ypos, 'SF NII/Ha', charsize=1.5, charthick=6, color=6, /normal
        endif

        if lineratio[i] gt 98. then begin
        xyouts, 0.655, ypos, 'Ambiguous NII/Ha', charsize=1.5, charthick=6, color=4, /normal
        endif

        if posOII[i] ge EWdelim[i] then begin
        xyouts, 0.93, ypos, 'High OII/Ha', charsize=1.5, charthick=6, color=5, /normal
        endif

        if posOII[i] lt EWdelim[i] then begin
        xyouts, 0.93, ypos, 'Low OII/Ha', charsize=1.5, charthick=6, color=6, /normal
        endif

	
	f606w = '/Volumes/Data2/orelse/lemaux/NIRSPEC/06042007/2007jun04/reduced/final_spec/final_data/mosaic/postage_stamps/ACScutout.f606w.' + mask[i] + '.' + slit[i] + '.fits'
        f814w = '/Volumes/Data2/orelse/lemaux/NIRSPEC/06042007/2007jun04/reduced/final_spec/final_data/mosaic/postage_stamps/ACScutout.f814w.' + mask[i] + '.' + slit[i] + '.fits'

                image606 = mrdfits(f606w,0, /silent)*(-1)
                image814 = mrdfits(f814w,0, /silent)*(-1)

                f606whdr = headfits(f606w)
                f814whdr = headfits(f814w)

		

 	print, mask[i] + '.' + slit[i]

	windowXSize = !d.x_size
        windowYSize = !d.y_size
        x0 = (startHST606) * windowXSize
        x00 = (startHST814) * windowXSize
	y0 = (ypos +  postageysize/2 - 0.134 + 0.2/14) * windowYSize
	xsize = postagexsize * windowXSize
        ysize = postageysize * windowYSize       


 
	  tv, bytscl(image606, min=-0.029, max=0.045), startHST606, ypos - postageysize/2, xsize=postagexsize, ysize=postageysize, /normal
                ;tvbox, slitsize, (334+taroffsetx)*xsize/668 + x0, (334+taroffsety)*ysize/668 + y0, color=2, angle = -slitPA, thick=4
                ;tvcircle, 1./0.03*xsize/668, LAExcenter606*xsize/668 + x0, LAEycenter606*ysize/668 + y0, color=6, thick=7
                tv, bytscl(image814, min=-0.0418, max=0.0551), startHST814, ypos - postageysize/2, xsize=postagexsize, ysize=postageysize, /normal

	endfor

	device, /close_file




	set_plot, 'PS'
        device, /color, xoffset=0.1, yoffset=0.2, xsize=7.5, ysize=10, /inches, filename='NIRSPECmosaic3.ps'

        for i=12, 18-1 do begin

        ypos = (5-(i-12))/5.5 + 0.04
        
	xyouts, 0.01, ypos+0.03, mask[i] + '.' + slit[i], charsize=1.3, charthick=5, /normal
	xyouts, 0.008, ypos, [i+1], charsize=1.5, charthick=6, color=1, /normal

        if lineratio[i] ge -0.22 and lineratio[i] lt 98. then begin
        xyouts, 0.655, ypos+0.03, 'AGN NII/Ha', charsize=1.5, charthick=6, color=5, /normal
        endif

        if lineratio[i] lt -0.22 then begin
        xyouts, 0.655, ypos+0.03, 'SF NII/Ha', charsize=1.5, charthick=6, color=6, /normal
        endif

        if lineratio[i] gt 98. then begin
        xyouts, 0.655, ypos+0.03, 'Ambiguous NII/Ha', charsize=1.5, charthick=6, color=4, /normal
        endif

        if posOII[i] ge EWdelim[i] then begin
        xyouts, 0.93, ypos+0.03, 'High OII/Ha', charsize=1.5, charthick=6, color=5, /normal
        endif

        if posOII[i] lt EWdelim[i] then begin
        xyouts, 0.93, ypos+0.03, 'Low OII/Ha', charsize=1.5, charthick=6, color=6, /normal
        endif
	

	f606w = '/Volumes/Data2/orelse/lemaux/NIRSPEC/06042007/2007jun04/reduced/final_spec/final_data/mosaic/postage_stamps/ACScutout.f606w.' + mask[i] + '.' + slit[i] + '.fits'
        f814w = '/Volumes/Data2/orelse/lemaux/NIRSPEC/06042007/2007jun04/reduced/final_spec/final_data/mosaic/postage_stamps/ACScutout.f814w.' + mask[i] + '.' + slit[i] + '.fits'

	
	image606 = mrdfits(f606w,0, /silent)*(-1)
        image814 = mrdfits(f814w,0, /silent)*(-1)

        f606whdr = headfits(f606w)
        f814whdr = headfits(f814w)

        print, mask[i] + '.' + slit[i]

        windowXSize = !d.x_size
        windowYSize = !d.y_size
        x0 = (startHST606) * windowXSize
        x00 = (startHST814) * windowXSize
        y0 = (ypos - 0.074  + 0.2/14) * windowYSize
        xsize = postagexsize * windowXSize
        ysize = postageysize * windowYSize



        tv, bytscl(image606, min=-0.029, max=0.045), startHST606, ypos - 0.05, xsize=postagexsize, ysize=postageysize, /normal

        tv, bytscl(image814, min=-0.0418, max=0.0551), startHST814, ypos - 0.05, xsize=postagexsize, ysize=postageysize, /normal
        

        endfor

	device, /close_file


	set_plot, 'PS'
        device, /color, xoffset=0.1, yoffset=0.2, xsize=7.5, ysize=10, /inches, filename='NIRSPECmosaic4.ps'

        for i=18, n_elements(mask)-1 do begin
 
        ypos = (3-(i-18))/5.5 + 0.07

	xyouts, 0.01, ypos+0.03, mask[i] + '.' + slit[i], charsize=1.3, charthick=5, /normal
	xyouts, 0.008, ypos, [i+1], charsize=1.5, charthick=6, color=1, /normal
	
        if lineratio[i] ge -0.22 and lineratio[i] lt 98. then begin
        xyouts, 0.655, ypos+0.03, 'AGN NII/Ha', charsize=1.5, charthick=6, color=5, /normal
        endif

        if lineratio[i] lt -0.22 then begin
        xyouts, 0.655, ypos+0.03, 'SF NII/Ha', charsize=1.5, charthick=6, color=6, /normal
        endif

        if lineratio[i] gt 98. then begin
        xyouts, 0.655, ypos+0.03, 'Ambiguous NII/Ha', charsize=1.5, charthick=6, color=4, /normal
        endif

        if posOII[i] ge EWdelim[i] then begin
        xyouts, 0.91, ypos+0.03, 'High OII/Ha', charsize=1.5, charthick=6, color=5, /normal
        endif

        if posOII[i] lt EWdelim[i] then begin
        xyouts, 0.91, ypos+0.03, 'Low OII/Ha', charsize=1.5, charthick=6, color=6, /normal
        endif



        f606w = '/Volumes/Data2/orelse/lemaux/NIRSPEC/06042007/2007jun04/reduced/final_spec/final_data/mosaic/postage_stamps/ACScutout.f606w.' + mask[i] + '.' + slit[i] + '.fits'
        f814w = '/Volumes/Data2/orelse/lemaux/NIRSPEC/06042007/2007jun04/reduced/final_spec/final_data/mosaic/postage_stamps/ACScutout.f814w.' + mask[i] + '.' + slit[i] + '.fits'

        image606 = mrdfits(f606w,0, /silent)*(-1)
        image814 = mrdfits(f814w,0, /silent)*(-1)

        f606whdr = headfits(f606w)
        f814whdr = headfits(f814w)


        print, mask[i] + '.' + slit[i]

        windowXSize = !d.x_size
        windowYSize = !d.y_size
        x0 = (startHST606) * windowXSize
        x00 = (startHST814) * windowXSize
        y0 = (ypos - 0.074  + 0.2/14) * windowYSize
        xsize = postagexsize * windowXSize
        ysize = postageysize * windowYSize



        tv, bytscl(image606, min=-0.029, max=0.045), startHST606, ypos - 0.05, xsize=postagexsize, ysize=postageysize, /normal

        tv, bytscl(image814, min=-0.0418, max=0.0551), startHST814, ypos - 0.05, xsize=postagexsize, ysize=postageysize, /normal


	endfor

        device, /close_file
end	
