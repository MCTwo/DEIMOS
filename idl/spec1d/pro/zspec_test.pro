;+
;  NAME:
;    zspec
; 
;  PURPOSE:
;    view reduced deep2 spectra and manually choose a best redshift.
;
;  CALLING SEQUENCE:
;    zspec, maskname [, /local_only, /nodeep, _extra=extra]
;
;  INPUTS:
;    maskname   - the masknumber to be viewed e.g. 2200
;
;  OPTIONAL INPUTS:
;    _extra = optional keywords to be passed to ATV and
;             fill_gap.pro. examples of possible parameters are listed
;             below.
;               min=min (set the minimum display value for ATV)
;               max=max (set the maximum display value for ATV)
;               boxsprof=1 (use the boxcar 1-d extraction - fill_gap).
;
;  KEYWORDS:
;    local_only - specify this keyword to force the program to ignore
;                 the non-local sky reductions. this must be set for
;                 all masks which lack non-local sky data.
;    nodeep     - specify this keyword for all non-DEEP2 masks.
;
;  OUTPUTS:
;    user can select to save the results to a fits file.  This file
;    also includes an ascii header with the same results for easy
;    access.
;
;  NOTES:
;    (1) the nodeep and local_only keywords are NOT required for the
;        KTRS masks. zspec will automatically detect these masks and
;        set nodeep/local_only accordingly.
;
;  REVISION HISTORY:
;	Created by DSM
;	July, 2003 - Various revisions requested by Faber and
;                    implemented by npk@astro.ucsc.edu
;       August, 2003 - renamed autosave files, redid logic for
;                    choosing which postage stamps to show in
;                    zspec_1d_plot and zspec_2d_plot - bjw
;	January, 2004 - Various revision requested by UCSC team
;			implemented by npk@astro.ucsc.edu
;       bugs: still needs overwrite protection on output
;             needs ability to smooth the 2-d atv window?
;
;-
;---------------------

FUNCTION zspec_1d_read, nonlocal, blue=blue, red=red
; Read in the 1d spectrum

  COMMON results_COMMON, results, slitcount, q_index
  COMMON directory_COMMON, fulldir
  COMMON twod_data_COMMON, fullfluxdata, fluxdata, extraction_r1, extraction_r2, extraction_offset, slit_display_options

  specfile = 'spec1d.'+results[slitcount].maskname+'.*' + $
             results[slitcount].objname+'.fits'
  specname = findfile(djs_filepath(specfile, root_dir=fulldir), count=ct)
  IF(ct EQ 0) THEN BEGIN
    mess_txt = 'File not found: ' + specfile
    res = dialog_message(mess_txt,/error)
  ENDIF

; -------- read in data with fill_gap
  IF n_elements(blue) eq 0 then begin			;ADDED BL 7/22/09, gives ability to call the blue response correction if blue is set
	if n_elements(red) eq 0 then begin		;ADDED BL 5/17/10, same thing as blue, but for spectra where lambda_c ~ 8000
		IF(nonlocal EQ 1) THEN $
	    	data = fill_gap(specname[0], /nlsky, /silent, _extra=extra,/tweak) $
	  	ELSE $
	    	data = fill_gap(specname[0], /silent, _extra=extra,/tweak);,  optimal = (boxcar eq 0))
	endif else begin
		IF(nonlocal EQ 1) THEN $
                data = fill_gap_reallyred(specname[0], /nlsky, /silent, _extra=extra,/tweak) $
                ELSE $
                data = fill_gap_reallyred(specname[0], /silent, _extra=extra,/tweak)
	endelse
  endif else begin
	IF(nonlocal EQ 1) THEN $
        data = fill_gap_blue(specname[0], /nlsky, /silent, _extra=extra,/tweak) $
        ELSE $
	data = fill_gap_blue(specname[0], /silent, _extra=extra,/tweak);,  optimal = (boxcar eq 0))
  endelse


  if size(data, /tname) eq 'STRUCT' then begin
  ; -------- correct A band
    hdr = HEADFITS(specname[0], ext=1, /SILENT)
    airmass = sxpar(hdr, 'AIRMASS')
    remove_telluric, data, airmass

  ;------ convert wavelengths from air to vacuum
    lambda = data.lambda
    AIRTOVAC, lambda
    data.lambda = float(lambda)
    extraction_r1 = data.r1
    extraction_r2 = data.r2
  endif else begin
    mess_txt = string(format='(%"%s\n%s\n%s\n%s")', $
    'The fill_gap routine was not able to construct a 1d spectrum', $
    'for this object, possibly because of a bad wavelength solution.', $
    'You should mark this object as Q=2 and enter a comment that',$
    'fill_gap failed for this object.')
    res = dialog_message(mess_txt,/error)
    data = {spec:fltarr(100), $
            lambda:fltarr(100), $
            ivar:fltarr(100), $
            objpos:0.0}
  endelse

  return, data
END



;---------------------------------------------------------------

FUNCTION zspec_1d_eigs
; Reads in the eigenvectors used by the spec1d pipeline to match z's

  deep_dir = getenv('IDLSPEC1D_DIR')
  IF deep_dir EQ '' THEN BEGIN
    mess_txt=string(format='(%"%s\n%s\n%s\n%s\n%s\n%s")', $
                    'The environment variable $IDLSPEC1D_DIR', $
                    'is not set.  This is used to specify the', $
                    'directory which contains the cvs deep IDL',  $
                    'routines.', $
                    'E.g. setenv $IDLSPEC1D_DIR ~dsm/cvs/deep/spec1d', $
                    ' - you should exit and restart IDL.')
    res = dialog_message(mess_txt,/error)
  ENDIF
  eigen_name1 = concat_dir(deep_dir, 'templates/spEigenStarDeep2.fits')
  eigen_name2 = concat_dir(deep_dir, 'templates/spDEEP.fits')
  eigen_name3 = concat_dir(deep_dir, 'templates/spEigenQSOdeep.fits')
  eig1 = mrdfits(eigen_name1, 0, h1, /silent)
  eig2 = mrdfits(eigen_name2, 0, h2, /silent)
  eig3 = mrdfits(eigen_name3, 0, h3, /silent)
; median smooth the absorption-line template spectrum.
  tmpsmth = djs_median(eig2[*,0], width=2500, boundary='reflect')
  npix = n_elements(eig2[*,0])
  loglam = sxpar(h2, 'COEFF0') + findgen(npix) * sxpar(h2, 'COEFF1')
  eig2[*,0] = eig2[*,0] - tmpsmth
; median smooth the K+A template.
  if n_elements(eig2[0,*] ge 3) then begin
      tmpsmth = djs_median(eig2[*,2], width=2500, boundary='reflect')
      eig2[*,2] = eig2[*,2] - tmpsmth
  endif
  eig = {star:eig1, $
         h1:h1, $
         galaxy:eig2, $
         h2:h2, $
         qso:eig3, $
         h3:h3, $
         loglam:loglam, $
         tmpsmth:tmpsmth $
        }

  return, eig

END

;-----------------------------------------------------------------------
function findpix, lambda, lambda0
  diff = abs(lambda - lambda0)
  minval = min(diff, minsub)
  minsub = (minsub - 1) > 0
  return, minsub
end

FUNCTION zspec_template, lambda, eigen

  COMMON results_COMMON, results, slitcount, q_index
  COMMON oned_data_COMMON, specdata, smoothspec1, smoothspec2, synspec

  IF(results[slitcount].nonlocal_select EQ 1) THEN BEGIN
    zans = results[slitcount].zresult_non[results[slitcount].selection]
  ENDIF ELSE BEGIN
    zans = results[slitcount].zresult_all[results[slitcount].selection]
  ENDELSE


  if strtrim(zans.tfile,2) eq '' then begin
     mess_txt = string(format='(%"%s\n%s\n%s")', $
    'Unable to construct a template corresponding to this spectrum.', $
    'Mark this object as Q=2 and add a comment to check that fill_gap', $
    'has actually returned a structure.')
    res = dialog_message(mess_txt,/error)
    return, fltarr(100)
  endif

  ; Determine whether this object is a star or galaxy
  ;--------------------------------------------------
  case strcompress(zans.tfile, /rem) of
      'spEigenStarDeep2.fits': begin
           eig = eigen.star
           h = eigen.h1
       end
       'spEigenQSOdeep.fits':begin
           eig = eigen.qso
           h = eigen.h3
       end
       else: begin
           eig = eigen.galaxy
           h = eigen.h2
       end
   endcase

  w = where(zans.tcolumn GE 0, nw)
  syn = 0

  ; Temporary fix to deal with stars
  IF nw LT 1 THEN BEGIN
    nw = 1
    w = [0]
    zans.tcolumn[w[0]] = 10
  ENDIF

  for i=0,nw-1 do begin
      jj = zans.tcolumn[w[i]]
      syn = syn + zans.theta[i] * eig[*,jj]
  endfor
  IF zans.npoly gt 0 THEN begin ;add in the polynomial terms
      npoints = n_elements(eig[*,0])
      parray = poly_array(npoints, zans.npoly)
      for i=0, zans.npoly-1 do syn = syn + zans.theta[nw+i]*parray[*,i]
  ENDIF
  if strcompress(zans[0].class, /rem) eq 'GALAXY' then begin
      cont = median(specdata.spec)
      npix = n_elements(specdata.lambda)
      restlambda = specdata.lambda / (1.0 + zans.z)
      minpix = findpix(eigen.loglam, alog10(restlambda[0]))
      maxpix = findpix(eigen.loglam, alog10(restlambda[npix-1]))
      tmp_med = median(eigen.tmpsmth[minpix:maxpix])
      syn = syn + eigen.tmpsmth * cont / tmp_med
  endif

  coeff0 = sxpar(h, 'COEFF0')
  coeff1 = sxpar(h, 'COEFF1')

  elambda = 10.d^(coeff0+dindgen((size(eig, /dimens))[0])*coeff1)
; ??? interpol is pretty stupid - use combine1fiber in idlspec2d (SDSS)
  synspec = interpol(syn, elambda, lambda/(zans.z+1))

  return, synspec
END

;------------------------------------------------------------------------

pro zspec_set_mobject, slitnames
  COMMON mobject_COMMON, mobject
  COMMON directory_COMMON, fulldir
; read in the bintab file.
  bfile = findfile(concat_dir(fulldir, '*bintab*.fits*'), count=nbin)
  if nbin eq 0 then message, 'No bintab file found!'
  bfile = bfile[0]
  fits_open, bfile, fcb
  extnames = fcb.extname
  fits_close, fcb
  dex = where(extnames eq 'DesiSlits', cnt)
  if cnt eq 0 then message, 'No DesiSlits Table found!'
  design = mrdfits(bfile, dex[0], /silent)
;  dex = where(extnames eq 'ObjectCat', cnt)
;  if cnt eq 0 then message, 'No ObjectCat Table found!'
;  objcat = mrdfits(bfile, dex[0], /silent)
  dex = where(extnames eq 'SlitObjMap', cnt)
  if cnt eq 0 then message, 'No ObjectCat Table found!'
  slitmap = mrdfits(bfile, dex[0], /silent)
; loop through the slitnames and determine the # of objects in each
; slit.
  nslits = n_elements(slitnames)
  for ii=0,nslits-1 do begin
      dex = where(long(design.slitname) eq long(slitnames[ii]), cnt)
      if cnt eq 0 then message, 'No entry in DesiSlits found!'
      dex = where(slitmap.dslitid eq design[dex[0]].dslitid, cnt)
      if cnt eq 0 then message, 'No entry in SlitObjMap found!'
      if cnt gt 1 then multi = cnt else multi = 1
      if n_elements(mobject) eq 0 then mobject = multi $
      else mobject = [mobject, multi]
  endfor
end
;------------------------------------------------------------------------



;------------------------------------------------------------------------
;------------------------------------------------------------------------

;                   Sub routines






PRO zspec_1d_plot, state, redshift

  COMMON oned_data_COMMON, specdata, smoothspec1, smoothspec2, synspec

  device, decomposed = 0
  tvlct, r, g, b, /get

  rtiny   = [0, 1, 0, 0, 0, 1, 1, 1]
  gtiny = [0, 0, 1, 0, 1, 0, 1, 1]
  btiny  = [0, 0, 0, 1, 1, 1, 0, 1]
  nt = n_elements(rtiny)-1

  rtiny   = [0, 1, 0, 0, 0, 1, 1, 1]
  gtiny = [0, 0, 1, 0, 1, 0, 1, 1]
  btiny  = [0, 0, 0, 1, 1, 1, 0, 1]
  tvlct, 255*rtiny, 255*gtiny, 255*btiny
  tvlct, [255], [255], [255], !d.table_size-1


  rlambda = specdata.lambda/(1.+redshift)
  IF(state.smooth_on EQ 1) THEN $
    ssyn = ivarsmooth(synspec,specdata.ivar,state.smooth1)
  IF(state.smooth_on EQ 2) THEN $
    ssyn = ivarsmooth(synspec,specdata.ivar,state.smooth2)

; Counter to increment which plot we are to draw
  nplot = 0

  ; OII emission feature
  ;-----------------------
;  wset, state.draw_1d_idx[0]
  wset, state.draw_1d_idx[nplot]
  cnt = WHERE(rlambda LE 3747 AND rlambda GE 3707, c)
  IF c GE 1 THEN BEGIN
    CASE state.smooth_on OF
      '0':  BEGIN
         plot, specdata.lambda[cnt], specdata.spec[cnt], $
              title='!6[OII]', charsize=0.8, background=7, color=0
         oplot, specdata.lambda[cnt], synspec[cnt], color=1, thick=2
     END
     '1': BEGIN
         plot, specdata.lambda[cnt], smoothspec1[cnt], $
               title='!6[OII]', charsize=0.8, background=7, color=0
         oplot, specdata.lambda[cnt], ssyn[cnt], color=1, thick=2
     END
     '2': BEGIN
         plot, specdata.lambda[cnt], smoothspec2[cnt], $
               title='!6[OII]', charsize=0.8, background=7, color=0
         oplot, specdata.lambda[cnt], ssyn[cnt], color=1, thick=2
     END
    ENDCASE
    nplot = nplot + 1
   ENDIF
;  ENDIF ELSE $
;   plot, [0,1], [0,1], charsize=0.8, background=7, color=0


  ; H and K absorption features
  ;-----------------------------
;  wset, state.draw_1d_idx[1]
  wset, state.draw_1d_idx[nplot]
  cnt = WHERE(rlambda LE 4000 AND rlambda GE 3900, c)
  IF c GE 1 THEN BEGIN
    CASE state.smooth_on OF
      '0':  BEGIN
         plot, specdata.lambda[cnt], specdata.spec[cnt], $
              title='H & K', charsize=0.8, background=7, color=0
         oplot, specdata.lambda[cnt], synspec[cnt], color=1, thick=2
     END
     '1': BEGIN
         plot, specdata.lambda[cnt], smoothspec1[cnt], $
               title='H & K', charsize=0.8, background=7, color=0
         oplot, specdata.lambda[cnt], ssyn[cnt], color=1, thick=2
     END
     '2': BEGIN
         plot, specdata.lambda[cnt], smoothspec2[cnt], $
               title='H & K', charsize=0.8, background=7, color=0
         oplot, specdata.lambda[cnt], ssyn[cnt], color=1, thick=2
     END
    ENDCASE
    nplot = nplot + 1
   ENDIF
;  ENDIF ELSE $
;   plot, [0,1], [0,1], charsize=0.8, background=7, color=0



; IF redshift GT 0.3 THEN BEGIN

  ; Hbeta feature
  ;----------------
;  wset, state.draw_1d_idx[2]
  wset, state.draw_1d_idx[nplot]
  cnt = WHERE(rlambda LE 4880 AND rlambda GE 4840, c)
  IF c GE 1 THEN BEGIN
    CASE state.smooth_on OF
      '0':  BEGIN
         plot, specdata.lambda[cnt], specdata.spec[cnt], $
              title='!6H!7b!6', charsize=0.8, background=7, color=0
         oplot, specdata.lambda[cnt], synspec[cnt], color=1, thick=2
     END
     '1': BEGIN
         plot, specdata.lambda[cnt], smoothspec1[cnt], $
               title='!6H!7b!6', charsize=0.8, background=7, color=0
         oplot, specdata.lambda[cnt], ssyn[cnt], color=1, thick=2
     END
     '2': BEGIN
         plot, specdata.lambda[cnt], smoothspec2[cnt], $
               title='!6H!7b!6', charsize=0.8, background=7, color=0
         oplot, specdata.lambda[cnt], ssyn[cnt], color=1, thick=2
     END
    ENDCASE
    nplot = nplot + 1
  ENDIF
;  ENDIF ELSE $
;   plot, [0,1], [0,1], charsize=0.8, background=7, color=0


  ; O[III] emission feature
  ;--------------------------
;  wset, state.draw_1d_idx[3]
  wset, state.draw_1d_idx[nplot]
  cnt = WHERE(rlambda LE 5030 AND rlambda GE 4940, c)
  IF c GE 1 THEN BEGIN
    CASE state.smooth_on OF
      '0':  BEGIN
         plot, specdata.lambda[cnt], specdata.spec[cnt], $
              title='[OIII]', charsize=0.8, background=7, color=0
         oplot, specdata.lambda[cnt], synspec[cnt], color=1, thick=2
     END
     '1': BEGIN
         plot, specdata.lambda[cnt], smoothspec1[cnt], $
               title='[OIII]', charsize=0.8, background=7, color=0
         oplot, specdata.lambda[cnt], ssyn[cnt], color=1, thick=2
     END
     '2': BEGIN
         plot, specdata.lambda[cnt], smoothspec2[cnt], $
               title='[OIII]', charsize=0.8, background=7, color=0
         oplot, specdata.lambda[cnt], ssyn[cnt], color=1, thick=2
     END
    ENDCASE
    nplot = nplot + 1
  ENDIF
;  ENDIF ELSE $
;   plot, [0,1], [0,1], charsize=0.8, background=7, color=0


; ENDIF ELSE BEGIN

  ; Halpha feature
  ;-----------------
; now we have to check to make sure we don't draw >4 plots
  IF nplot LE 3 THEN BEGIN
;  wset, state.draw_1d_idx[2]
  wset, state.draw_1d_idx[nplot]
  cnt = WHERE(rlambda LE 6600 AND rlambda GE 6540, c)
  IF c GE 1 THEN BEGIN
    CASE state.smooth_on OF
      '0':  BEGIN
         plot, specdata.lambda[cnt], specdata.spec[cnt], $
              title='!6H!7a!6', charsize=0.8, background=7, color=0
         oplot, specdata.lambda[cnt], synspec[cnt], color=1,thick=2
     END
     '1': BEGIN
         plot, specdata.lambda[cnt], smoothspec1[cnt], $
               title='!6H!7a!6', charsize=0.8, background=7, color=0
         oplot, specdata.lambda[cnt], ssyn[cnt], color=1,thick=2
     END
     '2': BEGIN
         plot, specdata.lambda[cnt], smoothspec2[cnt], $
               title='!6H!7a!6', charsize=0.8, background=7, color=0
         oplot, specdata.lambda[cnt], ssyn[cnt], color=1,thick=2
     END
    ENDCASE
    nplot = nplot + 1
  ENDIF
  ENDIF
;  ENDIF ELSE $
;   plot, [0,1], [0,1], charsize=0.8, background=7, color=0

  ; S emission feature
  ;---------------------
  IF nplot LE 3 THEN BEGIN
;  wset, state.draw_1d_idx[3]
  wset, state.draw_1d_idx[nplot]
  cnt = WHERE(rlambda LE 6745 AND rlambda GE 6700, c)
  IF c GE 1 THEN BEGIN
    CASE state.smooth_on OF
      '0':  BEGIN
         plot, specdata.lambda[cnt], specdata.spec[cnt], $
              title='[SII]', charsize=0.8, background=7, color=0
         oplot, specdata.lambda[cnt], synspec[cnt], color=1, thick=2
     END
     '1': BEGIN
         plot, specdata.lambda[cnt], smoothspec1[cnt], $
               title='[SII]', charsize=0.8, background=7, color=0
         oplot, specdata.lambda[cnt], ssyn[cnt], color=1, thick=2
     END
     '2': BEGIN
         plot, specdata.lambda[cnt], smoothspec2[cnt], $
               title='[SII]', charsize=0.8, background=7, color=0
         oplot, specdata.lambda[cnt], ssyn[cnt], color=1, thick=2
     END
    ENDCASE
    nplot = nplot + 1
  ENDIF
  ENDIF
;  ENDIF ELSE $
;   plot, [0,1], [0,1], charsize=0.8, background=7, color=0

;  ENDELSE

; fill up the remaining plots, if any, to make sure no old ones are left

  for ii=nplot,3 do begin
    wset, state.draw_1d_idx[ii]
    plot, [0,1], [0,1], charsize=0.8, background=7, color=0
  endfor

END








PRO zspec_splot, state, redshift, first

  COMMON zspec_lines_COMMON, linenames, linewaves, telluricnames, telluricwaves
  COMMON oned_data_COMMON, specdata, smoothspec1, smoothspec2, synspec

  rlambda = specdata.lambda/(1.+redshift)
  IF(state.smooth_on EQ 1) THEN $
    ssyn = ivarsmooth(synspec,specdata.ivar,state.smooth1)
  IF(state.smooth_on EQ 2) THEN $
    ssyn = ivarsmooth(synspec,specdata.ivar,state.smooth2)



 ; Splot routines.
 ; I have taken these straight from pspec2...
 ;---------------------------------------------
  medspec = median(specdata.spec)
  stdspec = stddev(specdata.spec)
  maxs = medspec + 3.*stdspec ;set limits not so affected by bad points
  mins = medspec - 3.*stdspec
  mins_new = mins - 0.1*(maxs-mins)

  lambdas = specdata.lambda
  VACTOAIR, lambdas

  CASE state.smooth_on OF
    '0':  BEGIN
       splot, lambdas, specdata.spec, $
         xtitle='!6Observed Wavelength (Ang)', $
         ytitle='Counts/hour',  yrange=[mins_new,maxs]
       soplot, lambdas, synspec, color=4
   END
   '1': BEGIN
       splot, lambdas, smoothspec1, $
         xtitle='!6Observed Wavelength (Ang)', $
         ytitle='Counts/hour',  yrange=[mins_new,maxs]
       soplot, lambdas, ssyn, color=4
   END
   '2': BEGIN
       splot, lambdas, smoothspec2, $
         xtitle='!6Observed Wavelength (Ang)', $
         ytitle='Counts/hour',  yrange=[mins_new,maxs]
       soplot, lambdas, ssyn, color=4
   END
  ENDCASE

  variance = specdata.ivar*0.
  good = where(specdata.ivar gt 0.)
  variance[good] = 1./specdata.ivar[good] <  1.e5
  maxv = MAX(variance)
  minv = MIN(variance)
  variance = 0.5*(maxs-mins)/(maxv-minv)*variance + mins
  SOPLOT, lambdas, variance, color=1
  min_lambda = MIN(lambdas)
  max_lambda = MAX(lambdas)


  y = [mins_new,maxs]
  for i = 0, n_elements(linewaves)-1 do begin
  	l = linewaves[i]*(1.+redshift)
	if l GT min_lambda AND l LT max_lambda then begin
		soplot, [l, l], y, linestyle=2,color=6
		sxyouts, l, 0.8 * maxs, linenames[i]
	endif
	if linenames[i] EQ 'OII' then i = i+2
  endfor

  for i = 0, n_elements(telluricwaves)-1 do begin
  	l = telluricwaves[i]
	if l GT min_lambda and l LT max_lambda then begin
		soplot, [l, l], y, linestyle=2,color=5,thick=3
		sxyouts, l, 0.9 * maxs, telluricnames[i]
	endif
  endfor


  IF(NOT arg_present(first)) THEN widget_control, state.base, /show


END

FUNCTION zspec_get_redshift
	COMMON results_COMMON, results, slitcount
	IF(results[slitcount].nonlocal_select EQ 1) THEN BEGIN
	      redshift = $
		results[slitcount].zresult_non[results[slitcount].selection].z
	ENDIF ELSE BEGIN
	       ; Normal sky subtraction options
	       redshift = $
	        results[slitcount].zresult_all[results[slitcount].selection].z
	ENDELSE
	return, redshift
END


PRO zspec_splot_show, event

  COMMON results_COMMON, results, slitcount, q_index

  widget_control, event.top, get_uvalue=state

  IF(results[slitcount].nonlocal_select EQ 1) THEN BEGIN
    redshift = $
             results[slitcount].zresult_non[results[slitcount].selection].z
  ENDIF ELSE BEGIN
  ; Normal sky subtraction options
    redshift = $
            results[slitcount].zresult_all[results[slitcount].selection].z
  ENDELSE

  first = 1
  zspec_splot, state, redshift, first
END



;------------------------------------------------------------------------

PRO zspec_2d_plot, state, redshift

  COMMON oned_data_COMMON, specdata, smoothspec1, smoothspec2, synspec
  COMMON twod_data_COMMON, fullfluxdata, fluxdata, extraction_r1, extraction_r2, extraction_offset, slit_display_options

  if(size(fluxdata, /tname) ne 'STRUCT') then return

; Counter to increment which plot we are to draw
; One has to hope that the logic chooses the same set of lines for 2d
; as it did in zspec_1d_plot; it should ...
  nplot = 0

  ; [OII]
  ;--------
;  wset, state.draw_2d_idx[0]
  wset, state.draw_2d_idx[nplot]

  lam_min  = (1.+redshift) * 3707
  lam_max  = (1.+redshift) * 3747
  oii_blam = where(fluxdata.blambda[*,0] gt lam_min and $
                   fluxdata.blambda[*,0] lt lam_max,nb)
  oii_rlam = where(fluxdata.rlambda[*,0] gt lam_min and $
                   fluxdata.rlambda[*,0] lt lam_max,nr)

  ; figure out which side, red or blue, is the one with OII
  side = 0
  IF (nr GT 1 AND nr GT nb) THEN BEGIN
    nlam1 = FINDGEN(200) * N_ELEMENTS(fluxdata.rflux[oii_rlam,0])/200
    nlam2 = FINDGEN(50)
    flux = INTERPOLATE(fluxdata.rflux[oii_rlam,*],nlam1,nlam2,/grid)
    fullrange = minmax(flux)

    if fullrange[1] lt (-5) then minscl = fullrange[0] else minscl =  (-5)
    if fullrange[0] gt (25) then maxscl = fullrange[1] else maxscl =  (25)

    tv, bytscl(flux, min = minscl, max = maxscl, top = 255b-8b) +8b
    nplot = nplot + 1
  ENDIF ELSE IF (nb GT 1 AND nb GT nr) THEN BEGIN
    nlam1 = FINDGEN(200) * N_ELEMENTS(fluxdata.bflux[oii_blam,0])/200
    nlam2 = FINDGEN(50)
    flux = INTERPOLATE(fluxdata.bflux[oii_blam,*],nlam1,nlam2,/grid)
    fullrange = minmax(flux)

    if fullrange[1] lt (-5) then minscl = fullrange[0] else minscl =  (-5)
    if fullrange[0] gt (25) then maxscl = fullrange[1] else maxscl =  (25)

    tv, bytscl(flux, min = minscl, max = maxscl, top = 255b-8b) +8b
    nplot = nplot + 1
  ENDIF
;  ENDIF ELSE $
;    tv, bytarr(200,50)

  ; H&K
  ;-------
;  wset, state.draw_2d_idx[1]
  wset, state.draw_2d_idx[nplot]
  lam_min  = (1.+redshift) * 3900
  lam_max  = (1.+redshift) * 4000
  oii_blam = where(fluxdata.blambda[*,0] gt lam_min and $
                   fluxdata.blambda[*,0] lt lam_max,nb)
  oii_rlam = where(fluxdata.rlambda[*,0] gt lam_min and $
                   fluxdata.rlambda[*,0] lt lam_max,nr)

  ; figure out which side, red or blue, is the one with OII
  side = 0
  IF (nr GT 1 AND nr GT nb) THEN BEGIN
    nlam1 = FINDGEN(200) * N_ELEMENTS(fluxdata.rflux[oii_rlam,0])/200
    nlam2 = FINDGEN(50)
    flux = INTERPOLATE(fluxdata.rflux[oii_rlam,*],nlam1,nlam2,/grid)
    fullrange = minmax(flux)

    if fullrange[1] lt (-5) then minscl = fullrange[0] else minscl =  (-5)
    if fullrange[0] gt (25) then maxscl = fullrange[1] else maxscl =  (25)

    tv, bytscl(flux, min = minscl, max = maxscl, top = 255b-8b) +8b
    nplot = nplot + 1
  ENDIF ELSE IF (nb GT 1 AND nb GT nr) THEN BEGIN
    nlam1 = FINDGEN(200) * N_ELEMENTS(fluxdata.bflux[oii_blam,0])/200
    nlam2 = FINDGEN(50)
    flux = INTERPOLATE(fluxdata.bflux[oii_blam,*],nlam1,nlam2,/grid)
    fullrange = minmax(flux)

    if fullrange[1] lt (-5) then minscl = fullrange[0] else minscl =  (-5)
    if fullrange[0] gt (25) then maxscl = fullrange[1] else maxscl =  (25)

    tv, bytscl(flux, min = minscl, max = maxscl, top = 255b-8b) +8b
    nplot = nplot + 1
  ENDIF
;  ENDIF ELSE $
;    tv, bytarr(200,50)



; wset, state.draw_2d_idx[2]
; tv, bytarr(200,50)
; wset, state.draw_2d_idx[3]
; tv, bytarr(200,50)

; IF redshift GT 0.3 THEN BEGIN
  ; Hbeta
  ;-------
;  wset, state.draw_2d_idx[2]
  wset, state.draw_2d_idx[nplot]
  lam_min  = (1.+redshift) * 4840
  lam_max  = (1.+redshift) * 4880
  oii_blam = where(fluxdata.blambda[*,0] gt lam_min and $
                   fluxdata.blambda[*,0] lt lam_max,nb)
  oii_rlam = where(fluxdata.rlambda[*,0] gt lam_min and $
                   fluxdata.rlambda[*,0] lt lam_max,nr)

  ; figure out which side, red or blue, is the one with the feature
  side = 0
  IF (nr GT 1 AND nr GT nb) THEN BEGIN
    nlam1 = FINDGEN(200) * N_ELEMENTS(fluxdata.rflux[oii_rlam,0])/200
    nlam2 = FINDGEN(50)
    flux = INTERPOLATE(fluxdata.rflux[oii_rlam,*],nlam1,nlam2,/grid)
    fullrange = minmax(flux)

    if fullrange[1] lt (-5) then minscl = fullrange[0] else minscl =  (-5)
    if fullrange[0] gt (25) then maxscl = fullrange[1] else maxscl =  (25)

    tv, bytscl(flux, min = minscl, max = maxscl, top = 255b-8b) +8b
    nplot = nplot + 1
  ENDIF ELSE IF (nb GT 1 AND nb GT nr) THEN BEGIN
    nlam1 = FINDGEN(200) * N_ELEMENTS(fluxdata.bflux[oii_blam,0])/200
    nlam2 = FINDGEN(50)
    flux = INTERPOLATE(fluxdata.bflux[oii_blam,*],nlam1,nlam2,/grid)
    fullrange = minmax(flux)

    if fullrange[1] lt (-5) then minscl = fullrange[0] else minscl =  (-5)
    if fullrange[0] gt (25) then maxscl = fullrange[1] else maxscl =  (25)

    tv, bytscl(flux, min = minscl, max = maxscl, top = 255b-8b) +8b
    nplot = nplot + 1
  ENDIF


  ; [OIII]
  ;-------
;  wset, state.draw_2d_idx[3]
  wset, state.draw_2d_idx[nplot]

  lam_min  = (1.+redshift) * 4940
  lam_max  = (1.+redshift) * 5030
  oii_blam = where(fluxdata.blambda[*,0] gt lam_min and $
                   fluxdata.blambda[*,0] lt lam_max,nb)
  oii_rlam = where(fluxdata.rlambda[*,0] gt lam_min and $
                   fluxdata.rlambda[*,0] lt lam_max,nr)

  ; figure out which side, red or blue, is the one with OIII
  side = 0
  IF (nr GT 1 AND nr GT nb) THEN BEGIN
    nlam1 = FINDGEN(200) * N_ELEMENTS(fluxdata.rflux[oii_rlam,0])/200
    nlam2 = FINDGEN(50)
    flux = INTERPOLATE(fluxdata.rflux[oii_rlam,*],nlam1,nlam2,/grid)
    fullrange = minmax(flux)

    if fullrange[1] lt (-5) then minscl = fullrange[0] else minscl =  (-5)
    if fullrange[0] gt (25) then maxscl = fullrange[1] else maxscl =  (25)

    tv, bytscl(flux, min = minscl, max = maxscl, top = 255b-8b) +8b
    nplot = nplot + 1
  ENDIF ELSE IF (nb GT 1 AND nb GT nr) THEN BEGIN
    nlam1 = FINDGEN(200) * N_ELEMENTS(fluxdata.bflux[oii_blam,0])/200
    nlam2 = FINDGEN(50)
    flux = INTERPOLATE(fluxdata.bflux[oii_blam,*],nlam1,nlam2,/grid)
    fullrange = minmax(flux)

    if fullrange[1] lt (-5) then minscl = fullrange[0] else minscl =  (-5)
    if fullrange[0] gt (25) then maxscl = fullrange[1] else maxscl =  (25)

    tv, bytscl(flux, min = minscl, max = maxscl, top = 255b-8b) +8b
    nplot = nplot + 1
  ENDIF

; ENDIF ELSE IF(redshift GT 0.0 AND redshift LT 0.3) THEN BEGIN

  ; Halpha
  ;--------
; now we have to check to make sure we don't draw >4 plots
  IF nplot LE 3 THEN BEGIN
;  wset, state.draw_2d_idx[2]
  wset, state.draw_2d_idx[nplot]
  lam_min  = (1.+redshift) * 6540
  lam_max  = (1.+redshift) * 6600
  oii_blam = where(fluxdata.blambda[*,0] gt lam_min and $
                   fluxdata.blambda[*,0] lt lam_max,nb)
  oii_rlam = where(fluxdata.rlambda[*,0] gt lam_min and $
                   fluxdata.rlambda[*,0] lt lam_max,nr)

  ; figure out which side, red or blue, is the one with Ha
  side = 0
  IF (nr GT 1 AND nr GT nb) THEN BEGIN
    nlam1 = FINDGEN(200) * N_ELEMENTS(fluxdata.rflux[oii_rlam,0])/200
    nlam2 = FINDGEN(50)
    flux = INTERPOLATE(fluxdata.rflux[oii_rlam,*],nlam1,nlam2,/grid)
    fullrange = minmax(flux)

    if fullrange[1] lt (-5) then minscl = fullrange[0] else minscl =  (-5)
    if fullrange[0] gt (25) then maxscl = fullrange[1] else maxscl =  (25)

    tv, bytscl(flux, min = minscl, max = maxscl, top = 255b-8b) +8b
    nplot = nplot + 1
  ENDIF ELSE IF (nb GT 1 AND nb GT nr) THEN BEGIN
    nlam1 = FINDGEN(200) * N_ELEMENTS(fluxdata.bflux[oii_blam,0])/200
    nlam2 = FINDGEN(50)
    flux = INTERPOLATE(fluxdata.bflux[oii_blam,*],nlam1,nlam2,/grid)
    fullrange = minmax(flux)

    if fullrange[1] lt (-5) then minscl = fullrange[0] else minscl =  (-5)
    if fullrange[0] gt (25) then maxscl = fullrange[1] else maxscl =  (25)

    tv, bytscl(flux, min = minscl, max = maxscl, top = 255b-8b) +8b
    nplot = nplot + 1
  ENDIF
  ENDIF


  ; [SII]
  ;-------
  IF nplot LE 3 THEN BEGIN
;  wset, state.draw_2d_idx[3]
  wset, state.draw_2d_idx[nplot]

  lam_min  = (1.+redshift) * 6700
  lam_max  = (1.+redshift) * 6745
  oii_blam = where(fluxdata.blambda[*,0] gt lam_min and $
                   fluxdata.blambda[*,0] lt lam_max,nb)
  oii_rlam = where(fluxdata.rlambda[*,0] gt lam_min and $
                   fluxdata.rlambda[*,0] lt lam_max,nr)

  ; figure out which side, red or blue, is the one with OII
  side = 0
  IF (nr GT 1 AND nr GT nb) THEN BEGIN
    nlam1 = FINDGEN(200) * N_ELEMENTS(fluxdata.rflux[oii_rlam,0])/200
    nlam2 = FINDGEN(50)
    flux = INTERPOLATE(fluxdata.rflux[oii_rlam,*],nlam1,nlam2,/grid)
    fullrange = minmax(flux)

    if fullrange[1] lt (-5) then minscl = fullrange[0] else minscl =  (-5)
    if fullrange[0] gt (25) then maxscl = fullrange[1] else maxscl =  (25)

    tv, bytscl(flux, min = minscl, max = maxscl, top = 255b-8b) +8b
    nplot = nplot + 1
  ENDIF ELSE IF (nb GT 1 AND nb GT nr) THEN BEGIN
    nlam1 = FINDGEN(200) * N_ELEMENTS(fluxdata.bflux[oii_blam,0])/200
    nlam2 = FINDGEN(50)
    flux = INTERPOLATE(fluxdata.bflux[oii_blam,*],nlam1,nlam2,/grid)
    fullrange = minmax(flux)

    if fullrange[1] lt (-5) then minscl = fullrange[0] else minscl =  (-5)
    if fullrange[0] gt (25) then maxscl = fullrange[1] else maxscl =  (25)

    tv, bytscl(flux, min = minscl, max = maxscl, top = 255b-8b) +8b
    nplot = nplot + 1
  ENDIF
  ENDIF

; ENDIF

; fill up the remaining plots, if any, to make sure no old ones are left

  for ii=nplot,3 do begin
    wset, state.draw_2d_idx[ii]
    tv, bytarr(200,50)
  endfor

END


;-------------------------------------------------------------------------




;-----------------------------------------------------------------------

FUNCTION untilt_slit, spec, lambda
	nrows = n_elements(spec.flux[0,*])
	for i = 0, nrows - 1 do begin
		spec.flux[*,i] = interpol(spec.flux[*,i], lambda[*,i], $
			spec.lambda0)
		lambda[*,i] = spec.lambda0
	endfor
	
	return, spec 
END
	
FUNCTION zspec_2d_read, spec1d, nonlocal, fullspec=fullspec
; Read in the 2d spectrum

  COMMON results_COMMON, results, slitcount, q_index
  COMMON directory_COMMON, fulldir
  COMMON twod_data_COMMON, fullfluxdata, fluxdata, extraction_r1, extraction_r2, extraction_offset, slit_display_options

  fname1 = 'slit.' + strcompress(results[slitcount].maskname, /remove_all) $
           + '.' + strcompress(results[slitcount].slitname, /remove_all) $
           + 'B.fits.gz'
  fname1 = concat_dir(fulldir,fname1)
  fname2 = 'slit.' + strcompress(results[slitcount].maskname, /remove_all) $
           + '.' + strcompress(results[slitcount].slitname, /remove_all) $
           + 'R.fits.gz'
  fname2 = concat_dir(fulldir,fname2)

  IF(nonlocal EQ 1) THEN BEGIN
    slit_blue = mrdfits(fname1, 3, /silent,status=status1)
    slit_red  = mrdfits(fname2, 3, /silent,status=status2)
  ENDIF ELSE BEGIN
    slit_blue = mrdfits(fname1, 1, /silent,status=status1)
    slit_red  = mrdfits(fname2, 1, /silent,status=status2)
  ENDELSE
  if status1 ne -1 then lambda_blue = lambda_eval(slit_blue) $
  else lambda_blue = lambda_eval(slit_red) * 0.0
  if status2 ne -1 then lambda_red  = lambda_eval(slit_red) $
  else lambda_red = lambda_eval(slit_blue) * 0.0


  ; Do a quick check to make sure things have actually worked
  ;------------------------------------------------------
  IF (status1 EQ -1) and (status2 eq -1) THEN BEGIN
    mess_txt = 'No slit files found for ' + $
      strcompress(results[slitcount].slitname, /rem) + '!'
    res = dialog_message(mess_txt, /error)
    return, 0
  ENDIF
  if (status1 eq -1) then begin
      mess_txt = fname1 + ' not found!'
      res = dialog_message(mess_txt, /error)
      slit_blue = {flux:(slit_red.flux * 0.0)}
  endif
  if (status2 eq -1) then begin
      mess_txt = fname2 + ' not found!'
      res = dialog_message(mess_txt, /error)
      slit_red = {flux:(slit_blue.flux * 0.0)}
  endif

  dim1 = size(slit_blue.flux, /dimensions)
  dim2 = size(slit_red.flux, /dimensions)
  width1 = dim1[1]
  width2 = dim2[1]

  objpos = fix(spec1d.objpos)
  if slit_display_options.untilt then begin
  	slit_blue = untilt_slit(slit_blue, lambda_blue)
  	slit_red = untilt_slit(slit_red, lambda_red)
  endif

  if keyword_set(fullspec) eq 1 then begin
  	extraction_offset = 0
  	data = {blambda:lambda_blue, $
  		rlambda:lambda_red, $
  		bflux: slit_blue.flux, $
  		rflux: slit_red.flux}
  	return, data
  endif


  ; Blue slit--
  extraction_offset = 0
  IF width1 GT 50 THEN BEGIN
      IF objpos GT 25 AND (objpos+25) LT width1 THEN BEGIN
          width1 = 50
          ext1_1 = objpos-25
          ext2_1 = objpos+25
      ENDIF ELSE IF objpos GT 25 AND (objpos+25) GE width1 THEN BEGIN
          ext2_1 = width1
          ext1_1 = width1-50
          width1 = 50
      ENDIF ELSE BEGIN
          width1 = 50
          ext1_1 = 0
          ext2_1 = width1
      ENDELSE
  ENDIF ELSE BEGIN
      ext1_1 = 0
      ext2_1 = 2*fix(width1/2)
  ENDELSE


  ; Red slit
  IF width2 GT 50 THEN BEGIN
      IF objpos GT 25 AND (objpos+25) LT width2 THEN BEGIN
          width2 = 50
          ext1_2 = objpos-25
          ext2_2 = objpos+25
      ENDIF ELSE IF objpos GT 25 AND (objpos+25) GE width2 THEN BEGIN
          ext2_2 = width2
          ext1_2 = width2-50
          width2 = 50
      ENDIF ELSE BEGIN
          width2 = 50
          ext1_2 = 0
          ext2_2 = width2
      ENDELSE
  ENDIF ELSE BEGIN
      ext1_2 = 0
      ext2_2 = 2*fix(width2/2)
  endelse


  flux_blue = fltarr(dim1[0],50)
  flux_red  = fltarr(dim2[0],50)
  start1  = 25-fix(width1/2)
  finish1 = 25+fix(width1/2)
  start2  = 25-fix(width2/2)
  finish2 = 25+fix(width2/2)


  flux_blue[*,start1:finish1-1] = slit_blue.flux[*,ext1_1:ext2_1-1]
  flux_red[*,start2:finish2-1]  = slit_red.flux[*,ext1_2:ext2_2-1]

  if ext1_1 ne ext1_2 then begin
  	print, "The 1d extraction begins at ", string(ext1_1), " on the red ", $
		"side, but at ", string(ext1_2), " on the blue side."
  endif
  if status1 ne -1 then lambda_blue = lambda_eval(slit_blue) $
  else lambda_blue = lambda_eval(slit_red) * 0.0
  if status2 ne -1 then lambda_red  = lambda_eval(slit_red) $
  else lambda_red = lambda_eval(slit_blue) * 0.0



  extraction_offset = start1
  data = {blambda:lambda_blue, $
          rlambda:lambda_red, $
          bflux:flux_blue, $
          rflux:flux_red $
         }

  return, data
END


;-----------------------------------------------------------------------

; Re-scan z:   M.Cooper  (12/4/03)

; ss1d = 1-d object spectrum structure containing tags for flux, ivar,
; and lambda, etc..
; z = the z value that the user suggests as correct.
; returns a zresult structure with the new redshift as found by
; zrefind.
; zold = a full zresult structure for the object...from which to grab
; the maskname, slitname, eigenfile, etc.
; OUTPUT: a new zresult structure equipped with new z value and new
; vdisp value and corresponding errors.


function zspec_zfind, z, zold, ss1d,  zmin=zmin, zmax=zmax, $
                      wmin=wmin, wmax=wmax, $
                      columns=columns, eigenfile=eigenfile, $
                      star=star, qso=qso

; make some basic definitions for parameters to be passed to
; zfind.pro.
  pspace = 1
  width = 3.0 * pspace
;  zguess = z
;  zmin = zguess - 0.05
;  zmax = zguess + 0.05
  pwidth = 15
; use the same template file and columns (eigenspectra) from that
; file.
  if n_elements(eigenfile) gt 0 then eigenfile = eigenfile[0] $
  else eigenfile = 'spDEEP.fits'
; append the needed path.
  spec1d_dir = getenv('IDLSPEC1D_DIR')
  eigendir = concat_dir(spec1d_dir, 'templates')

; don't just take the old tcolumn info.....try fitting with both templates.
  if n_elements(columns) gt 0 then columns = columns $
  else columns =  [0,1,2]

  if n_elements(star) gt 0 then star = star ge 1 else star = 0
  if n_elements(qso) gt 0 then qso = qso ge 1 else qso = 0

; use the same type of polynomial fitting: npoly defines the number
; of polynomials to fit to the continuum.
  npoly = zold.npoly
; define the number of redshifts to find.
  nfind = 1

; find the redshift value...
  if star then begin
      shdr = headfits(concat_dir(eigendir, eigenfile))
      nstar = sxpar(shdr, 'NAXIS2') > 1
      subclass = strarr(nstar)
      for istar=0,nstar-1 do $
        subclass[istar] = $
        strtrim( sxpar(shdr, 'NAME'+strtrim(string(istar),2)), 2)
      zresult = zfind_star(ss1d, /linear_lambda, $
                           eigenfile=eigenfile, eigendir=eigendir, $
                           npoly=npoly, subclass=subclass, $
                           zmin=zmin, zmax=zmax, pspace=pspace, $
                           nfind=nfind, width=5*pspace, $
                           wvmin=wvmin, wvmax=wvmax)
      zresult = zresult[0]
  endif

  if qso then $
    zresult = zfind_qso(ss1d, /linear_lambda, eigendir=eigendir, $
                        eigenfile=eigenfile, npoly=npoly, $
                        zmin=zmin, zmax=zmax, pspace=pspace, $
                        nfind=nfind, width=5*pspace, $
                        wvmin=wvmin, wvmax=wvmax)

  if qso eq 0 and star eq 0 then $
    zresult = zfind(ss1d, eigenfile=eigenfile, eigendir=eigendir, $
                    columns=columns, npoly=npoly, $
                    zmin=zmin, zmax=zmax, $
                    nfind=nfind, pspace=pspace, $
                    width=5*pspace, /linear_lambda, $
                    wvmin=wmin, wvmax=wmax)

  if columns[0] eq 1 then begin
      theta = fltarr(10)
      theta[1] = zresult.theta[0]
      zresult.theta = theta
      tcolumn = lonarr(10)
      tcolumn[1] = zresult.tcolumn[0]
      zresult.tcolumn = tcolumn
  endif

  if columns[0] eq 2 then begin
      theta = fltarr(10)
      theta[2] = zresult.theta[0]
      zresult.theta = theta
      tcolumn = lonarr(10)
      tcolumn[2] = zresult.tcolumn[0]
      zresult.tcolumn = tcolumn
  endif


; now redo the vdisp determination.
  vdispfit, ss1d, zobj=zresult.z, sigma=sigma, sigerr=sigerr
  zresult.vdisp = sigma
  zresult.vdisp_err = sigerr

; copy over some of the info from the old zresult structure.
  zresult.maskname = zold.maskname
  zresult.slitname = zold.slitname
  zresult.objname = zold.objname
  zresult.date = zresult.date


; make sure to set the rchi2diff value to 0.0; it is now meaningless.
  zresult.rchi2diff = 0.0
;  zresult.comment = ''
;  zresult.zquality = 0B
;  zresult.rotcurve = 0B

; return the new zresult structure.
  return, zresult

end



;------------------------------------------------------------------------
;------------------ Event based routines ---------------------------

PRO zspec_selectz, event

  COMMON results_COMMON, results, slitcount, q_index
  COMMON oned_data_COMMON, specdata, smoothspec1, smoothspec2, synspec
  widget_control, event.top, get_uvalue=state

  results[slitcount].selection = event.index
  synspec  = zspec_template(specdata.lambda,state.eigenvec)

;  results[slitcount].zresult =
;  results[slitcount].zresult_all[event.index]
  if(results[slitcount].nonlocal_select eq 1) then $
    redshift = results[slitcount].zresult_non[results[slitcount].selection].z
  if(results[slitcount].nonlocal_select eq 0) then $
    redshift = results[slitcount].zresult_all[results[slitcount].selection].z

  zspec_1d_plot, state, redshift
  IF(xregistered('splot')) THEN zspec_splot, state, redshift
  zspec_2d_plot, state, redshift
  zspec_showtext, state
  if (xregistered('atv')) then zspec_update_atv, state

END


;----------------------------------------------------------------



FUNCTION zspec_set_smooth_event, event
  COMMON zspec_set_smooth_COMMON, selection

  widget_control, event.ID, get_uvalue=uval

  IF(uval EQ -1) THEN  $
    widget_control, event.top, /destroy $
  ELSE $
    selection = event.value

  return, 0
END



PRO zspec_set_smooth, event
;  Smooth the spectrum and template

  COMMON zspec_set_smooth_COMMON, selection
  COMMON oned_data_COMMON, specdata, smoothspec1, smoothspec2, synspec
  COMMON results_COMMON, results, slitcount, q_index



  widget_control, event.top, get_uvalue=state, /no_copy
  widget_control, event.ID, get_uvalue=uval


  ; These are the smoothings to choose from..
  ;-----------------------------------------
  smooth = [3,10,15,25,30,50,100,150]
  if(uval eq 'ONE') then $
    current = where(smooth EQ state.smooth1)
  if(uval eq 'TWO') then $
    current = where(smooth EQ state.smooth2)


  selection = current  ; By default we leave the smoothing unchanged.



  ; Create a pop-up widget to ask the user for the smoothing.
  ; Execution will stop until a smoothing is selected.  Once
  ; selected the smoothing level is then stored in the
  ; state variable.
  ;--------------------------------------------------------
  base = widget_base(/column, title='Select smoothing', $
                     group_leader=event.top, /modal)
  values = ['3', $
            '10', $
            '15', $
	    '25', $
	    '30', $
	    '50', $
	    '100', $
	    '150']
  smooth_list = CW_BGROUP(base, values, /row, /frame, $
                          label_left='Smoothing', $
                          event_func='zspec_set_smooth_event', $
                          /return_index, /exclusive, set_value=current[0])
  button = widget_button(base, value='Done', uvalue=-1, $
                         event_func='zspec_set_smooth_event', xsize=200)
  widget_control, base, /realize
  xmanager, 'zspec_set_smooth', base



  ; Once the widget is destroyed we save the new smoothing
  ;-------------------------------------------------------
  if(uval eq 'ONE') then $
    state.smooth1 = smooth[selection]
  if(uval eq 'TWO') then $
    state.smooth2 = smooth[selection]
  smoothspec1 = ivarsmooth(specdata.spec,specdata.ivar,state.smooth1)
  smoothspec2 = ivarsmooth(specdata.spec,specdata.ivar,state.smooth2)


  IF(results[slitcount].nonlocal_select EQ 1) THEN BEGIN
    redshift = results[slitcount].zresult_non[results[slitcount].selection].z
  ENDIF ELSE BEGIN
    redshift = results[slitcount].zresult_all[results[slitcount].selection].z
  ENDELSE

  zspec_1d_plot, state, redshift
  IF(xregistered('splot')) then zspec_splot, state, redshift

  widget_control, event.top, set_uvalue=state, /no_copy

END


;-------------------------------------------------------------------------


PRO zspec_smooth, event
; Smooth the 1d spectrum and re-display it

  COMMON results_COMMON, results, slitcount, q_index

  widget_control, event.TOP, get_uvalue=state, /no_copy
  widget_control, event.ID, get_uvalue=uval

  IF(uval eq 'OFF') THEN $
    state.smooth_on = 0
  IF(uval eq 'ONE') THEN $
    state.smooth_on = 1
  IF(uval eq 'TWO') THEN $
    state.smooth_on = 2

  IF(results[slitcount].nonlocal_select EQ 1) THEN BEGIN
    redshift = results[slitcount].zresult_non[results[slitcount].selection].z
  ENDIF ELSE BEGIN
    redshift = results[slitcount].zresult_all[results[slitcount].selection].z
  ENDELSE

  zspec_1d_plot, state, redshift
  if(xregistered('splot')) THEN zspec_splot, state, redshift

  widget_control, event.TOP, set_uvalue=state, /no_copy


END





;-------------------------------------------------------------------------



PRO zspec_showtext, state

  COMMON results_COMMON, results, slitcount, q_index
  COMMON mobject_COMMON, mobject

  widID = state.maintext
  IF(results[slitcount].nonlocal_select EQ 1) THEN $
    zres = results[slitcount].zresult_non[results[slitcount].selection] $
  ELSE $
    zres = results[slitcount].zresult_all[results[slitcount].selection]

  nobj_in_slit = string(mobject[slitcount], format='(I1.1)')
  topdir = getenv('D2_RESULTS')

  bintabs =findfile(topdir +'/' + results[slitcount].maskname + '/*/' + results[slitcount].maskname + '.bintabs.fits')	;Added BCL, 5/21/10 to report mag in zspec, used to find the bintabs file
  bintabsread = mrdfits(bintabs[0], 1, /silent)	;read in the bintabs file
 
  bintabscount = where(strcompress(bintabsread.object,/remove_all) eq results[slitcount].objname)
  if bintabscount[0] ge 0 then imag = strcompress(string(bintabsread[bintabscount].mag), /remove_all) else imag = 'N/A'

  string1 = 'Objno: ' + results[slitcount].objname + ',  Slit no: ' $
            + results[slitcount].slitname +', Nobj: ' + nobj_in_slit
  string2 = '------------------------------------'
  IF(results[slitcount].nonlocal_select EQ 1) THEN $
    string3 = 'Current settings:   (non local sky)' $
  ELSE $
    string3 = 'Current settings:'
  IF(results[slitcount].zmanual EQ 1 AND $
     results[slitcount].selection EQ 9) THEN $
    string4 ='    zmanual = '+string(format='(F9.5)', zres.z) $
  ELSE string4 = '    zbest = ' + string(format='(F9.5)', zres.z)
  string5 = string4 + ',   Q   = ' + string(format='(I3)', results[slitcount].zquality)
  string6 = '  comment = ' + results[slitcount].comment
  string7 = 'Other Info:'
  string8 = '   delta chi^2 = ' + string(zres.rchi2diff * zres.dof)
  string9 = '   eigenvalues = ' + string(format = '(2f8.2)', zres.theta[0:1])
  string10 = '   A star = ' + string(format = '(2f8.2)', zres.theta[2]) + ',   i mag = ' + imag 


  widget_control, widID, set_value=string1
  widget_control, widID, set_value=string2, /append
  widget_control, widID, set_value=string3, /append
  widget_control, widID, set_value=string5, /append
  widget_control, widID, set_value=string6, /append
  widget_control, widID, set_value=string2, /append
  widget_control, widID, set_value=string7, /append
  widget_control, widID, set_value=string8, /append
  widget_control, widID, set_value=string9, /append
  widget_control, widID, set_value=string10, /append
  ;widget_control, widID, set_value=string11, /append


END

;-----------------------------------------------------------------------
FUNCTION zspec_check_qualities
	COMMON results_COMMON, results, slitcount, q_index

	roi = where(results.zquality eq 0, count)

	stop_quit = 1

	if count gt 0 then begin
		n = string(count, format='(i3)')
		res = dialog_message(["Warning: You have " + n + $
			" objects with zquality = 0","Would you like to exit?"],$
			/Question)
		if res eq 'Yes' then stop_quit = 0
	endif

return, stop_quit
end


;-----------------------------------------------------------------------
PRO zspec_save_event, event

  COMMON zspec_save_COMMON, savename, cancel, dumpname

  widget_control, event.ID, get_uvalue=uval
  IF(uval EQ -1) THEN BEGIN
    cancel = 1
    widget_control, event.TOP, /destroy
  ENDIF ELSE BEGIN
    widget_control, uval, get_value=savename
    widget_control, event.TOP, /destroy
  ENDELSE

END



PRO zspec_save, event
; If the save file name has not already been specified then prompt
; the user for a file name.  Otherwise just save the results.

  COMMON results_COMMON, results, slitcount, q_index
  COMMON zspec_save_COMMON, savename, cancel, dumpname
  COMMON mask_comment_COMMON, maskcomment

  cancel = 0

  widget_control, event.TOP, get_uvalue=state
  base = widget_base(/column, title='Enter file name', $
                       group_leader=event.top, /modal)
  label = widget_label(base,value='Enter the name of file to be saved:')
  label = widget_label(base,value='     (use a .fits suffix)')
  text  = widget_text(base, xsize=50, uvalue=0L, $
                        event_PRO='zspec_save_event', $
                        value=savename, /editable)
  widget_control, text, set_uvalue=text
  button = widget_button(base,value='Save', uvalue=text, $
                           event_PRO='zspec_save_event')
  button = widget_button(base,value='Dont save', uvalue='-1', $
                           event_PRO='zspec_save_event')
  widget_control, base, /realize
  xmanager, 'zspec_save', base

  IF(cancel EQ 1) THEN BEGIN
    return
  ENDIF

   t = findfile(dumpname[0])
   if t ne '' then begin
   	res = dialog_message('Are you sure you want to overwrite?',/Question)
   	if res eq 'No' then return
   endif


  dumpname = $                    ; Dump file for restarting zspec
            strmid(savename, 0, strlen(savename)-5) + '.sav'
  zspec_save_file, event
END


PRO zspec_save_file, event;, savename, dumpname
  COMMON results_COMMON, results, slitcount, q_index
  COMMON zspec_save_COMMON, savename, cancel, dumpname
  COMMON mask_comment_COMMON, maskcomment

; check for all cases where z_err = 999.9. for these objects, make
; sure that zquality = 2. also for all objects with z_err = 0.0, set
; the zquality to zero.
  zman = where(results.zresult.z_err eq 999.9, zmanct)
  if zmanct gt 0 then begin
      tmp = results[zman].zresult
      tmp.zquality = 2
      results[zman].zresult = tmp
      results[zman].zquality = 2
  endif
  zsc = where(results.zresult.z_err eq 0.0, zsct)
  if zsct gt 0 then begin
      tmp = results[zsc].zresult
      tmp.zquality = 0
      results[zsc].zresult = tmp
      results[zsc].zquality = 0
  endif

; now save the output...
  save, results, filename=dumpname[0]
  hdr = strarr(n_elements(results)+1)
  hdr[0] = 'Mask:   ' + results[0].maskname
  FOR i=0, n_elements(results)-1 DO BEGIN
    hdr[i+1] = string(format='(A5,A12,3I4,F13.5,I4,A3,A0)', $
                    results[i].slitname, $
                    results[i].objname, $
                    results[i].selection, $
                    results[i].zmanual, $
                    results[i].nonlocal_select, $
                    results[i].zresult.z, $
                    results[i].zquality, $
                    '  ', $
                    results[i].comment $
                    )
  ENDFOR
  results.zresult.comment = results.comment
  mwrfits, results.zresult, savename, hdr, /create
  mwrfits, maskcomment, savename
END



;-----------------------------------------------------------------------

PRO zspec_comment_event, event
  COMMON zspec_comment_COMMON, comment

  widget_control, event.ID, get_uvalue=uval
  widget_control, uval, get_value=comment
  widget_control, event.TOP, /destroy

END



PRO zspec_comment, event
; Prompt the user for a comment about the object

  COMMON results_COMMON, results, slitcount, q_index
  COMMON zspec_comment_COMMON, comment

  widget_control, event.TOP, get_uvalue=state

  ; Create a child widget
  ;-----------------------
  base = widget_base(/column, title='Enter comment', $
                     group_leader=event.top, /modal)
  label = widget_label(base,value='Enter your comment below:')
  text  = widget_text(base, xsize=50, uvalue=0L, $
                      event_PRO='zspec_comment_event', $
		      value = results[slitcount].comment, $
                      /editable)
  widget_control, text, set_uvalue=text
  button = widget_button(base,value='Done', uvalue=text, $
                       event_PRO='zspec_comment_event')

  widget_control, base, /realize
  xmanager, 'zspec_comment', base


  results[slitcount].comment = comment
  zspec_showtext, state

END


PRO zspec_comment_mask_event, event
	COMMON mask_comment_COMMON, maskcomment

	widget_control, event.ID, get_uvalue=uval
	widget_control, uval, get_value=maskcomment
	widget_control, event.TOP, /destroy
END

PRO zspec_comment_mask, event
	COMMON mask_comment_COMMON, maskcomment

	widget_control, event.TOP, get_uvalue=state

	base = widget_base(/column, title='Enter Mask Comment:',$
		group_leader=event.TOP, /modal)
	label = widget_label(base, value='Enter a Mask Comment')
	text = widget_text(base,xsize=50,ysize=25,uvalue=0L, $
		event_PRO='zspec_comment_mask_event', $
		value=maskcomment, /editable)
	widget_control, text, set_uvalue=text
	button = widget_button(base, value='Done', uvalue=text, $
		event_PRO='zspec_comment_mask_event')
	widget_control, base, /realize
END



PRO zspec_atv_overplot_lines, lambdabottom, lambdatop, linewaves, linenames, nbins, height, color
	; iteratively go through the above values
	; and plop them down on the ATV

	right = 1
  	FOR i=0,n_elements(linenames)-1 DO BEGIN
		l = linewaves[i]

      		; convert wavelength into an index:
      		nl0 = where(lambdabottom ge l-.5 and lambdabottom le l+.5, count)
      		if count lt 1 then continue
      		nl0 = nl0[(n_elements(nl0)-1)/2]
      		nlbottom = nl0[0]

      		nl0 = where(lambdatop ge l-.5 and lambdatop le l+.5, count)
      		if count lt 1 then continue
      		nl0 = nl0[(n_elements(nl0)-1)/2]
      		nltop = nl0[0]

      		; figure out the x & y positions
      		x = fix((float(nlbottom)/nbins-nlbottom/nbins)*nbins)
      		y = height*(7-nlbottom/nbins)-1

		if right then $
      			atvxyouts, x+10, y, linenames[i], color=color, charsize=2 $
		else $
      			atvxyouts, x-10-10*strlen(linenames[i]), y, linenames[i], color=color, charsize=2
		right = not right

      		atvplot, fltarr(7)+x, findgen(7)+y, color=color, thick=2
      		x = fix((float(nltop)/nbins-nl0/nbins)*nbins)
      		y = height*(7-nltop/nbins)
      		atvplot, fltarr(7)+x, height-4-findgen(7)+y, color=color,thick=2

      		if linenames[i] eq 'OII' then i = i+2
		if linenames[i] eq 'CA' then i = i+2
  	ENDFOR
END


;-----------------------------------------------------------------


PRO zspec_update_atv, state

  COMMON twod_data_COMMON, fullfluxdata, fluxdata, extraction_r1, extraction_r2, extraction_offset, slit_display_options
  COMMON zspec_lines_COMMON, linenames, linewaves, telluricnames, telluricwaves
  COMMON results_COMMON, results, slitcount, q_index

  redshift = zspec_get_redshift()


  if slit_display_options.fullslit then begin
	zflux = fullfluxdata
  endif else begin
  	zflux = fluxdata
  endelse

  nbins = (size(zflux.blambda, /dimensions))[0]/4

  
  height = n_elements(zflux.bflux[0,*])+4 ; 4 for the four rows of lam info
  blambda0 = zflux.blambda[*,0]
  rlambda0 = zflux.rlambda[*,0]

  ; Split spectrum into 8 segments and then recombine
  flux1 = zflux.bflux[0:nbins-1,*]
  flux2 = zflux.bflux[nbins:2*nbins-1,*]
  flux3 = zflux.bflux[2*nbins:3*nbins-1,*]
  flux4 = zflux.bflux[3*nbins:4*nbins-1,*]
  lam1  = blambda0[0:nbins-1]
  lam2  = blambda0[nbins:2*nbins-1]
  lam3  = blambda0[2*nbins:3*nbins-1]
  lam4  = blambda0[3*nbins:4*nbins-1]
  flux_blue = [[flux4],[lam4],[lam4],[lam4],[lam4],[flux3],[lam3],[lam3],[lam3],[lam3], $
  	[flux2],[lam2],[lam2],[lam2],[lam2],[flux1],[lam1],[lam1],[lam1],[lam1]]

  flux1 = zflux.rflux[0:nbins-1,*]
  flux2 = zflux.rflux[nbins:2*nbins-1,*]
  flux3 = zflux.rflux[2*nbins:3*nbins-1,*]
  flux4 = zflux.rflux[3*nbins:4*nbins-1,*]
  lam1  = rlambda0[0:nbins-1]
  lam2  = rlambda0[nbins:2*nbins-1]
  lam3  = rlambda0[2*nbins:3*nbins-1]
  lam4  = rlambda0[3*nbins:4*nbins-1]
  flux_red = [[flux4],[lam4],[lam4],[lam4],[lam4],[flux3],[lam3],[lam3],[lam3],[lam3], $
  	[flux2],[lam2],[lam2],[lam2],[lam2],[flux1],[lam1],[lam1],[lam1],[lam1]]

  atv, [[flux_red],[flux_blue]], min=-50, max=100, /stretch ;_extra=extra, /stretch


  hlambda = (size(zflux.blambda))[2] < (size(zflux.rlambda))[2]
  lambdabottom = [blambda0, rlambda0]
  lambdatop = [zflux.blambda[*,hlambda-1], zflux.rlambda[*,hlambda-1]]
  offset = (height - (size(zflux.blambda, /dimensions))[1])/2>0
  zspec_atv_overplot_lines, lambdabottom, lambdatop, linewaves*(1.+redshift), $
  	linenames, nbins, height, 6
  zspec_atv_overplot_lines, lambdabottom, lambdatop, telluricwaves, telluricnames, $
  	nbins, height, 5

      ;

      ; Draw Ticks to indicate Extraction Window

  offset = abs(extraction_offset)


  for i=0, 7 do begin
      base = fltarr(10) + height*i + offset - 1
  	; left side
      atvplot, findgen(10), base + extraction_r1, color=1, thick=2
      atvplot, findgen(10), base + extraction_r2, color=1, thick=2
      	; right side
      atvplot, 1024-findgen(10), base + extraction_r1, color=1, thick=2
      atvplot, 1024-findgen(10), base + extraction_r2, color=1, thick=2
  endfor


  IF(arg_present(state)) THEN $
    widget_control, state.base, /show

END




PRO zspec_atv, event
; Show full 2d spectrum in an ATV window

  zspec_update_atv

END


;-------------------------------------------------------------------------

FUNCTION zspec_setz, event, state=state

  COMMON results_COMMON, results, slitcount, q_index

; check the z_by_hand keyword.
  if n_elements(state) gt 0 then z_by_hand = 1 else z_by_hand = 0

  if z_by_hand eq 0 then $
    widget_control, event.TOP, get_uvalue=state

  qual = [-2,-1,1,2,3,4]

  results[slitcount].qindex = event.value
  results[slitcount].zquality = qual[event.value]

  IF(results[slitcount].nonlocal_select EQ 1) THEN BEGIN
    results[slitcount].zresult = $
	results[slitcount].zresult_non[results[slitcount].selection]
  ENDIF ELSE BEGIN
    results[slitcount].zresult = $
	results[slitcount].zresult_all[results[slitcount].selection]
  ENDELSE


  results[slitcount].zresult.zquality = qual[event.value]
  results[slitcount].zresult.comment  = results[slitcount].comment
  zspec_showtext, state
  IF(state.autonext EQ 1) THEN $
    zspec_GOTO_next, state

END


;---------------------------------------------------------------------

PRO zspec_select_nonlocal, event

  COMMON results_COMMON, results, slitcount, q_index
  COMMON oned_data_COMMON, specdata, smoothspec1, smoothspec2, synspec
  COMMON twod_data_COMMON, fullfluxdata, fluxdata, extraction_r1, extraction_r2, extraction_offset, slit_display_options

  widget_control, event.TOP, get_uvalue=state



  IF(results[slitcount].nonlocal_select EQ 1) THEN BEGIN
    widget_control, state.nonlocal, set_value='Show non-local sky'
    results[slitcount].nonlocal_select = 0

    widget_control, /hourglass
    widget_control, state.redshift_select, $
                  set_value = string(format='(F9.5)', $
                                     results[slitcount].zresult_all.z)
;    widget_control, state.redshift_select, $
;                  set_droplist_select=results[slitcount].selection
    results[slitcount].selection = 0
    zspec_showtext, state
    redshift = results[slitcount].zresult_all[results[slitcount].selection].z

  ENDIF ELSE BEGIN
    widget_control, state.nonlocal, set_value='Show local sky'
    results[slitcount].nonlocal_select = 1

    widget_control, /hourglass
    widget_control, state.redshift_select, $
                  set_value = string(format='(F9.5)', $
                                     results[slitcount].zresult_non.z)
;    widget_control, state.redshift_select, $
;                  set_droplist_select=results[slitcount].selection
    results[slitcount].selection = 0
    zspec_showtext, state
    redshift = results[slitcount].zresult_non[results[slitcount].selection].z
  ENDELSE


    specdata = zspec_1d_read(results[slitcount].nonlocal_select)
    smoothspec1 = ivarsmooth(specdata.spec,specdata.ivar,state.smooth1)
    smoothspec2 = ivarsmooth(specdata.spec,specdata.ivar,state.smooth2)
    synspec  = zspec_template(specdata.lambda,state.eigenvec)
    zspec_1d_plot, state, redshift
    IF(xregistered('splot')) THEN  zspec_splot, state, redshift
    fluxdata = zspec_2d_read(specdata,results[slitcount].nonlocal_select)
    if slit_display_options.fullslit then begin
            fullfluxdata = zspec_2d_read(specdata,results[slitcount].nonlocal_select, /full)
    endif

    zspec_2d_plot, state, redshift
    if(xregistered('atv')) then zspec_update_atv, state

end


;-----------------------------------------------------------------------
FUNCTION zspec_set_fullspec_event, event
  COMMON zspec_manualz_COMMON, wavelength, feature, agree, lambdarange, fullspec, tmpdex

  widget_control, event.ID, get_uvalue=uval
  IF(uval EQ -1) THEN  $
    widget_control, event.top, /destroy $
  ELSE $
    fullspec = event.value
  return, fullspec
END

FUNCTION zspec_set_tmp_choice, event
  COMMON zspec_manualz_COMMON, wavelength, feature, agree, lambdarange, fullspec, tmpdex

  widget_control, event.ID, get_uvalue=uval
  IF(uval EQ -1) THEN  $
    widget_control, event.top, /destroy $
  ELSE $
    tmpdex = event.value
  return, tmpdex
END

PRO zspec_manualz_event, event

  COMMON zspec_manualz_COMMON, wavelength, feature, agree, lambdarange, fullspec, tmpdex

  widget_control, event.ID, get_uvalue=uval
  widget_control, event.TOP, get_uvalue=zstate

  IF(uval EQ 'DONE' OR uval EQ 'TEXT' OR uval EQ 'TEXT_WIDTH') THEN BEGIN
    widget_control, zstate.text, get_value=wave
    widget_control, zstate.textw, get_value=range
    wavelength = float(wave[0])
    lambdarange = float(range[0])
    widget_control, event.TOP, /destroy
  ENDIF ELSE IF(uval EQ 'LIST') THEN BEGIN
    feature = event.index
  ENDIF

END


PRO zspec_manualz_event2, event

  COMMON zspec_manualz_COMMON, wavelength, feature, agree, lambdarange, fullspec, tmpdex

  widget_control, event.ID, get_uvalue=uval
  IF(uval EQ 'YES') THEN BEGIN
    agree = 1
    widget_control, event.TOP, /destroy
  ENDIF ELSE BEGIN
    agree = 0
    widget_control, event.top, /destroy
  ENDELSE

END

PRO zspec_z_by_hand_event, event
	COMMON zspec_z_by_hand_COMMON, zres

	widget_control, event.ID, get_uvalue=uval
  	widget_control, event.TOP, get_uvalue=zstate

	IF(uval EQ 'DONE' OR uval EQ 'TEXT') THEN BEGIN
		widget_control, zstate.text, get_value=zmes
		zres = float(zmes[0])
		widget_control, event.TOP, /destroy
	ENDIF
END

PRO zspec_z_by_hand, event
        COMMON zspec_z_by_hand_COMMON, zres
  	COMMON results_COMMON, results, slitcount, q_index
  	COMMON oned_data_COMMON, specdata, smoothspec1, smoothspec2, synspec

	widget_control, event.TOP, get_uvalue=state
	base = widget_base(/column, title='Enter a Redshift', $
		group_leader = event.TOP, /modal)
	label = widget_label(base, value = 'Enter a Redshift')
	text = widget_text(base, xsize = 50, uvalue = 'TEXT', $
		event_PRO='zspec_z_by_hand_event', /editable)

  	button = widget_button(base,value='Done', uvalue='DONE')
  	zstate = {text:text}
  	widget_control, base, set_uvalue=zstate

  	widget_control, base, /realize
  	xmanager, 'zspec_z_by_hand', base

  	IF(NOT finite(zres)) THEN BEGIN
    		res = dialog_message('You did not enter a good z value', /error)
    		return
  	ENDIF
  	zres_new = zspec_zfind(zres, results[slitcount].zresult, specdata, $
                         zmin=zres, zmax=zres)

; force the error value for the "by-hand" redshifts to be very large.
        zres_new.z_err = 999.9
        zres_new.vdisp = 999.9
        zres_new.vdisp_err = 999.9
        zres_new.zquality = 2

        if(results[slitcount].nonlocal_select EQ 1) then $
    		results[slitcount].zresult_non[9] = zres_new $
  	else  $
    		results[slitcount].zresult_all[9] = zres_new
  	results[slitcount].selection = 9
  	results[slitcount].zmanual = 1

  	IF(state.qonly EQ 1) THEN BEGIN
    		pos = where(q_index EQ slitcount, ct)
    		IF(pos GT 0) THEN slitcount = q_index[pos[0]-1] $
    		ELSE slitcount = -1
  	ENDIF ELSE IF(slitcount GT 0) THEN $
    		slitcount = slitcount-1 $
  	ELSE slitcount = -1

; tell the user that this object should be set as a Q=2!
        mess_txt = string(format='(%"%s\n%s\n%s\n%s")', $
                          'You have just entered a redshift manually, ' + $
                          'without performing a detailed fit to the', $
                          'profile. For this reason the Q has been ' + $
                          'automatically set to Q=2. Do not', $
                          'change this to Q=3 or Q=4!')
        res = dialog_message(mess_txt)


; force the zquality to be 2!
        slitcount = slitcount + 1
        foo = {value:3}
        foo = zspec_setz(foo, state=state)
        slitcount = slitcount - 1
;        results[slitcount].zquality = 2
;        results[slitcount].zresult.zquality = 2

        zspec_GOTO_next, state

END

PRO zspec_manualz, event

  COMMON zspec_manualz_COMMON, wavelength, feature, agree, lambdarange, fullspec, tmpdex
  COMMON results_COMMON, results, slitcount, q_index
  COMMON oned_data_COMMON, specdata, smoothspec1, smoothspec2, synspec
  COMMON zspec_lines_COMMON, linenames, linewaves, telluricnames, telluricwaves

  lambdarange = 20.0

  widget_control, event.TOP, get_uvalue=state

  base = widget_base(/column, title='Enter a redshift', $
                     group_leader=event.top, /modal)
  label = widget_label(base,value='Enter the (air) wavelength of the feature')
  label = widget_label(base,value='below and then select it from the list')
  text  = widget_text(base, xsize=50, uvalue='TEXT', $
                      event_PRO='zspec_manualz_event', $
                      /editable)
  label = widget_label(base,value='Enter a wavelength range in Angstroms')
  textw = widget_text(base, xsize=50, uvalue='TEXT_WIDTH', $
  			value = strtrim(string(lambdarange), 1), $
                      event_PRO='zspec_manualz_event', $
                      /editable)

  values = strarr(n_elements(linewaves))
  for i = 0, n_elements(linewaves)-1 do begin
  	values[i] = string(linenames[i],format='(a10)') + string(linewaves[i])
  endfor
; note that the LyA, CIV lines are not in the template wavelength
; range. and note that Hdelta, Hgamma, NeIII, etc. are not in the
; template.
  values[0] = string('(bad) Lya', format='(a10)') + '      N/A'
  values[1] = string('(bad) CIV', format='(a10)') + '      N/A'
  values[9] = string('(*) NeIII', format='(a10)') + string(linewaves[9])
  values[10] = string('(*) HeI', format='(a10)') + string(linewaves[10])
  values[13] = string('(*) He', format='(a10)') + string(linewaves[13])
  values[14] = string('(*) Hd', format='(a10)') + string(linewaves[14])
  values[16] = string('(*) Hg', format='(a10)') + string(linewaves[16])

  select = widget_droplist(base,  value=values, $
                           title='Line selection (air):', $
                           event_pro='zspec_manualz_event', $
                           /dynamic_resize, /frame, uvalue='LIST')
  answers = ['No', 'Yes']
  fullspec = 1
  fullspecval = CW_BGROUP(base, answers, /row, /frame, $
                       label_left='Fit over entire spectrum?', $
                       event_func='zspec_set_fullspec_event', $
                       /return_index, /exclusive, set_value=1)

  tmps = ['Galaxy', 'Abs-Line', 'Emi-Line', 'A-star', 'Star', 'QSO']
  tmpdex = 0
  tmpval = CW_BGROUP(base, tmps, /row, /frame, $
                     label_left='Template(s):', $
                     event_func='zspec_set_tmp_choice', $
                     /return_index, /exclusive, set_value=0)

  button = widget_button(base,value='Done', uvalue='DONE')
  zstate = {text:text, $
  	    textw:textw, $
            select:select $
           }
  feature = 0  ; Default to [OII]
  widget_control, base, set_uvalue=zstate

  widget_control, base, /realize
  xmanager, 'zspec_manualz', base

  lines = linewaves
  redshift = (wavelength/lines[feature]) - 1.0
  print, 'First guess z = ', redshift

  if lambdarange eq 0.0 then lambdarange = 100.0
  zrange =  lambdarange / lines[feature]
  zmin = redshift[0] - (zrange/2.0)
  zmax = redshift[0] + (zrange/2.0)

  lam_min = lines[feature] * (1.0 + zmin)
  lam_max = lines[feature] * (1.0 + zmax)
  wvbuf =  10.0
  cntl = where(specdata.lambda GE (lam_min - wvbuf) AND $
              specdata.lambda LE (lam_max + wvbuf), ct)
  IF(ct LT 1) THEN BEGIN
    res=dialog_message('Error: You selected an invalid wavelength.', /error)
    return
  ENDIF



  base = widget_base(/column, title='Is this OK?', $
                     group_leader=event.top, /modal)
  label = widget_label(base,value='I am going to scan between the dotted')
  label = widget_label(base,value='lines.  Is this OK?')
  draw  = widget_draw(base, xsize=250, ysize=200)
  button = widget_button(base, value='Yes', uvalue='YES')
  button = widget_button(base, value='No', uvalue='NO')

  widget_control, base, /realize
  widget_control, draw, get_value=draw_ID

  wset, draw_ID

  plot, specdata.lambda[cntl], specdata.spec[cntl], xticklayout=0, $
    xtickinterval=lambdarange/3.0
  miny = min(specdata.spec)-1000.0
  maxy = max(specdata.spec)+1000.0
;  lam_min = wavelength[0]*(1.+redshift[0]-0.01)/(1.+redshift[0])
;  lam_max = wavelength[0]*(1.+redshift[0]+0.01)/(1.+redshift[0])
  oplot, [lam_min, lam_min], [miny, maxy], line=2, color=1, thick=2
  oplot, [lam_max, lam_max], [miny, maxy], line=2, color=1, thick=2

  agree = 0
  xmanager, 'zspec_manualz', base, event_handler='zspec_manualz_event2'


  IF(agree EQ 0) THEN return

  PRINT, 'Now fitting...'
  widget_control, /hourglass

  if tmpdex eq 4 then begin
      eigenfile = 'spEigenStarDeep2.fits'
      star = 1
  endif else star = 0
  if tmpdex eq 5 then begin
      eigenfile = 'spEigenQSOdeep.fits'
      qso = 1
  endif else qso = 0
  if tmpdex eq 1 then begin
      eigenfile = 'spDEEP.fits'
      columns = 0
  endif
  if tmpdex eq 2 then begin
      eigenfile = 'spDEEP.fits'
      columns = 1
  endif
  if tmpdex eq 3 then begin
      eigenfile = 'spDEEP.fits'
      columns = 2
  endif
  if tmpdex eq 0 then begin
      eigenfile = 'spDEEP.fits'
      columns = [0,1,2]
  endif

  if fullspec then $
    zres_new = zspec_zfind(redshift[0], results[slitcount].zresult, $
                           specdata, zmin=zmin, zmax=zmax, $
                           columns=columns, eigenfile=eigenfile, $
                           qso=qso, star=star) $
  else $
    zres_new = zspec_zfind(redshift[0], results[slitcount].zresult, $
                           specdata, zmin=zmin, zmax=zmax, $
                           wmin=lam_min, wmax=lam_max, $
                           columns=columns, eigenfile=eigenfile, $
                           qso=qso, star=star)
;print, 'COL: ', columns
;print, 'THETA: ', zres_new.theta
  print, 'New redshift = ', zres_new.z

  IF(NOT finite(zres_new.theta[0])) THEN BEGIN
    res = dialog_message('Fit failed, please try again perhaps with a different wavelength selection.', /error)
    return
  ENDIF


  ; Create a new template spectrum of the fit
  ;--------------------------------------------------
  eigen = state.eigenvec
  IF(zres_new.tfile EQ 'spEigenStarDeep2.fits') THEN BEGIN
    eig = eigen.star
    h   = eigen.h1
  ENDIF ELSE BEGIN
    eig = eigen.galaxy
    h   = eigen.h2
  ENDELSE

  w = where(zres_new.tcolumn GE 0, nw)
  syn = 0

  ; Temporary fix to deal with stars
  IF nw LT 1 THEN BEGIN
    nw = 1
    w = [0]
    zres_new.tcolumn[w[0]] = 10
  ENDIF

  for i=0,nw-1 do begin
     jj = zres_new.tcolumn[w[i]]
     syn = syn+zres_new.theta[i]*eig[*,jj]
  endfor

  IF zres_new.npoly gt 0 THEN begin ;add in the polynomial terms
      npoints = n_elements(eig[*,0])
      parray = poly_array(npoints, zres_new.npoly)
      for i=0,zres_new.npoly-1 do $
        syn = syn + zres_new.theta[nw+i]*parray[*,i]
  ENDIF
  if tmpdex le 2 then begin
      cont = median(specdata.spec)
      npix = n_elements(specdata.lambda)
      restlambda = specdata.lambda / (1.0 + zres_new.z)
      minpix = findpix(eigen.loglam, alog10(restlambda[0]))
      maxpix = findpix(eigen.loglam, alog10(restlambda[npix-1]))
      tmp_med = median(eigen.tmpsmth[minpix:maxpix])
      syn = syn + eigen.tmpsmth * cont / tmp_med
  endif

  coeff0 = sxpar(h, 'COEFF0')
  coeff1 = sxpar(h, 'COEFF1')

  elambda = 10.d^(coeff0+dindgen((size(eig, /dimens))[0])*coeff1)
; ??? interpol is pretty stupid - use combine1fiber in idlspec2d (SDSS)
  newsynspec = interpol(syn, elambda, specdata.lambda/(zres_new.z+1))


  base = widget_base(/column, title='Is this OK?', $
                     group_leader=event.top, /modal)
  label = widget_label(base,value='This is the best fit template')
  label = widget_label(base,value='z = ' + $
                       strcompress(zres_new.z[0], /rem) )
  label = widget_label(base,value='        Is this OK?          ')
  draw  = widget_draw(base, xsize=250, ysize=200)
  text  = widget_text(base,ysize=7,xsize=40)

  button = widget_button(base, value='Yes', uvalue='YES')
  button = widget_button(base, value='No', uvalue='NO')

  widget_control, base, /realize
  widget_control, draw, get_value=draw_ID
  widget_control, text, set_value='If you select YES then this fit'
  widget_control, text, set_value='will replace the 10th auto fit.', /append
  widget_control, text, set_value='You can then switch between it', /append
  widget_control, text, set_value='and the other fits using the drop', /append
  widget_control, text, set_value='list as before, but you will not', /append
  widget_control, text, set_value='be able to recover the previous', /append
  widget_control, text, set_value='10th entry.', /append

  wset, draw_ID

  plot, specdata.lambda[cntl], specdata.spec[cntl], background=7, $
    color=0, xticklayout=0, xtickinterval=lambdarange/3.0
  oplot, specdata.lambda[cntl], newsynspec[cntl], color=1, thick=3

  agree = 0
  xmanager, 'zspec_manualz', base, event_handler='zspec_manualz_event2'

  IF(agree EQ 1) THEN begin
      print, 'You agreed!'
;      if tmpdex eq 1 or tmpdex eq 2 then begin
;          mess_txt = string(format='(%"%s\n%s\n%s\n%s")', $
;                            'You have just fit for a redshift using only ' + $
;                            '1 template. For this reason, please include', $
;                            'the "iffy" comment for this object!')
;          res = dialog_message(mess_txt)
;      endif
; if fit was done with only 1 template (emission-line or
; absorption-line only), then set the other eigenvalue to a flagvalue:
; -1.0.
      if tmpdex eq 1 then zres_new.theta[1:2] = -1.0
      if tmpdex eq 2 then begin
          zres_new.theta[0] = -1.0
          zres_new.theta[2] = -1.0
      endif
      if tmpdex eq 3 then zres_new.theta[0:1] = -1.0
  endif ELSE BEGIN
    print, 'you disagreed...'
    return
  ENDELSE

  if(results[slitcount].nonlocal_select EQ 1) then $
    results[slitcount].zresult_non[9] = zres_new $
  else $
    results[slitcount].zresult_all[9] = zres_new
  results[slitcount].selection = 9
  results[slitcount].zmanual = 1

  ; We now pretend to go back to the previous spectrum and then start the
  ; goto_next routine.  This way the screen will be refreshed with the
  ; new entry.
  IF(state.qonly EQ 1) THEN BEGIN
    pos = where(q_index EQ slitcount, ct)
    IF(pos GT 0) THEN slitcount = q_index[pos[0]-1] $
    ELSE slitcount = -1
  ENDIF ELSE IF(slitcount GT 0) THEN $
    slitcount = slitcount-1 $
  ELSE slitcount = -1

  zspec_GOTO_next, state

END









;------------------------------------------------------------------------

PRO zspec_jump_event, event

  COMMON zspec_jump_COMMON, destination

  widget_control, event.ID, get_uvalue=uval
  widget_control, uval, get_value=destination
  widget_control, event.TOP, /destroy

END

 ;     Untilt the slit as requested
PRO zspec_untiltslit, event
 	COMMON oned_data_COMMON, specdata, smoothspec1, smoothspec2, synspec
 	COMMON twod_data_COMMON, fullfluxdata, fluxdata, extraction_r1, $
		extraction_r2, extraction_offset, slit_display_options
 	COMMON twod_data_HACK_COMMON, hackfluxdata 
 	COMMON results_COMMON, results, slitcount, q_index
 
	slit_display_options.untilt = not slit_display_options.untilt
	if slit_display_options.fullslit then begin
		fullfluxdata = zspec_2d_read(specdata, $
			results[slitcount].nonlocal_select, /fullspec)
	endif
	fluxdata = zspec_2d_read(specdata, results[slitcount].nonlocal_select)


	widget_control, event.top, get_uvalue=state
	if (xregistered('atv')) then zspec_update_atv, state
END

;     Show the whole slit
PRO zspec_showfullslit, event

	COMMON twod_data_COMMON, fullfluxdata, fluxdata, extraction_r1, $
		extraction_r2, extraction_offset, slit_display_options
	COMMON oned_data_COMMON, specdata, smoothspec1, smoothspec2, synspec
	COMMON results_COMMON, results, slitcount, q_index

; npk: Major hack
	COMMON twod_data_HACK_COMMON, hackfluxdata 
	widget_control, event.top, get_uvalue=state

	slit_display_options.fullslit = not slit_display_options.fullslit

	if slit_display_options.fullslit then begin
		fullfluxdata = zspec_2d_read(specdata, $
			results[slitcount].nonlocal_select, /fullspec)
	endif

	if (xregistered('atv')) then zspec_update_atv, state
END



PRO zspec_jump, event
; Jump to a specified slitnumber in file

  COMMON zspec_jump_COMMON, destination
  COMMON results_COMMON, results, slitcount, q_index
  COMMON oned_data_COMMON, specdata, smoothspec1, smoothspec2, synspec
  COMMON twod_data_COMMON, fullfluxdata, fluxdata, extraction_r1, extraction_r2, extraction_offset, slit_display_options


  widget_control, event.TOP, get_uvalue=state

  IF(state.qonly EQ 1) THEN BEGIN
    mess_txt = string(format='(%"%s\n%s\n%s\n%s\n%s")', $
                      'This routine is not compatible with the your', $
                      'selection of the SHOW ONLY Q=? option.', $
                      'If you need to jump around this file you will', $
                      'need to restart zspec and not select the ', $
                      ' SHOW ONLY Q=? option from the view menu.')
    res = dialog_message(mess_txt, /error, dialog_parent=event.top)
    return
  ENDIF


  base = widget_base(/column, title='Enter slit number', $
                     group_leader=event.top, /modal)
  label = widget_label(base,value='Enter destination slit number below:')
  label = widget_label(base,value='(Note: not all slits are available)')
  text  = widget_text(base, xsize=50, uvalue=0L, $
                      event_PRO='zspec_jump_event', $
                      /editable)
  widget_control, text, set_uvalue=text
  button = widget_button(base,value='Done', uvalue=text, $
                       event_PRO='zspec_jump_event')

  widget_control, base, /realize
  xmanager, 'zspec_jump', base


  ; Check to see this slit is available, i.e. hasn't been removed before
  cnt = where(fix(results.slitname) EQ fix(destination[0]), c)
  IF c LT 1 THEN BEGIN
    mess_txt = string(format='(%"%s\n%s")', $
                      'Sorry, but that slit is not available in this file.', $
                      'Please check that it is not a sky slit or guide star.')
    res = dialog_message(mess_txt)
  ENDIF ELSE BEGIN
  ; Otherwise read in this spectrum and display it.
    slitcount = fix(cnt[0])

    widget_control, /hourglass
    IF(results[slitcount].nonlocal EQ 0) THEN $
      widget_control, state.nonlocal, sensitive=0 $
    ELSE $
      widget_control, state.nonlocal, sensitive=1

    IF(results[slitcount].nonlocal_select EQ 1) THEN BEGIN
      widget_control, state.redshift_select, $
                    set_value = string(format='(F9.5)', $
                                       results[slitcount].zresult_non.z)
      widget_control, state.redshift_select, $
                    set_droplist_select=results[slitcount].selection
      zspec_showtext, state
      redshift = $
              results[slitcount].zresult_non[results[slitcount].selection].z
          widget_control, state.nonlocal, set_value='Show local sky'
    ENDIF ELSE BEGIN
    ; Normal sky subtraction options
    widget_control, state.redshift_select, $
                  set_value = string(format='(F9.5)', $
                                     results[slitcount].zresult_all.z)
    widget_control, state.redshift_select, $
                  set_droplist_select=results[slitcount].selection
    zspec_showtext, state
    redshift = $
            results[slitcount].zresult_all[results[slitcount].selection].z
          widget_control, state.nonlocal, set_value='Show non-local sky'
    ENDELSE
    specdata = zspec_1d_read(results[slitcount].nonlocal_select)
    smoothspec1 = ivarsmooth(specdata.spec,specdata.ivar,state.smooth1)
    smoothspec2 = ivarsmooth(specdata.spec,specdata.ivar,state.smooth2)
    synspec  = zspec_template(specdata.lambda,state.eigenvec)
    IF(xregistered('splot') GT 0) THEN $
       zspec_splot, state, redshift
    zspec_1d_plot, state, redshift
    fluxdata = $
            zspec_2d_read(specdata,results[slitcount].nonlocal_select)
	if slit_display_options.fullslit then begin
		fullfluxdata = zspec_2d_read(specdata,results[slitcount].nonlocal_select, /full)
	endif

    zspec_2d_plot, state, redshift
    IF(xregistered('atv') GT 0) THEN $
       zspec_update_atv, state
ENDELSE


END



;--------------------------------------------------------------------------


PRO zspec_version_event, event
  widget_control, event.top, /destroy
END


PRO zspec_version, event

  base  = widget_base(/column,title='ZSPEC',group_leader=event.top)
  title = widget_label(base, value='ZSPEC Beta (v1_0)', font='lucidasans-24')
  text1 = widget_label(base, value='Programmed by: D.S Madgwick')
  text2 = widget_label(base, value='  dsm@astron.berkeley.edu  ')

  done  = widget_button(base, value='Close', xsize=60)

  widget_control, base, /realize
  xmanager, 'zspec_version', base

END





;----------------------------------------------------------------------

PRO zspec_autonext, event

  widget_control, event.TOP, get_uvalue=state, /no_copy

  IF(state.autonext EQ 1) THEN $
    state.autonext = 0 $
  ELSE $
    state.autonext = 1

  widget_control, event.TOP, set_uvalue=state, /no_copy
END

;-----------------------------------------------------------------------
PRO zspec_goto_next, state
; Automatically move to next spectrum
; To use this routine you need to select 'Automatic next' from the
; view menu.

  COMMON results_COMMON, results, slitcount, q_index
  COMMON oned_data_COMMON, specdata, smoothspec1, smoothspec2, synspec
  COMMON twod_data_COMMON, fullfluxdata, fluxdata, extraction_r1, extraction_r2, extraction_offset, slit_display_options

  nslit = n_elements(results)

  IF(state.qonly EQ 1) THEN BEGIN
    IF(slitcount = -1) THEN slitcount=q_index[0] $
    ELSE BEGIN
      nq = n_elements(q_index)
      pos = where(q_index EQ slitcount, ct)
      IF(pos LT nq-1) THEN slitcount = q_index[pos[0]+1] $
      ELSE BEGIN
        res = dialog_message('No more slits to view', $
                             dialog_parent=event.TOP)
        return
      ENDELSE
  ENDELSE
  ENDIF ELSE IF(slitcount LT nslit-1) THEN $
    slitcount = slitcount+1 $
  ELSE BEGIN
    res = dialog_message('No more slits to view', $
                         dialog_parent=event.TOP)
    return
  ENDELSE


  widget_control, /hourglass
  IF(results[slitcount].nonlocal EQ 0) THEN $
    widget_control, state.nonlocal, sensitive=0 $
  ELSE $
    widget_control, state.nonlocal, sensitive=1

  IF(results[slitcount].nonlocal_select EQ 1) THEN BEGIN
    widget_control, state.redshift_select, $
                  set_value = string(format='(F9.5)', $
                                     results[slitcount].zresult_non.z)
    widget_control, state.redshift_select, $
                    set_droplist_select=results[slitcount].selection
    zspec_showtext, state
    redshift = $
             results[slitcount].zresult_non[results[slitcount].selection].z
          widget_control, state.nonlocal, set_value='Show local sky'
      ENDIF ELSE BEGIN
     ; Normal sky subtraction options
    widget_control, state.redshift_select, $
               set_value = string(format='(F9.5)', $
                                  results[slitcount].zresult_all.z)
    widget_control, state.redshift_select, $
                  set_droplist_select=results[slitcount].selection
    zspec_showtext, state
    redshift = $
            results[slitcount].zresult_all[results[slitcount].selection].z
          widget_control, state.nonlocal, set_value='Show non-local sky'
  ENDELSE

  specdata = zspec_1d_read(results[slitcount].nonlocal_select)
  smoothspec1 = ivarsmooth(specdata.spec,specdata.ivar,state.smooth1)
  smoothspec2 = ivarsmooth(specdata.spec,specdata.ivar,state.smooth2)
  synspec  = zspec_template(specdata.lambda,state.eigenvec)
  IF(xregistered('splot') GT 0) THEN $
     zspec_splot, state, redshift
  zspec_1d_plot, state, redshift
  fluxdata = $
            zspec_2d_read(specdata,results[slitcount].nonlocal_select)

	if slit_display_options.fullslit then begin
		fullfluxdata = zspec_2d_read(specdata,results[slitcount].nonlocal_select, /full)
	endif

  zspec_2d_plot, state, redshift
  IF(xregistered('atv') GT 0) THEN $
     zspec_update_atv, state

END




;-----------------------------------------------------------------------


FUNCTION zspec_showonly_event, event
  COMMON zspec_showonly_COMMON, selection

  widget_control, event.ID, get_uvalue=uval

  IF(uval EQ -1) THEN  $
    widget_control, event.top, /destroy $
  ELSE IF(uval EQ -2) THEN BEGIN
    selection = -1
    widget_control, event.top, /destroy
  ENDIF ELSE $
    selection = event.value

  return, 0
END


PRO zspec_showonly, event

  COMMON zspec_showonly_COMMON, selection
  COMMON results_COMMON, results, slitcount, q_index

  widget_control, event.TOP, get_uvalue=state, /no_copy


  options = [-1,0,1,2,3,4]
  IF(state.qonly EQ 1) THEN $
    selection = state.qonly_q+1 $
  ELSE $
    selection = 1



  ; Create a pop-up widget to ask the user for the Q.
  ;--------------------------------------------------------
  base = widget_base(/column, title='Select Q to view', $
                     group_leader=event.top, /modal)
  values = ['-1', $
            '0', $
            '1', $
            '2', $
            '3', $
            '4']
  q_list = CW_BGROUP(base, values, /row, /frame, $
                          label_left='Show only Q = ', $
                          event_func='zspec_showonly_event', $
                          /return_index, /exclusive, set_value=selection)
  button = widget_button(base, value='Done', uvalue=-1, $
                         event_func='zspec_showonly_event', xsize=200)
  button = widget_button(base, value='Cancel', uvalue=-2, $
                         event_func='zspec_showonly_event', xsize=200)

  widget_control, base, /realize
  xmanager, 'zspec_showonly', base


  ; Specify the selection (if one was made)
  IF(selection EQ -1) THEN $
    return

    prevq = state.qonly_q 

  state.qonly = 1
  state.qonly_q = options[selection]
  print, 'You selected Q = ', options[selection]


  q_index = where(results.zquality EQ options[selection], ct)
  print, ct, ' objects with specified Q.'

  if ct eq 0 then begin
  	null = dialog_message(['There were no objects with that specificied Zquality',$
	          'Zspec will ignore your request'])
	state.qonly = 0
	state.qonly_q = prevq
  	widget_control, event.TOP, set_uvalue=state, /no_copy
	return
  endif

  slitcount = -1
  zspec_GOTO_next, state
  widget_control, event.TOP, set_uvalue=state, /no_copy

END




;-----------------------------------------------------------------------

PRO zspec_open, event

  COMMON results_COMMON, results, slitcount, q_index
  COMMON mask_comment_COMMON, maskcomment

  widget_control, event.TOP, get_uvalue=state

  fname = dialog_pickfile(filter='*.sav', dialog_parent=event.top)
  IF(fname eq '') THEN $
    return
  restore, fname

  fname = strmid(fname,0,strlen(fname)-3) + 'fits'
  maskcomment = string(mrdfits(fname, 2, /silent))



  slitcount = -1
  zspec_GOTO_next, state
END


;-------------------------------------------------------------------------


PRO zspec_help, event

  mess_txt = string(format='(%"%s\n%s\n%s")',  $
                    'Sorry.  I have not set this up yet... ',  $
                    'please refer to the online documentation:',  $
                    'http://deep.Berkeley.EDU/~marc/deep/private/qa/zspec/zspec.html')

  res = dialog_message(mess_txt, dialog_parent=event.top)

END


;-----------------------------------------------------------------------

PRO zspec_event, event
;  Main event handling routine

  COMMON results_COMMON, results, slitcount, q_index
  COMMON oned_data_COMMON, specdata, smoothspec1, smoothspec2, synspec
  COMMON twod_data_COMMON, fullfluxdata, fluxdata, extraction_r1, extraction_r2, extraction_offset, slit_display_options
  COMMON zspec_save_COMMON, savename, cancel, dumpname

  widget_control, event.TOP, get_uvalue=state
  widget_control, event.ID, get_uvalue=uval

  CASE uval OF
    'NEXT' : BEGIN
    	temp = savename
	temp2 = dumpname
;    	savename = '~/zspec.autosave.fits'
;	dumpname = '~/zspec.autosave.sav'
        savename = 'bkup.' + temp
        dumpname = 'bkup.' + temp2
    	zspec_save_file, event
    	savename = temp
	dumpname = temp2
        nslit = n_elements(results)

        IF(state.qonly EQ 1) THEN BEGIN
          nq = n_elements(q_index)
          pos = where(q_index EQ slitcount, ct)
          IF(pos LT nq-1) THEN slitcount = q_index[pos[0]+1] $
          ELSE BEGIN
            res = dialog_message('No more slits to view', $
                                   dialog_parent=event.TOP)
            return
          ENDELSE
        ENDIF ELSE IF(slitcount LT nslit-1) THEN $
        slitcount = slitcount+1 $
        ELSE BEGIN
          res = dialog_message('No more slits to view', $
                               dialog_parent=event.TOP)
          return
        ENDELSE

        widget_control, /hourglass
        IF(results[slitcount].nonlocal EQ 0) THEN $
            widget_control, state.nonlocal, sensitive=0 $
        ELSE $
           widget_control, state.nonlocal, sensitive=1

        IF(results[slitcount].nonlocal_select EQ 1) THEN BEGIN
        ; If user has already selected a non-local sky based
        ; result then show this instead of the normal reduction
          widget_control, state.redshift_select, $
                 set_value = string(format='(F9.5)', $
                                     results[slitcount].zresult_non.z)
          widget_control, state.redshift_select, $
                            set_droplist_select=results[slitcount].selection
          zspec_showtext, state
          redshift = $
              results[slitcount].zresult_non[results[slitcount].selection].z
          widget_control, state.nonlocal, set_value='Show local sky'
        ENDIF ELSE BEGIN
        ; Normal sky subtraction options
          widget_control, state.redshift_select, $
                  set_value = string(format='(F9.5)', $
                                     results[slitcount].zresult_all.z)
          widget_control, state.redshift_select, $
                            set_droplist_select=results[slitcount].selection
          zspec_showtext, state
          redshift = $
              results[slitcount].zresult_all[results[slitcount].selection].z
          widget_control, state.nonlocal, set_value='Show non-local sky'
        ENDELSE
        specdata = zspec_1d_read(results[slitcount].nonlocal_select)
        smoothspec1 = ivarsmooth(specdata.spec,specdata.ivar,state.smooth1)
        smoothspec2 = ivarsmooth(specdata.spec,specdata.ivar,state.smooth2)
        synspec  = zspec_template(specdata.lambda,state.eigenvec)
        IF(xregistered('splot') GT 0) THEN $
          zspec_splot, state, redshift
        zspec_1d_plot, state, redshift
        fluxdata = $
                  zspec_2d_read(specdata,results[slitcount].nonlocal_select)
	if slit_display_options.fullslit then begin
		fullfluxdata = zspec_2d_read(specdata,results[slitcount].nonlocal_select, /full)
	endif

        zspec_2d_plot, state, redshift
        IF(xregistered('atv') GT 0) THEN $
          zspec_update_atv, state
    END
    'BACK' : BEGIN
    	temp = savename
	temp2 = dumpname
;    	savename = '~/zspec.autosave.fits'
;	dumpname = '~/zspec.autosave.sav'
        savename = 'bkup.' + temp
        dumpname = 'bkup.' + temp2
    	savename = temp
	dumpname = temp2
    	zspec_save_file, event
        IF(state.qonly EQ 1) THEN BEGIN
          pos = where(q_index EQ slitcount, ct)
          IF(pos GT 0) THEN slitcount = q_index[pos[0]-1] $
          ELSE BEGIN
            res = dialog_message('Cannot go back any further in file', $
                                   dialog_parent=event.TOP)
            return
          ENDELSE
        ENDIF ELSE IF(slitcount GT 0) THEN $
        slitcount = slitcount-1 $
        ELSE BEGIN
          res = dialog_message('Cannot go back any further in file', $
                               dialog_parent=event.TOP)
          return
        ENDELSE



        widget_control, /hourglass
        IF(results[slitcount].nonlocal EQ 0) THEN $
          widget_control, state.nonlocal, sensitive=0 $
        ELSE $
          widget_control, state.nonlocal, sensitive=1


        IF(results[slitcount].nonlocal_select EQ 1) THEN BEGIN
        ; If user has already selected a non-local sky based
        ; result then show this instead of the normal reduction
          widget_control, state.nonlocal, set_value='Show local sky'
          widget_control, state.redshift_select, $
                  set_value = string(format='(F9.5)', $
                                     results[slitcount].zresult_non.z)
          widget_control, state.redshift_select, $
                            set_droplist_select=results[slitcount].selection
          zspec_showtext, state
          redshift = $
              results[slitcount].zresult_non[results[slitcount].selection].z
        ENDIF ELSE BEGIN
        ; Normal sky subtraction options
          widget_control, state.nonlocal, set_value='Show non-local sky'
          widget_control, state.redshift_select, $
                  set_value = string(format='(F9.5)', $
                                     results[slitcount].zresult_all.z)
          widget_control, state.redshift_select, $
                            set_droplist_select=results[slitcount].selection
          zspec_showtext, state
          redshift = $
              results[slitcount].zresult_all[results[slitcount].selection].z
        ENDELSE
        specdata = zspec_1d_read(results[slitcount].nonlocal_select)
        smoothspec1 = ivarsmooth(specdata.spec,specdata.ivar,state.smooth1)
        smoothspec2 = ivarsmooth(specdata.spec,specdata.ivar,state.smooth2)
        synspec  = zspec_template(specdata.lambda,state.eigenvec)
        IF(xregistered('splot') GT 0) THEN $
          zspec_splot, state, redshift
        zspec_1d_plot, state, redshift
        fluxdata = zspec_2d_read(specdata,results[slitcount].nonlocal_select)
	if slit_display_options.fullslit then begin
		fullfluxdata = zspec_2d_read(specdata,results[slitcount].nonlocal_select, /full)
	endif

        zspec_2d_plot, state, redshift
        IF(xregistered('atv') GT 0) THEN $
          zspec_update_atv, state
    END
    'DONE' : BEGIN
        zspec_save, event
        widget_control, event.TOP, /destroy
        IF(xregistered('splot')) THEN splot_shutdown
        IF(xregistered('atv')) THEN atv_shutdown
        return
    END
  ENDCASE

END



;-----------------------------------------------------------------------

PRO zspec_startup, blue=blue, red=red
;  Carries out all the initialisation of the widgets etc.

  COMMON results_COMMON, results, slitcount, q_index
  COMMON oned_data_COMMON, specdata, smoothspec1, smoothspec2, synspec
  COMMON twod_data_COMMON, fullfluxdata, fluxdata, extraction_r1, extraction_r2, extraction_offset, slit_display_options

  base = widget_base(/column,title='ZSPEC  Beta v1', mbar=menu)


  ; These are the main menus.  They are pretty much all repeated in
  ; various buttons in the widget itself, but are included here for
  ; completeness.
  ;-----------------------------------------------

  ; File menu
  menu1 = widget_button(menu, value='File', /menu, xsize=200)
    button = widget_button(menu1, value='Save', $
                           event_pro='zspec_save')
    button = widget_button(menu1, value='Open', $
                           event_PRO='zspec_open')
    button = widget_button(menu1, value='Quit', uvalue='DONE')

  ; View menu
  menu2 = widget_button(menu, value='View', /menu, xsize=200)
    button = widget_button(menu2, value='Show 1d spectrum', $
                           event_PRO='zspec_splot_show')
    button = widget_button(menu2, value='Show 2d spectrum', $
                           event_pro='zspec_atv')
    button = widget_button(menu2, value='Set smoothing 1', $
                           event_PRO='zspec_set_smooth', uvalue='ONE')
    button = widget_button(menu2, value='Set smoothing 2', $
                           event_PRO='zspec_set_smooth', uvalue='TWO')
    button = widget_button(menu2, value='Jump to slit', $
                           event_PRO='zspec_jump')
    button = widget_button(menu2, value='Show only Q=?', $
                           event_PRO='zspec_showonly')
    button = widget_button(menu2, value='Show Full Slit', $
                               event_PRO='zspec_showfullslit')
    button = widget_button(menu2, value='Untilt Slit', $
				event_PRO='zspec_untiltslit')
    button = widget_button(menu2, value='Enter Mask Comments', $ 
    				event_PRO='zspec_comment_mask')


  ; Redshift menu
  menu3 = widget_button(menu, value='Redshifts', /menu, xsize=200)
    button = widget_button(menu3, value='Fit for z (preferred)', $
                           event_PRO='zspec_manualz')
    button = widget_button(menu3, value='Force z (last resort!)', $
                           event_PRO='zspec_z_by_hand')
    button = widget_button(menu3, value='Enter comment', $
                           event_PRO='zspec_comment')
    button = widget_button(menu3, value='Auto jump to next', $
                           event_PRO='zspec_autonext')

  ; Help menu
  menu4 = widget_button(menu, value='Help', /menu, xsize=200)
    button = widget_button(menu4, value='Show help', $
                           event_PRO='zspec_help')
    button = widget_button(menu4, value='Version', $
                           event_PRO='zspec_version')




  ; Line feature plots
  draw_1d_wid = lonarr(4)
  draw_2d_wid = lonarr(4)
  draw_base1 = WIDGET_BASE(base,/ROW)
  draw_base2 = WIDGET_BASE(base,/ROW)
  FOR i=0, 3 DO BEGIN
    draw_1d_wid[i] = widget_draw(draw_base1, xsize=200, ysize=200)
    draw_2d_wid[i] = widget_draw(draw_base2, xsize=200, ysize=50)
  ENDFOR


  ; Horizontal base widget
  ;-----------------------
  basex = widget_base(base, /ROW)




  ; This tab holds all the textual information available
  ; to the user.  This can be extended upon simply by adding
  ; more tabs.  The widget IDs of each text field are stored according
  ; to the names used below in the state structure.
  ;------------------------------------------------------------
  maintext = widget_text(basex, xsize=40, ysize=10, frame=1)


  basex2 = widget_base(basex, /column)






  ; Redshifts:
  ; This section consists of a drop list for selecting z's
  ; and a text list of the actual redshifts themselves.
  ; The currently selected z is stored in the state variable
  ; state.zselect
  ;-----------------------------------------------------
  base4 = widget_base(basex2, /row)
  values = ['Best', '2nd', '3rd', '4th','5th','6th','7th','8th','9th','10th']
  redshift_select = widget_droplist(base4,  value=values, $
                                  title='Redshift selection:', $
                                  event_pro='zspec_selectz', $
                                  /dynamic_resize, /frame)



  ; Non-local sky:
  ; Switches between the local and non-local sky subtracted
  ; spectra (if these are available).
  ;--------------------------------------------------------
  nonlocal = widget_button(base4, value='Show non-local sky', xsize=200, $
                           event_PRO='zspec_select_nonlocal')






  ; The user can select the redshift quality flag using these buttons.
  ; Multiple button presses simply overwrite the previous selected
  ; value
  ;-------------------------------------------------------

  base6 = widget_base(basex2, /ROW)
  values = ['Q=-2', 'Star', 'Junk', 'Q = 2', 'Q = 3', 'Q = 4']
  redshift_quality = CW_BGROUP(base6, values, /row, /frame, $
                        label_left='Quality Selection', $
                        event_func='zspec_setz', /return_index)




  ; Comments:
  ; Comments about each spectrum can be manually entered in this text
  ; field.  The easiest way to submit this info is by pressing return
  ; after typing.  However, in order to allow the user to push a
  ; botton instead, the ID of the text field is stored in state.comment
  ;--------------------------------------------------------
  button_comment = widget_button(basex2, value='Enter comment', $
                                 event_PRO='zspec_comment', $
                                 xsize=300)



  ; Smoothing:
  ; Pressing this button turns the smoothing on or off.  The level
  ; of smoothing is specified via the View menu.
  ;--------------------------------------------------------
  sbase = widget_base(basex2,/row)
  button_smooth0 = widget_button(sbase, value='Smooth off', $
                                 event_PRO='zspec_smooth', $
                                 xsize=120, uvalue='OFF')
  button_smooth1 = widget_button(sbase, value='Smooth 1', $
                                 event_PRO='zspec_smooth', $
                                 xsize=120, uvalue='ONE')
  button_smooth2 = widget_button(sbase, value='Smooth 2', $
                                 event_PRO='zspec_smooth', $
                                 xsize=120, uvalue='TWO')







  ; Forwards and Backwards:
  ; These buttons display either the next or previous slit in the mask
  ; file.
  ;--------------------------------------------------------
  base_x3 = widget_base(base, /row)
  next    = widget_button(base_x3, value='Next', uvalue='NEXT', xsize=230)
  back    = widget_button(base_x3, value='Back', uvalue='BACK', xsize=230)
  finish  = widget_button(base_x3, value='Done', uvalue='DONE', xsize=230)



  widget_control, base, /realize

  ; After the widget is realised we now need to get the widget
  ; values of each drawable area (for plotting the 1d and 2d spectra)
  ;-------------------------------------------------------
  draw_1d_idx = lonarr(4)
  draw_2d_idx = lonarr(4)
  FOR i=0, 3 DO BEGIN
    widget_control, draw_1d_wid[i], get_value=draw_1d
    draw_1d_idx[i] = draw_1d
    widget_control, draw_2d_wid[i], get_value=draw_2d
    draw_2d_idx[i] = draw_2d
  ENDFOR

  eig = zspec_1d_eigs()

  state = {base:base, $
           maintext:maintext, $ ; Textual widget for object
           draw1d:draw_1d_wid, $       ; Widget IDs of the 1d plots
           draw2d:draw_2d_wid, $       ; Widget IDs of the 2d plots
           draw_1d_idx:draw_1d_idx, $
           draw_2d_idx:draw_2d_idx, $
           zselect:0, $                ; The currently selected z index
           smooth1:3, $                 ; The smoothing currently selected
           smooth2:15, $
           smooth_on:0, $              ; Specifies whether curr on/off
           eigenvec:eig, $       ; Contains eigenvectors to create synspec
           redshift_select:redshift_select, $
;           redshift_text:redshift_text, $
           redshift_quality:redshift_quality, $
           nonlocal:nonlocal, $
           autonext:0, $
           qonly:0, $
           qonly_q:0 $
           }

	slit_display_options = {fullslit: 0, untilt: 0}


  widget_control, base, set_uvalue=state




  ; Display the first spectrum etc
  ;--------------------------------

  widget_control,/hourglass
  IF(results[slitcount].nonlocal EQ 0) THEN $
    widget_control, state.nonlocal, sensitive=0

  IF(results[slitcount].nonlocal_select EQ 1) THEN BEGIN
    widget_control, state.nonlocal, set_value='Show local sky'
    widget_control, state.redshift_select, $
                    set_value = string(format='(F9.5)', $
                                     results[slitcount].zresult_non.z)
    widget_control, state.redshift_select, $
                    set_droplist_select=results[slitcount].selection
    zspec_showtext, state
    redshift = results[slitcount].zresult_non[results[slitcount].selection].z
  ENDIF ELSE BEGIN
    widget_control, state.redshift_select, $
                    set_value = string(format='(F9.5)', $
                                     results[slitcount].zresult_all.z)
    widget_control, state.redshift_select, $
                    set_droplist_select=results[slitcount].selection
    zspec_showtext, state
    redshift = results[slitcount].zresult_all[results[slitcount].selection].z
  ENDELSE

  if n_elements(blue) ne 0 then begin
  specdata = zspec_1d_read(results[slitcount].nonlocal_select, /blue)		;CHANGED BL 7/22, if blue is set, use blue repsonse curve and fill gap
  endif else begin
  specdata = zspec_1d_read(results[slitcount].nonlocal_select)
  endelse  

  synspec  = zspec_template(specdata.lambda,state.eigenvec)
  smoothspec1 = ivarsmooth(specdata.spec,specdata.ivar,state.smooth1)
  smoothspec2 = ivarsmooth(specdata.spec,specdata.ivar,state.smooth2)
;  zspec_splot, state, redshift
  zspec_1d_plot, state, redshift
  fluxdata = zspec_2d_read(specdata,results[slitcount].nonlocal_select)

  if slit_display_options.fullslit then begin
	fullfluxdata = zspec_2d_read(specdata,results[slitcount].nonlocal_select, /full)
  endif
  zspec_2d_plot, state, redshift
  maskcomment = ''

  xmanager, 'zspec', base

END




;-------------------------------------------------------------------------

PRO zspec_splash_event, event

  COMMON splash_COMMON, proceed

  widget_control, event.id, get_uvalue=uval
  IF(uval EQ 'EXIT') THEN BEGIN
    proceed = 0
    widget_control, event.TOP, /destroy
  ENDIF ELSE BEGIN
    proceed = 1
    widget_control, event.TOP, /destroy
  ENDELSE

END


FUNCTION zspec_splash

  COMMON splash_COMMON, proceed

  base = widget_base(/column,title='Welcome')
  title = widget_label(base, value='ZSPEC Beta (v1_0)', font='lucidasans-24')
  text1 = widget_label(base, value='Programmed by: D.S Madgwick')
  text2 = widget_label(base, value='  dsm@astron.berkeley.edu  ')
  blank = widget_label(base, value='  ')
  blank = widget_label(base, value='  ')

  text3 = widget_text(base,xsize=55,ysize=16)

  var1 = getenv('D2_RESULTS')
  var2 = getenv('IDLSPEC1D_DIR')

  string01 = 'Note:'
  string02 = '-----'
  blank = ''
  string1 = 'This program makes use of two environment variables.'
  string2 = 'Your current settings for these are:'
  string3 = '  $D2_RESULTS = ' + var1
  string4 = '  $IDLSPEC1D_DIR   = ' + var2
  string5 = 'If these are not correct then please exit the program'
  string6 = 'and set them to their correct values before proceeding.'
  string7 = '  E.g.  setenv IDLSPEC1D_DIR ~dsm/cvs/deep/spec1d'
  string8 = 'If you need more information please refer to the '
  string9 = 'documentation for this program.'

  widget_control, text3, set_value=string01
  widget_control, text3, set_value=string02, /append
  widget_control, text3, set_value=blank, /append
  widget_control, text3, set_value=string1, /append
  widget_control, text3, set_value=string2, /append
  widget_control, text3, set_value=blank, /append
  widget_control, text3, set_value=string3, /append
  widget_control, text3, set_value=string4, /append
  widget_control, text3, set_value=blank, /append
  widget_control, text3, set_value=string5, /append
  widget_control, text3, set_value=string6, /append
  widget_control, text3, set_value=blank, /append
  widget_control, text3, set_value=string7, /append
  widget_control, text3, set_value=blank, /append
  widget_control, text3, set_value=string8, /append
  widget_control, text3, set_value=string9, /append


  done  = widget_button(base, value='Exit zspec', uvalue='EXIT', xsize=60)
  but2  = widget_button(base, value='Continue', uvalue='CONT', xsize=60)

  widget_control, base, /realize
  xmanager, 'zspec_splash', base

  return, proceed
END

;---------------------------------------------------------------------------
;---------------------------------------------------------------------------



PRO zspec_test, mname, local_only=local_only, nodeep=nodeep, orelse=orelse, blue=blue, _extra=extra, red=red
; Main routine

; check if this is a KTRS mask. if so, set nodeep and local_only.
 mname = strcompress(mname, /rem)
  if strlen(mname) ge 4 then begin
      if strupcase(strmid(mname, 0, 4)) eq 'KTRS' or $
         strupcase(strmid(mname, 0, 5)) eq 'GOODS' then begin
          local_only = 1
          nodeep = 1
      endif
  endif

 ; check if this is a ORELSE mask, if so, set nodeep and local_only
 if n_elements(orelse) gt 0 then begin 
 	local_only=1
  	nodeep=1
 endif
 if n_elements(blue) gt 0 then begin
 	blue=1
 endif
; check local_only and notdeep keywords.
  if n_elements(local_only) gt 0 then local_only = local_only[0] ge 1 $
  else local_only = 0
  nlsky = 1 - local_only
  if n_elements(nodeep) gt 0 then nodeep =  nodeep[0] ge 1 $
  else nodeep =  0



  !x.style=1
  COMMON results_common, results, slitcount
  COMMON directory_COMMON, fulldir
  COMMON zspec_save_COMMON, savename, cancel, dumpname
  COMMON zspec_lines_COMMON, linenames, linewaves, telluricnames, telluricwaves

  linenames = ['Lya', 'CIV', 'FeII', 'FeII','MGII', 'NeV','OII', 'OII', $
               'OII', 'NeIII', 'HeI', 'CaK', 'CaH', $
               'He', 'Hd', 'G Band', 'Hg',  'Hb', $
               'OIII', 'OIII', 'Mgb', 'NaD', 'OI', 'NII', $
               'Ha', 'NII', 'SII', 'SII', 'Ca', $
               'Ca', 'Ca']
; note that LyA and CIV are not in the spec1d template. for this
; reason, the manual z code sets these wavelength values to "N/A" when
; displayed in the "Fit for z" widget.
  linewaves = [1215.67, 1549.5, 2379.0, 2593.0, 2799.11, 3426., 3727.55, 3726.03, $
               3728.82, 3869., 3889., 3933.7, 3968.5, $
               3970.07, 4101.73, 4303.6, 4340.5, 4861.3, $
               4958.9, 5006.9, 5173., 5893.0, 6300., 6548.1, $
               6562.8, 6583.4, 6716.4, 6730.8, 8498.062, $
               8542.144, 8662.170]

  telluricnames = ['B-Band', 'A-Band']
  telluricwaves = [6860, 7600]

  ; Check if the maskname to use has been specified, if not
  ; then exit program
  ;------------------------------------------
  IF (N_PARAMS() LT 1) THEN BEGIN
    mess_txt = 'You need to specify the mask name to run zspec.  E.g. IDL> zspec, 3206'
    res = dialog_message(mess_txt,/error)
   return
  ENDIF




  proceed = zspec_splash()
  IF(proceed EQ 0) THEN return

  topdir = getenv('D2_RESULTS')

  ; Now we start to read in the data
  ;-----------------------------------
  if nodeep eq 1 then maskname = strcompress(mname,  /remove_all) $
   else maskname = STRING(mname,FORMAT='(I4.4)')

  path = topdir+'/zresult'

  filelist = findfile(path + '/zresult.' + maskname $
                      + '*.fits*', count=ct)
  IF ct NE 1 THEN BEGIN
    IF ct LT 1 THEN $
      mess_txt = 'ERROR: No zresult file found, please select the one you would like to use.'
    IF ct gt 1 then $
      mess_txt = string(format='(%"%s\n%s")', $
              'There were more than one zresults files corresponding to', $
              'your mask.  Please select the one you would like to use.')

    res = dialog_message(mess_txt)
    fname = dialog_pickfile(filter='*.fits',path=concat_dir(topdir,'/zresult'))
  endif else $
    fname = filelist[0]

  zbest = mrdfits(fname, 1, /silent) ; 1st HDU has the best fits
  zall  = mrdfits(fname, 2, /silent) ; 2nd HDU has everything
  if nlsky then znonl = mrdfits(fname, 4, /silent) ; Non-local sky results


  ; We use the obj_info file to automatically determine the
  ; directory name for the slit files etc.  If this file isn't found
  ; or there are more then one option, the user can do this manually.
  ;----------------------------------------------

  filelist=findfile(topdir +'/' + maskname + '/*/obj_info*.fits', count=nobjfile)
  ;filelist=findfile(topdir +'/' + maskname + '/obj_info*.fits', count=nobjfile)
  IF nobjfile NE 1 THEN BEGIN
    mess_txt = string(format='(%"%s\n%s\n%s\n%s\n%s")', $
    'I was not able to automatically determine the path to the ',  $
    'slit* and spec1d* files, possibly because there were more ', $
    'than one possibilities.  Would you like to do this manually? ', $
    '(All you need to do is select any fits file in the directory', $
    '  you want me to use)')
    res = dialog_message(mess_txt,/question)
    IF(res EQ 'Yes') THEN begin
      res = dialog_pickfile(get_path=fulldir, filter='*.fits', path=topdir, $
                            title='Pick any file in the directory you want to use')
      print, fulldir
    ENDIF ELSE $
      return
  ENDIF else begin
    filelist=filelist[0]
    dirlist = strsplit(filelist, '/', /extract)
    ndir = n_elements(dirlist)
    date = dirlist[ndir-2]
    fulldir = concat_dir(concat_dir(topdir,maskname),date)
    print, fulldir
  ENDELSE





  ; We now have a list of the mask slits (in zbest), and we have the
  ; directory name for all corresponding spec1d files etc.
  ;--------------------------------------------------------

  slitcount = 0
  slitlist = zbest.slitname
  slitnames = string(slitlist, FORMAT='(I3.3)')
  nslit = N_ELEMENTS(zbest.slitname)

  objlist  = zbest.objname



  ; Remove serendips from the list
  ;--------------------------------

  serendip_count = 0
  FOR i=0, nslit-1 DO BEGIN
    IF(strmid(objlist[i],0,1) EQ 's') THEN $
      serendip_count = serendip_count + 1
  ENDFOR
  print, 'No of serendips: ', serendip_count

  serendip = 0
  objnames = strarr(nslit - serendip_count)
  slnames  = strarr(nslit - serendip_count)
  FOR i=0, nslit-1 DO BEGIN
    IF(STRMID(objlist[i],0,1) NE 's') THEN BEGIN
      objnames[serendip] = objlist[i]
      slnames[serendip] = slitnames[i]
      serendip = serendip + 1
    ENDIF
  ENDFOR

  slitnames = slnames
  nslit = N_ELEMENTS(objnames)



  ; Remove allignment stars and sky slits
  ;-----------------------------------------
  counter = 0
  aln_stars =  0
  sky_stuff =  0
  objnames_new =  STRARR(nslit)
  slitnames_new =  STRARR(nslit)

  info_file = findfile(fulldir+'/obj_info*.fits', count=ct)

  IF(ct NE 1) THEN BEGIN
    mess_txt = string(format='(%"%s\n%s\n%s\n%s")', $
    'I could not find the obj_info*.fits file which specifies which', $
    'slits are guide stars and sky slits.  Would you like to specify', $
    'it manually?  If you choose not to then these objects will be', $
    'included in the objects viewed by zspec.')
    res = dialog_message(mess_txt,/question)
    IF(res EQ 'Yes') THEN BEGIN
      info_file = dialog_pickfile(path=fulldir,filter='*.fits')
      ct = 1
      info = [info_file]
    ENDIF
  ENDIF
  IF(ct EQ 1) THEN BEGIN
    info =  MRDFITS(info_file[0], 1, /silent)
    FOR i=0,  nslit-1 DO BEGIN
      cnt = WHERE(info.objno EQ LONG(objnames[i]), c)
      IF c NE 2 THEN PRINT, 'Wrong number of entries in infofile for obj: ',  objnames[i]
      IF c EQ 0 THEN BEGIN
          print, 'Error:  No entry in obj_info for: ', objnames[i], '!!!'
          print, '   ... object skipped.'
          CONTINUE
      ENDIF
      obj_type =  info[cnt[0]].objtype
      IF obj_type EQ 'P' THEN BEGIN
        objnames_new[counter] =  objnames[i]
        slitnames_new[counter] =  slitnames[i]
        counter =  counter+1
      ENDIF
      IF obj_type EQ 'S' THEN BEGIN
        sky_stuff= sky_stuff + 1
        PRINT,  objnames[i],  ' is sky'
      ENDIF
      IF obj_type EQ 'A' THEN BEGIN
        aln_stars= aln_stars + 1
        PRINT,  objnames[i],  ' is a star'
      ENDIF
    ENDFOR

    PRINT,  aln_stars,  ' allignment stars removed.'
    PRINT,  sky_stuff,  ' sky slits removed'

    cnt =  WHERE(objnames_new NE '',  c)
    objnames =  objnames_new[cnt]
    slitnames =  slitnames_new[cnt]
    nslit =  N_ELEMENTS(objnames)
    PRINT,  c,  ' objects to be redshifted'
  ENDIF

  if n_elements(blue) ne 0 then begin		;ADDED BL 7/22, one time warning about which response curve is being used to display 1d data 
 	print, 'Warning: using blue fill gap and response correction, central wavelength should be less than 7300 A'
  endif else begin
	if n_elements(red) ne 0 then begin
		print, 'Warning: using reallyred fill gap and response correction, central wavelength should be greater than 7900 A'
	endif else begin
		print, 'Warning: using standard red fill gap and response correction, central wavelength should be between 7300 A and 7900 A'
	endelse
  endelse 

; check for each slit if there is more than 1 object in the slit.
  zspec_set_mobject, slitnames


  ; For each object we determine if it has a spectrum with non-local
  ; sky subtraction.  The user can select to view the corresponding
  ; spectra for these objects.
  ; To see if a spectrum has nonlocal I use the total degrees of
  ; freedom of the template fits.  If this fit has been performed then
  ; we would expect this to be non-zero.
  ;------------------------------------------------------
  nonlocallist = intarr(nslit)
  znonl2 = zall
  FOR i=0, nslit-1 DO BEGIN
    if nlsky then $
      cnt1 = where(strtrim(znonl.objname) EQ strtrim(objnames[i]), ct1)
    cnt2 = where(strtrim(zall.objname) EQ strtrim(objnames[i]), ct2)
    if nlsky then check = total(znonl[cnt1].dof) else check = 0
    IF check GT 0 THEN BEGIN
      znonl2[cnt2] = znonl[cnt1]
      nonlocallist[i] = 1
    ENDIF ELSE $
      nonlocallist[i] = 0
  ENDFOR


  ; All of the data we want to eventually store is contained in zall
  ; we simply need to pick which bits of it are the correct values.
  ; We create a structure of these correct values (and other useful
  ; information) that will be the main tool for navigating through
  ; the file and storing z qualities etc.
  ;-----------------------------------------------------------

  z1  = zbest[0]  ; An example zresults structure
  z10 = replicate(z1,10)
  stuff = {maskname:maskname, $
           slitname:'', $
           objname:'', $
           selection:0, $       ; 0-9 selection from zall
           nonlocal:0, $        ; 1 if a non-local spec exists
           nonlocal_select:0, $  ; 1 if the non-local spec is selected
           zmanual:0, $         ; 1 if a manual-z has been determined
           zquality:0, $
           qindex:0, $
           comment:'', $
           zresult:z1, $       ; This will hold the best structure
           zresult_all:z10, $  ; This holds all 10 structures from spec1d
           zresult_non:z10 $   ; This holds the 10 non-local structures
          }

  results = replicate(stuff, nslit)
  FOR i=0, nslit-1 DO BEGIN
    results[i].slitname = slitnames[i]
    results[i].objname  = strtrim(objnames[i])
    dex = where(strtrim(zbest.objname) eq strtrim(objnames[i]), num)
    if num ne 1 then message, 'ERROR with zbest variable!' $
    else results[i].zresult = zbest[dex[0]] 
;    results[i].zresult  = zbest[i]
    cnt = where(strtrim(zall.objname) EQ strtrim(objnames[i]), ct)
    IF ct EQ 10 THEN results[i].zresult_all = zall[cnt] $
    ELSE BEGIN
      print, 'Error: missing entries for object', objnames[i]
      results[i].zresults_all[0:ct-1]=zall[cnt]
    ENDELSE
    IF nonlocallist[i] EQ 1 THEN BEGIN
      results[i].nonlocal = 1
      cnt = where(strtrim(znonl2.objname) EQ strtrim(objnames[i]), ct)
      results[i].zresult_non[0:9] = znonl2[cnt]
  ENDIF ELSE $
      results[i].zresult_non = zall[cnt]
  ENDFOR



; define the date to be used in the output file.
  dex = where(strlen(strcompress(zbest.date, /rem)) eq 10, ngood)
  if ngood eq 0 then message, 'No good dates in zresult file!'
  savdate = strcompress(zbest[dex[0]].date, /rem)
  savename = 'zspec.' + getenv('LOGNAME') + '.' +  $
    maskname + '.' + savdate + '.fits'    ; name of the output file
  dumpname = 'zspec.' + getenv('LOGNAME') + '.' + $
    maskname + '.' + savdate + '.sav'           ; name of the output file


  ; Start the widget
  ;-------------------
  if n_elements(blue) then begin	;ADDED BL 7/22, allows for blue response curve 
  zspec_startup, /blue
  endif else begin
  zspec_startup
  endelse

  cnt = where(results.zquality EQ -1, cn1)
  cnt = where(results.zquality EQ 0, c0)
  cnt = where(results.zquality EQ 1, c1)
  cnt = where(results.zquality EQ 2, c2)
  cnt = where(results.zquality EQ 3, c3)
  cnt = where(results.zquality EQ 4, c4)

  IF(cn1 GT 0) THEN print, cn1, ' objects with Q = -1'
  IF(c0 GT 0)  THEN print, c0, ' objects with Q = 0'
  IF(c1 GT 0)  THEN print, c1, ' objects with Q = 1'
  IF(c2 GT 0)  THEN print, c2, ' objects with Q = 2'
  IF(c3 GT 0)  THEN print, c3, ' objects with Q = 3'
  IF(c4 GT 0)  THEN print, c4, ' objects with Q = 4'

  good = c3+c4
  all = c1+c2+c3+c4
  IF(all GT 0) THEN BEGIN
      perc = string(format='(F5.1)', float(good)/all * 100.0)
      print, '-------------------------'
      print, '  Completeness of galaxy redshifts checked: ', perc, '%'
  ENDIF

END








