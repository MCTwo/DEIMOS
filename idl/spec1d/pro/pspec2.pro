;+
; NAME:
;  pspec2
;
; PURPOSE:
;  Routine for plotting spectra from DEIMOS 1D spectro outputs
;
; CALLING SEQUENCE:
;  pspec2, maskname,[/curdir], [/noremove], [oldfile='file', 
;                                                zquality={0,1,2,3,4}]  
;
; INPUTS:
;  maskname  -   mask name (string or integer)
;        NOTE: if maskname is forgotten the program defaults to 3206
;
;
;  Other routines contained in this file are:
;    MY_REBIN()    - Rebins the observed spectrum, so that it is easier
;                    to see what's going on.
;    MY_PLOT_SPEC  - Opens an splot window for viewing the entire spectrum,
;                    also makes separate plots of the main emission/
;                    absorption features that should be present.
;    MY_SHOW_TEXT  - Prints out the main information for each object,
;                    including the redshift, reduced chi**2 and eigenvalues
;    MY_PRINT_ERROR - Simple error writer
;
; KEYWORDS:
;    curdir  --- set this keyword if you want to read spec1d and slit
;                files from current directory.  In this case the
;                zresult subdirectory should be under the current directory.
;    noremove -- stops the program from trying to read a obj_info
;                file, which is only available for v1_0 and later
;                versions of the pipeline
;    oldfile --- allows you to read in the results from a previously 
;                saved session of pspec2 for that mask.  Simply specify
;                the name of the file you created at the end of your
;                last pspec2 session.
;    zquality -- set this keyword to an integer value corresponding to
;                the Q of the galaxies you'd like to review in this mask
;                Note that without the `oldfile' keyword this quantity is
;                of course redundant.  Having zquality set should allow you
;                to review (and alter) only those galaxies which currently
;                have Q=zquality.
;
; COMMENTS:
;
;  The main routine is PSPEC2, this decides on the maskfile to open and then
;  creates a widget.  It then plots out the first spectrum in the 
;  requested maskfile and then switches over control to the event manager.
;
;  NOTE: All of the spectra and associated information is read with the
;  DEIMOS_READSPEC and DEIMOS_READOBJ routines.
;
;  The next most important routine is PSPEC2_EVENT which controls what
;  happens when a button is pushed.  Note that the splot window is 
;  essentially independent of the main window, so you can do whatever you
;  like with it without upsetting the main program (which will generate
;  a new splot window when you select another spectrum).
;
;
;
; FUTURE WORK:  (still to be done)
;  2)  Plot 2D spectrum (perhaps call an ATV window)
;  3)  Add quality selection buttons so that the program creates an
;      output file with final redshifts and their associated reliabilities
;  4)  Add some mechanism for adding comments to the output files.
;  5)  Add the chi**2(lag) plots somewhere so that the user can visually
;      compare various redshifts.
;
;
;
; REVISION HISTORY
;  v0.1  DSM  18-Oct-02  - based upon pspec.pro (D.Finkbeiner and D.Schlegel)
;--------------------------------------------------------



FUNCTION read_z, maskname
;  Read in all redshifts for a given mask.
;*************

common pspec2_dirs, topdir0, fulldir

;  topdir = getenv('D2_RESULTS')
;  if topdir eq '' then print, 'You should set $D2_RESULTS'
  topdir = getenv('D2_RESULTS')

  filelist = findfile(topdir+'/zresult/zresult.'+maskname+'*.fits*', count=ct)
  if ct gt 1 then begin
     print, 'More than one zresult file for that mask'
     print, filelist
     read, 'Enter file number you want (starting with 0)', num
     fname = filelist[num]
     zresfile = num
  endif else fname = filelist[0]
; fname = 'zresult.'+maskname+'.*fits'
  zdir = concat_dir(topdir, 'zresult')
  zall = mrdfits(fname, 2, /silent) ; 2nd HDU has everything 


  return, zall

END




PRO select_z, zall, objnm, zres=zres, synspec=synspec, $
              index=index, lambda=lambda
  common pspec2_dirs, topdir0, fulldir

;  Return all redshifts for a *given* object number.
;
;  If you want to return a best-fitting synspec then include index number
;  and lambda, as well as synspec keyword.
;
;**********************

  ; Remove serendips
  objlist = zall.objname
  FOR i=0, N_ELEMENTS(objlist)-1 DO BEGIN
    IF( STRMID(objlist[i],0,1) EQ 's') THEN $
       objlist[i] = '0'
  ENDFOR

  w = where(long(objlist) EQ long(objnm), nslit)
  if nslit NE 5 then begin
;     print, 'Wrong no. of redshifts, object number ', objnm
;     class = zall[w].class
;     cnt = WHERE(class EQ 'GALAXY',c)
;     w = w[cnt]
  endif 
;  zres = zall[w[0:4]]   ; take the first 5 z's for now
  zres = zall[w]

  ; Create an array of templates
  ;++++++++++++++++++++++++++

  deep_dir = getenv('DEEP_DIR')
  if deep_dir eq '' then print, 'You should set $DEEP_DIR'

  files  = zres.tfile
  tcols  = zres.tcolumn
  thetas = zres.theta
  npolys = zres.npoly
  reds   = zres.z


  IF arg_present(synspec) THEN BEGIN
    file  = files[index]
    tcol  = tcols[*,index]
    theta = thetas[*,index]
    npoly = npolys[index]
    red   = reds[index]
 
    eigen_name = concat_dir(concat_dir(deep_dir, 'spec1d/templates'), $
                             file)
    eigen_name = strcompress(eigen_name, /rem)
    eig = mrdfits(eigen_name, 0, h, /silent)
    w = where(tcol GE 0, nw)
    syn = 0

  ; Temporary fix to deal with stars
  IF nw LT 1 THEN BEGIN
    nw = 1
    w = [0]
    tcol[w[0]] = 10
  ENDIF

     for i=0, nw-1 do begin
       jj = tcol[w[i]]
       syn = syn+theta[i]*eig[*, jj]
     endfor 
     IF npoly gt 0 THEN begin ;add in the polynomial terms
        npoints = n_elements(eig[*, 0])
        parray = poly_array(npoints, npoly)
        for i=0, npoly-1 do syn = syn + theta[nw+i]*parray[*, i]
     ENDIF

     coeff0 = sxpar(h, 'COEFF0')
     coeff1 = sxpar(h, 'COEFF1')
     
     elambda = 10.d^(coeff0+dindgen((size(eig, /dimens))[0])*coeff1)
; ??? interpol is pretty stupid - use combine1fiber in idlspec2d (SDSS)
     synspec = interpol(syn, elambda, lambda/(red+1))
  ENDIF


END





PRO MY_SHOW_TEXT, zans, zres
;
;  Display the relevant information for this spectrum
;
;************************

COMMON TEXT_STUFF, text_widg, text_widgz, pgal
COMMON SPEC_STUFF, maskname, slitnames, slitcount, nslit, objnames, zall

  fir_line = 'Mask name: ' + zans.maskname + '  Slitname: ' + zans.slitname
  sec_line = 'Object name: ' + zans.objname + '  P(gal) = ' + STRING(pgal[slitcount],FORMAT='(F5.2)')
  thr_line = 'Redshift: ' + STRING(zans.z, FORMAT='(F9.5)') + ' +/- ' + STRING(zans.z_err, FORMAT='(E9.1)')
  for_line = 'Delta chi^2: ' + STRING(zans.rchi2diff * zans.dof)
  fif_line = 'Eigenvalues: ' + STRING(zans.theta[0:4])
  six_line =  'Velocity dispersion: ' + STRING(zans.vdisp)
  sep_line = '-------------------------------------------------'


  redshifts = zres.z
  
  WIDGET_CONTROL, text_widg, SET_VALUE=fir_line
  WIDGET_CONTROL, text_widg, SET_VALUE=sec_line, /APPEND
  WIDGET_CONTROL, text_widg, SET_VALUE=thr_line, /APPEND
  WIDGET_CONTROL, text_widg, SET_VALUE=for_line, /APPEND
  WIDGET_CONTROL, text_widg, SET_VALUE=six_line, /APPEND
  WIDGET_CONTROL, text_widg, SET_VALUE=sep_line, /APPEND
;  WIDGET_CONTROL, text_widg, SET_VALUE=fif_line, /APPEND

  WIDGET_CONTROL, text_widg, SET_VALUE='Eigenvalues--absorption and emission: '+ STRING(zans.theta[0:1], FORMAT='(2F9.1)'), /append
  WIDGET_CONTROL, text_widg, SET_VALUE='Eigenvalues--continuum: ' + $
    STRING(zans.theta[2:4], FORMAT='(3F9.1)') , /append

;  WIDGET_CONTROL, text_widg, SET_VALUE=sep_line, /APPEND

 WIDGET_CONTROL, text_widgz, SET_VALUE=STRING(redshifts, FORMAT='(F5.3,9F7.3)')


END
;*****************************************************







PRO MY_PRINT_ERROR, err_text
;
;  Display a simple error message in the text window
;
;*****************

COMMON TEXT_STUFF, text_widg

  WIDGET_CONTROL, text_widg, SET_VALUE=err_text

END







;******************************************************
PRO MY_PLOT_SPEC
;
; Plot a spectrum and its variance in the main plotting window
;
;
COMMON DRAW_STUFF, drawID1, drawID2, drawID3, drawID4
COMMON SPECTRUM, lambda, spec, ivar, synspec, redshift, oldspec, oldsynspec, nsmooth
COMMON SPEC_STUFF, maskname, slitnames, slitcount, nslit, objnames, zall

  rlambda = lambda/(1.+redshift)




; initialize colors, just in case
; 0 - black
; 1 - red
; 2 - green
; 3 - blue
; 4 - cyan
; 5 - magenta
; 6 - yellow
; 7 - white
  
; the following line is voodoo, necessary to make IDL actually use the
; color table.

  device, decomposed = 0

  loadct, 0

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


  ; OII emission feature
  ;+++++++++++++++++++++++
  WSET, drawID1

  cnt = WHERE(rlambda LE 3747 AND rlambda GE 3707, c)
  IF c GE 1 THEN BEGIN
   PLOT, lambda[cnt], spec[cnt], title='O[II]', CHARSIZE=0.8, $
      BACKGROUND = 7, COLOR = 0
   OPLOT, lambda[cnt], synspec[cnt], COLOR=3;DJS_ICOLOR('blue')
  ENDIF
  IF c LT 1 THEN $
   PLOT, [0,1], [0,1]  , CHARSIZE=0.8, $
      BACKGROUND = 7, COLOR = 0


  ; H and K absorption features
  ;++++++++++++++++++++++++++++
  WSET, drawID2
  cnt = WHERE(rlambda LE 4000 AND rlambda GE 3900, c)
  IF c GE 1 THEN BEGIN
   PLOT, lambda[cnt], spec[cnt], title='H and K', CHARSIZE=0.8, $
      BACKGROUND = 7, COLOR = 0
   OPLOT, lambda[cnt], synspec[cnt], COLOR=3
  ENDIF
  IF c LT 1 THEN $
   PLOT, [0,1], [0,1], CHARSIZE=0.8, $  
      BACKGROUND = 7, COLOR = 0


 IF redshift GT 0.3 THEN BEGIN
  ; Hbeta feature
  ;++++++++++++++++
  WSET, drawID3
  cnt = WHERE(rlambda LE 4880 AND rlambda GE 4840, c)
  IF c GE 1 THEN BEGIN
   PLOT, lambda[cnt], spec[cnt], title='Hbeta', CHARSIZE=0.8, $
      BACKGROUND = 7, COLOR = 0
   OPLOT, lambda[cnt], synspec[cnt], COLOR=3
  ENDIF
  IF c LT 1 THEN $
   PLOT, [0,1], [0,1], title='Hbeta' , CHARSIZE=0.8, $ 
      BACKGROUND = 7, COLOR = 0

  ; O[III] emission feature
  ;+++++++++++++++++++++++++
  WSET, drawID4
  cnt = WHERE(rlambda LE 5030 AND rlambda GE 4940, c)
  IF c GE 1 THEN BEGIN
   PLOT, lambda[cnt], spec[cnt], title='O[III]', CHARSIZE=0.8, $
      BACKGROUND = 7, COLOR = 0
   OPLOT, lambda[cnt], synspec[cnt], COLOR=3
  ENDIF
  IF c LT 1 THEN $
   PLOT, [0,1], [0,1], title='O[III]'  , CHARSIZE=0.8, $
      BACKGROUND = 7, COLOR = 0
 ENDIF ELSE BEGIN

  ; Halpha feature
  ;++++++++++++++++
  WSET, drawID3
  cnt = WHERE(rlambda LE 6600 AND rlambda GE 6540, c)
  IF c GE 1 THEN BEGIN
   PLOT, lambda[cnt], spec[cnt], title='Halpha', CHARSIZE=0.8, $
      BACKGROUND = 7, COLOR = 0
   OPLOT, lambda[cnt], synspec[cnt], COLOR=3
  ENDIF
  IF c LT 1 THEN $
   PLOT, [0,1], [0,1], title='Hbeta' , CHARSIZE=0.8, $ 
      BACKGROUND = 7, COLOR = 0

  ; S emission feature
  ;+++++++++++++++++++++++++
  WSET, drawID4
  cnt = WHERE(rlambda LE 6745 AND rlambda GE 6700, c)
  IF c GE 1 THEN BEGIN
   PLOT, lambda[cnt], spec[cnt], title='[SII]', CHARSIZE=0.8, $
      BACKGROUND = 7, COLOR = 0
   OPLOT, lambda[cnt], synspec[cnt], COLOR=3
  ENDIF
  IF c LT 1 THEN $
   PLOT, [0,1], [0,1], title='O[III]'  , CHARSIZE=0.8, $
      BACKGROUND = 7, COLOR = 0


 ENDELSE




  ; Rebin the spectrum (don't need this much resolution!!)
  ;+++++++++++++++++++++++++++++++++++++++++++++++++++++++
;  spec    = 5.0*MY_REBIN(spec,5)
;  lambda  = MY_REBIN(lambda,5)
;  synspec = 5.0*MY_REBIN(synspec,5)


  ; Reset the plotting limits to fit variance etc in bottom of plot
  ;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;  mins = MIN(spec)
;  maxs = MAX(spec)
  medspec = median(spec)
  stdspec = stddev(spec)
  maxs = medspec + 3.*stdspec ;set limits not so affected by bad points
  mins = medspec - 3.*stdspec
  mins_new = mins - 0.1*(maxs-mins)
  SPLOT, lambda, spec, xtitle='!6Observed Wavelength (Ang)', $
    ytitle='Counts',  yrange=[mins_new,maxs], title='Object: '+objnames[slitcount]
  SOPLOT, lambda, synspec, COLOR=3

  ; Overplot the variance in the bottom
  ;+++++++++++++++++++++++++++++++++++++

  variance = ivar*0.
  good = where(ivar gt 0.)
  variance[good] = 1./ivar[good] <  1.e5
  maxv = MAX(variance)
  minv = MIN(variance)
;  variance = 0.1*(maxs-mins)/(maxv-minv)*variance + mins
  variance = 0.5*(maxs-mins)/(maxv-minv)*variance + mins
  SOPLOT, lambda, variance, color=1;DJS_ICOLOR('red')
  min_lambda = MIN(lambda)
  max_lambda = MAX(lambda)

  oii = 3726.05*(1.+redshift)
  oii2 =  3728.8*(1.+redshift)
  k   = 3933.4*(1.+redshift)
  h   = 3969.17*(1.+redshift)
  hd  = 4102.8*(1.+redshift)
  g   = 4303.36*(1.+redshift)
  hg  = 4340*(1.+redshift)
  hb  = 4861.32*(1.+redshift)
  o31 = 4958.911*(1.+redshift)
  o32 = 5006.843*(1.+redshift)
  mg  = 5174*(1.+redshift)
  nad = 5892.36*(1.+redshift)
  ha  = 6562.801*(1.+redshift)
  sii = 6716.44*(1.+redshift)
    
  y = [mins_new,maxs]
  IF oii GT min_lambda AND oii LT max_lambda THEN $
  SOPLOT, [oii,oii], y, LINESTYLE=2, color=6;DJS_ICOLOR('yellow')
  IF oii2 GT min_lambda AND oii2 LT max_lambda THEN $
  SOPLOT, [oii2,oii2], y, LINESTYLE=2, color=6;DJS_ICOLOR('yellow')
  IF k GT min_lambda AND k LT max_lambda THEN $
  SOPLOT, [k,k], y, LINESTYLE=2, color=2;DJS_ICOLOR('green')
  IF h GT min_lambda AND h LT max_lambda THEN $
  SOPLOT, [h,h], y, LINESTYLE=2, color=2;DJS_ICOLOR('green')
  IF hd GT min_lambda AND hd LT max_lambda THEN $
  SOPLOT, [hd,hd], y, LINESTYLE=2, color=6;DJS_ICOLOR('yellow')
  IF g GT min_lambda AND g LT max_lambda THEN $
  SOPLOT, [g,g], y, LINESTYLE=2, color=2;DJS_ICOLOR('green')
  IF hg GT min_lambda AND hg LT max_lambda THEN $
  SOPLOT, [hg,hg], y, LINESTYLE=2, color=6;DJS_ICOLOR('yellow')
  IF hb GT min_lambda AND hb LT max_lambda THEN $
  SOPLOT, [hb,hb], y, LINESTYLE=2, color=6;DJS_ICOLOR('yellow')
  IF o31 GT min_lambda AND o31 LT max_lambda THEN $
  SOPLOT, [o31,o31], y, LINESTYLE=2, color=6;DJS_ICOLOR('yellow')
  IF o32 GT min_lambda AND o32 LT max_lambda THEN $
  SOPLOT, [o32,o32], y, LINESTYLE=2, color=6;DJS_ICOLOR('yellow')
  IF mg GT min_lambda AND mg LT max_lambda THEN $
  SOPLOT, [mg,mg], y, LINESTYLE=2, color=2;DJS_ICOLOR('green')
  IF nad GT min_lambda AND nad LT max_lambda THEN $
  SOPLOT, [nad,nad], y, LINESTYLE=2, color=2;DJS_ICOLOR('green')
  IF ha GT min_lambda AND ha LT max_lambda THEN $
  SOPLOT, [ha,ha], y, LINESTYLE=2, color=6;DJS_ICOLOR('yellow')
  IF sii GT min_lambda AND sii LT max_lambda THEN $
  SOPLOT, [sii,sii], y, LINESTYLE=2, color=6;DJS_ICOLOR('yellow')
  


  SXYOUTS, oii2, 0.9*maxs, '[OII]'
  SXYOUTS, k, 0.9*maxs, 'K'
  SXYOUTS, h, 0.9*maxs, 'H'
  SXYOUTS, hd, 0.9*maxs, 'Hd'
  SXYOUTS, g, 0.9*maxs, 'g'
  SXYOUTS, hg, 0.9*maxs, 'Hg'
  SXYOUTS, hb, 0.9*maxs, 'Hb'
  SXYOUTS, o31, 0.9*maxs, '[OIII]'
  SXYOUTS, mg, 0.9*maxs, 'MgB'
  SXYOUTS, nad, 0.9*maxs, 'NaD'
  SXYOUTS, ha, 0.9*maxs, 'Ha'
  SXYOUTS, sii, 0.9*maxs, '[SII]'


END
;*************************************************










PRO pspec2_event, ev
; 
;  Event handling procedure
;
;****************************************

COMMON RESULTS, fin_redshifts, fin_quality, fin_comments, fin_canned, finprevious_data
COMMON SPEC_STUFF, maskname, slitnames, slitcount, nslit, objnames, zall
COMMON SPECTRUM,  lambda, spec, ivar, synspec, redshift, oldspec, oldsynspec, nsmooth
COMMON BUTTON_LABELS, redshift_list, canned_list ; Must update the list names!
COMMON SKIP_STUFF, indx_pos, zqq, indx_cnt


; Set the graphics device to our widget
 WIDGET_CONTROL, ev.top, GET_UVALUE=drawID1
 WSET, drawID1

; Which button has been pressed
;******************************
 WIDGET_CONTROL, ev.id, GET_UVALUE=uval



; handle the possible values of 'uval' via a case statement
 CASE uval of 


;*************************   NEXT BUTTON  *******************************
; Plot the next spectrum
;************
 'NEXT':  begin
  IF zqq EQ 0 THEN BEGIN
    slitcount = slitcount + 1    ; Move to the next slit
    IF slitcount GE nslit THEN BEGIN
      slitcount = slitcount - 1
      GOTO, ERROR1
    ENDIF
  ENDIF

   ; Check if we're only looking at certain Q's
   IF zqq EQ 1 THEN BEGIN
    indx_cnt = indx_cnt + 1
    IF indx_cnt GE nslit THEN BEGIN
      indx_cnt = indx_cnt + 1
      GOTO, ERROR1
    ENDIF
    slitcount = indx_pos[indx_cnt]
   ENDIF

 
   DEIMOS_READOBJ, maskname, objnames[slitcount], lambda=lambda, $
       spec=spec,ivar=ivar, zans=zans, synspec=synspec, fulldir=fulldir
   SELECT_Z, zall, objnames[slitcount], zres=zres

   redshift = zans.z
   oldspec = spec
   oldsynspec = synspec
   IF(nsmooth GT 1) THEN BEGIN
     spec = IVARSMOOTH(spec,ivar,nsmooth)
     synspec = IVARSMOOTH(synspec,ivar,nsmooth)
   ENDIF
   MY_PLOT_SPEC
   MY_SHOW_TEXT, zans, zres
   postage_stamp, maskname, slitnames[slitcount], redshift

red_str = 'Current settings:  Q = ' + $
    STRING(fin_quality[slitcount], FORMAT='(I1)') + $
    '  z = ' + STRING(fin_redshifts[slitcount], FORMAT='(F5.3)') + $
    '   ' + fin_canned[slitcount] + '  ' + fin_comments[slitcount]

WIDGET_CONTROL, finprevious_data, SET_VALUE=red_str
WIDGET_CONTROL, redshift_list, SET_DROPLIST_SELECT=0
WIDGET_CONTROL, canned_list, SET_DROPLIST_SELECT=0

 END
;********************************************************************




;*************************  BACK BUTTON  *******************************
; Plot the previous spectrum
;*****************
 'BACK' : BEGIN 
    slitcount = slitcount - 1    ; Move to the previous slit
;   PRINT, slitcount

    IF slitcount LT 0 THEN BEGIN
      slitcount = 1
      GOTO, ERROR2
    ENDIF

   ; Check if we're only looking at certain Q's
   IF zqq EQ 1 THEN BEGIN
    indx_cnt = indx_cnt - 1
    IF indx_cnt LT 0 THEN BEGIN
      indx_cnt = 0
      GOTO, ERROR2
    ENDIF
    slitcount = indx_pos[indx_cnt]
   ENDIF






    DEIMOS_READOBJ, maskname, objnames[slitcount], lambda=lambda, $
   spec=spec,ivar=ivar, zans=zans, synspec=synspec, fulldir=fulldir
   SELECT_Z, zall, objnames[slitcount], zres=zres

   redshift = zans.z
   oldspec = spec
   oldsynspec = synspec
   IF(nsmooth GT 1) THEN BEGIN
     spec = IVARSMOOTH(spec,ivar,nsmooth)
     synspec = IVARSMOOTH(synspec,ivar,nsmooth)
   ENDIF
   MY_PLOT_SPEC
   MY_SHOW_TEXT, zans, zres
   postage_stamp, maskname, slitnames[slitcount], redshift

red_str = 'Current settings:  Q = ' + $
    STRING(fin_quality[slitcount], FORMAT='(I1)') + $
    '  z = ' + STRING(fin_redshifts[slitcount], FORMAT='(F5.3)') + $
    '   ' + fin_canned[slitcount] + '  ' + fin_comments[slitcount]

WIDGET_CONTROL, finprevious_data, SET_VALUE=red_str
WIDGET_CONTROL, redshift_list, SET_DROPLIST_SELECT=0
WIDGET_CONTROL, canned_list, SET_DROPLIST_SELECT=0

 END
;********************************************************************




; If done is selected the session is ended
;*******************************************
 'DONE': BEGIN
    WIDGET_CONTROL, ev.top, /DESTROY
  END
 else: junk = 0 ;junk statement
 ENDCASE

;*******************************************
;            ERROR HANDLING
;*******************************************

  IF 1 EQ 0 THEN BEGIN
    ERROR1: PRINT, 'ERROR: maximum nslit exceeded!'
    MY_PRINT_ERROR, 'Thats all folks!' 
  ENDIF


  IF 1 EQ 0 THEN BEGIN
    ERROR2: PRINT, 'ERROR: selected nslit is less than 0!'
    MY_PRINT_ERROR, 'ERROR: selected nslit is less than 0!!!'
  ENDIF


END












PRO pspec2_selectz, event
;  Change the selected redshift
;********************************

COMMON SPEC_STUFF, maskname, slitnames, slitcount, nslit, objnames, zall
COMMON SPECTRUM,  lambda, spec, ivar, synspec, redshift, oldspec, oldsynspec, nsmooth

button_index = event.index

   SELECT_Z, zall, objnames[slitcount], zres=zres, synspec=synspec, $
              index=button_index, lambda=lambda
   redshifts = zres.z
   redshift = redshifts[button_index]

   oldsynspec = synspec
   IF(nsmooth GT 1) THEN $
     synspec = IVARSMOOTH(synspec,ivar,nsmooth)
   MY_PLOT_SPEC
   MY_SHOW_TEXT, zres[button_index], zres
   postage_stamp, maskname, slitnames[slitcount], redshift
END










PRO pspec2_smooth, event
;  Change the selected smoothing
;*******************************

COMMON SPECTRUM,  lambda, spec, ivar, synspec, redshift, oldspec, oldsynspec, nsmooth

button_index = event.index

IF button_index EQ 0 THEN BEGIN
   spec = oldspec
   synspec = oldsynspec
   nsmooth = 1
   MY_PLOT_SPEC
ENDIF ELSE BEGIN
 kernel = [3,5,7,11,13]
 kk = kernel[button_index-1]

 spec    = ivarsmooth(oldspec,ivar,kk)
 synspec = ivarsmooth(oldsynspec,ivar,kk)
 nsmooth = kk
 MY_PLOT_SPEC
ENDELSE

END









PRO pspec2_jump, event
; Jump to a new slit number
;*******************************

COMMON SPEC_STUFF, maskname, slitnames, slitcount, nslit, objnames, zall
COMMON SPECTRUM,  lambda, spec, ivar, synspec, redshift, oldspec, oldsynspec, nsmooth
COMMON BUTTON_LABELS, redshift_list, canned_list ; Must update the list names!
COMMON RESULTS, fin_redshifts, fin_quality, fin_comments, fin_canned, finprevious_data


WIDGET_CONTROL, event.id, GET_VALUE = slitnum

cnt = WHERE(FIX(slitnames) EQ FIX(slitnum[0]), c)
IF c GE 1 THEN BEGIN
  slitcount = cnt[0]
   objn = objnames[slitcount]

   DEIMOS_READOBJ, maskname, objn[0], lambda=lambda, $
   spec=spec,ivar=ivar, zans=zans, synspec=synspec, fulldir=fulldir
   SELECT_Z, zall, objn[0], zres=zres

   redshift = zans.z
   oldspec = spec
   oldsynspec = synspec
   IF(nsmooth GT 1) THEN BEGIN
     spec = ivarsmooth(spec,ivar,nsmooth)
     synspec = ivarsmooth(synspec,ivar,nsmooth)
   ENDIF
   MY_PLOT_SPEC
   MY_SHOW_TEXT, zans, zres
   postage_stamp, maskname, slitnames[slitcount], redshift
ENDIF
IF c LE 0 THEN BEGIN
  MY_PRINT_ERROR, 'Invalid slit selection...'
ENDIF

WIDGET_CONTROL, event.id, SET_VALUE = ''


red_str = 'Current settings:  Q = ' + $
    STRING(fin_quality[slitcount], FORMAT='(I1)') + $
    '  z = ' + STRING(fin_redshifts[slitcount], FORMAT='(F5.3)') + $
    '   ' + fin_canned[slitcount] + '  ' + fin_comments[slitcount]

WIDGET_CONTROL, finprevious_data, SET_VALUE=red_str
WIDGET_CONTROL, redshift_list, SET_DROPLIST_SELECT=0
WIDGET_CONTROL, canned_list, SET_DROPLIST_SELECT=0

END









FUNCTION pspec2_setz, event
;  Set the redshift quality flag for this objects
;**********************************************

COMMON SPEC_STUFF, maskname, slitnames, slitcount, nslit, objnames, zall
COMMON RESULTS, fin_redshifts, fin_quality, fin_comments, fin_canned, finprevious_data
COMMON SPECTRUM,  lambda, spec, ivar, synspec, redshift, oldspec, oldsynspec, nsmooth

qq = [-1,1,2,3,4]

    fin_redshifts[slitcount] = redshift
    fin_quality[slitcount] = qq[event.value]

red_str = 'Current settings:  Q = ' + $
    STRING(fin_quality[slitcount], FORMAT='(I1)') + $
    '  z = ' + STRING(fin_redshifts[slitcount], FORMAT='(F5.3)') + $
    '   ' + fin_canned[slitcount] + '  ' + fin_comments[slitcount]

WIDGET_CONTROL, finprevious_data, SET_VALUE=red_str

END







PRO pspec2_comment, event
;  Store a comment about a spectrum
;**************************************

COMMON RESULTS, fin_redshifts, fin_quality, fin_comments, fin_canned, finprevious_data
COMMON SPEC_STUFF, maskname, slitnames, slitcount, nslit, objnames, zall

WIDGET_CONTROL, event.id, GET_VALUE = comment
fin_comments[slitcount] = comment
WIDGET_CONTROL, event.id, SET_VALUE = ''


red_str = 'Current settings:  Q = ' + $
    STRING(fin_quality[slitcount], FORMAT='(I1)') + $
    '  z = ' + STRING(fin_redshifts[slitcount], FORMAT='(F5.3)') + $
    '   ' + fin_canned[slitcount] + '  ' + fin_comments[slitcount]

WIDGET_CONTROL, finprevious_data, SET_VALUE=red_str


END







PRO pspec2_canned, event
;  Store a 'canned' comment as a single character
;************************************************

COMMON RESULTS, fin_redshifts, fin_quality, fin_comments, fin_canned, finprevious_data
COMMON SPEC_STUFF, maskname, slitnames, slitcount, nslit, objnames, zall

comment_index = event.index

; If 'none' is selected, wipe all present
; otherwise append comment character
;***********************************
IF comment_index EQ 0 THEN BEGIN
  fin_canned[slitcount] = '-'
ENDIF ELSE BEGIN
  comments = ['o','a','e','r']
  IF fin_canned[slitcount] EQ '-' THEN $
    fin_canned[slitcount] = comments[comment_index-1] $
  ELSE $
   fin_canned[slitcount] = fin_canned[slitcount] + comments[comment_index-1]
ENDELSE

red_str = 'Current settings:  Q = ' + $
    STRING(fin_quality[slitcount], FORMAT='(I1)') + $
    '  z = ' + STRING(fin_redshifts[slitcount], FORMAT='(F5.3)') + $
    '   ' + fin_canned[slitcount] + '  ' + fin_comments[slitcount]

WIDGET_CONTROL, finprevious_data, SET_VALUE=red_str


END





PRO pspec2_atv, event
; Bring up an atv window to view 2d spectrum
;*****************************************

COMMON SPEC_STUFF, maskname, slitnames, slitcount, nslit, objnames, zall
common pspec2_dirs, topdir0, fulldir



   if not(keyword_set(fulldir)) then begin
       if (topdir0 NE '.') then begin
           filelist=findfile(getenv('D2_RESULTS') + '/' + $
                             maskname + '/*/obj_info*.fits', count=nobjfile)
           if nobjfile eq 0 then $
             message, 'No obj_info.' + maskname + '.fits file found!!'
           if nobjfile gt 1 then begin
               print, 'More than one reduction for mask: '
               print, filelist
               read, 'Enter file number you want (starting with 0)', num
               filelist = filelist[num]
           endif else filelist=filelist[0]
           dirlist = strsplit(filelist, '/', /extract)
           ndir = n_elements(dirlist)
           date = dirlist[ndir-2]
           fulldir = concat_dir((topdir+'/'+maskname),date)
       endif       
   endif



  fname1 = 'slit.'+strcompress(maskname, /remove_all)+'.'+strcompress(slitnames[slitcount], /remove_all)+'B.fits.gz'
  fname1 = concat_dir(fulldir,fname1)

  fname2 = 'slit.'+strcompress(maskname, /remove_all)+'.'+strcompress(slitnames[slitcount], /remove_all)+'R.fits.gz'
  fname2 = concat_dir(fulldir,fname2)
 
  spec_blue = MRDFITS(fname1,1,/silent,status=status1)
;  lam_blue  = LAMBDA_EVAL(spec_blue)
  spec_red  = MRDFITS(fname2,1,/silent,status=status2)
;  lam_red   = LAMBDA_EVAL(spec_red)
IF status1 NE -1 THEN BEGIN
  nbins = N_ELEMENTS(spec_blue.lambda0)/4
ENDIF ELSE BEGIN
  nbins = N_ELEMENTS(spec_red.lambda0)/4
ENDELSE
  ; Divide spectrum into four parts
IF status1 NE -1 THEN BEGIN
  flux1 = spec_blue.flux[0:nbins-1,*]
  flux2 = spec_blue.flux[nbins:2*nbins-1,*]
  flux3 = spec_blue.flux[2*nbins:3*nbins-1,*]
  flux4 = spec_blue.flux[3*nbins:4*nbins-1,*]
  lam1  = spec_blue.lambda0[0:nbins-1]
  lam2  = spec_blue.lambda0[nbins:2*nbins-1]
  lam3  = spec_blue.lambda0[2*nbins:3*nbins-1]
  lam4  = spec_blue.lambda0[3*nbins:4*nbins-1]
  flux_blue = [[flux4],[lam4],[lam4],[flux3],[lam3],[lam3],[flux2],[lam2],[lam2],[flux1],[lam1],[lam1]]
ENDIF ELSE BEGIN
  flux_blue = REPLICATE(0,50,nbins)
ENDELSE

IF status2 NE -1 THEN BEGIN
  flux1 = spec_red.flux[0:nbins-1,*]
  flux2 = spec_red.flux[nbins:2*nbins-1,*]
  flux3 = spec_red.flux[2*nbins:3*nbins-1,*]
  flux4 = spec_red.flux[3*nbins:4*nbins-1,*]
  lam1  = spec_red.lambda0[0:nbins-1]
  lam2  = spec_red.lambda0[nbins:2*nbins-1]
  lam3  = spec_red.lambda0[2*nbins:3*nbins-1]
  lam4  = spec_red.lambda0[3*nbins:4*nbins-1]
  flux_red = [[flux4],[lam4],[lam4],[flux3],[lam3],[lam3],[flux2],[lam2],[lam2],[flux1],[lam1],[lam1]]
ENDIF ELSE BEGIN
  flux_red = REPLICATE(0,50,nbins)
ENDELSE

;  min_b = MIN(spec_blue.flux)
;  min_r = MIN(spec_red.flux)
;  max_b = MAX(spec_blue.flux)
;  max_r = MAX(spec_red.flux)

IF status1 NE -1 AND status2 NE -1 THEN $
  atv, [[flux_red],[flux_blue]], min=-5, max=25
IF status1 EQ -1 THEN $
  atv, flux_red, min=-5,max=25
IF status2 EQ -1 THEN $
  atv, flux_blue, min=-5,max=25
 ;min=MIN([min_b,min_r]), max=MAX([max_b,max_r])

END








;***********
;  2D Plotting: A.Coil
;
; plots the OII doublet from the 2D spectrum for a given mask and slitn
; have to be in the reduction directory to execute this

pro postage_stamp,  mask,  slitn,  redshift
; ex: postage_stamp, 3206, '141', 1.0, 3727 (for OII)
; linewave must equal 1,2,3 or 4

COMMON DRAW_STUFF_2D, drawID5, drawID6, drawID7, drawID8
common pspec2_dirs, topdir0, fulldir

;tvscl, REPLICATE(0.0,200,100)



;  topdir = getenv('D2_RESULTS')
;  if topdir eq '' then print, 'You should set $D2_RESULTS'
;  topdir = concat_dir(topdir,'v0_9')
   if not(keyword_set(fulldir)) then begin
       topdir = topdir0
       if (topdir NE '.') then begin
           filelist=findfile(getenv('D2_RESULTS') + '/' + $
                             mask + '/*/obj_info*.fits', count=nobjfile) 
           if nobjfile eq 0 then $
             message, 'No obj_info.' + mask + '.fits file found!!'
           if nobjfile gt 1 then begin
               print, 'More than one reduction for mask: '
               print, filelist
               read, 'Enter file number you want (starting with 0)', num
               filelist = filelist[num]
           endif else filelist=filelist[0]
           dirlist = strsplit(filelist, '/', /extract)
           ndir = n_elements(dirlist)
           date = dirlist[ndir-2]
           fulldir = concat_dir((topdir+'/'+mask),date)
       endif                  
   endif


  fname1 = 'slit.'+strcompress(mask, /remove_all)+'.'+strcompress(slitn, /remove_all)+'B.fits.gz'
  fname1 = concat_dir(fulldir,fname1)  

  fname2 = 'slit.'+strcompress(mask, /remove_all)+'.'+strcompress(slitn, /remove_all)+'R.fits.gz'
  fname2 = concat_dir(fulldir,fname2)  

; read in the 2d spectra, blue and red
slit_blue = mrdfits(fname1, 1, /silent,status=status1)
slit_red  = mrdfits(fname2, 1, /silent,status=status2)
IF (status1 EQ -1) THEN BEGIN
  PRINT, fname1, ' not found, check ATV for 2d images'
  RETURN
ENDIF
IF (status2 EQ -1) THEN BEGIN
  PRINT, fname2, ' not found, check ATV for 2d images'
  RETURN
ENDIF

; get the 2d wavelength solution from slit.lambdax and slit.dlam
lambda_blue=lambda_eval(slit_blue)
lambda_red=lambda_eval(slit_red)







; now find where the object is spatially - in spec1d file
spec1dfile = findfile(concat_dir(fulldir,'spec1d.'+strcompress(mask, /remove_all)+'.'+strcompress(slitn, /remove_all)+'.*.fits'), count = ct)


; if ct is greater than 1 then there's more than one object on this
; slit?
spec1d = mrdfits(spec1dfile[0], 1, /silent)

; pull out the object plus a little region around it - use 5*FWHM
;spatially = findgen(fix(spec1d.fwhm*5))+round(spec1d.objpos)-2.5*round(spec1d.fwhm)
spatially = FINDGEN(50) + ROUND(spec1d.objpos) - 25




;***++++++++++++++++++++++++++++++
;  Check to see if each line is present.  If it is
;  then plot it.

; [OII]
;*******
 WSET, drawID5

 tv, bytarr(200, 100)
 lambda_min  = (1.+redshift)*3707
 lambda_max = (1.+redshift)*3747

; find region around OII 
oii_lambda_blue = where(lambda_blue gt lambda_min and lambda_blue lt lambda_max,nb)
oii_lambda_red = where(lambda_red gt lambda_min and lambda_red lt lambda_max,nr)

; figure out which side, red or blue, is the one with OII
 blue = 0
 red = 0
 if nb gt 0 then blue = 1
 if nr gt 0  then red = 1


if blue eq 1 then begin
  nrows = (size(lambda_blue))[2]
  ncols = n_elements(oii_lambda_blue)/float(nrows)
  oii_lambda_blue = oii_lambda_blue[0:ncols-1]
  flux = slit_blue.flux[oii_lambda_blue, *]
endif
if red eq 1 then begin
  nrows = (size(lambda_red))[2]
  ncols = n_elements(oii_lambda_red)/float(nrows)
  oii_lambda_red = oii_lambda_red[0:ncols-1]
  flux = slit_red.flux[oii_lambda_red, *]
endif
if red eq 0 and blue eq 0 then begin
  GOTO, NEXT
endif


flux = flux[*,spatially]
nlam1 = FINDGEN(200) * N_ELEMENTS(flux[*,0])/200
nlam2 = FINDGEN(100)
flux = INTERPOLATE(flux,nlam1,nlam2,/grid)

;IF MAX(flux) GT 0 THEN BEGIN
; flux = HIST_EQUAL(flux)
; flux = SMOOTH(flux,3)
;ENDIF


fullrange = minmax(flux)

if fullrange[1] lt (-5) then minscl = fullrange[0] else minscl =  (-5)
if fullrange[0] gt (25) then maxscl = fullrange[1] else maxscl =  (25)

tv, bytscl(flux, min = minscl, max = maxscl, top = 255b-8b) +8b

NEXT:



; H+K
;*******
 WSET, drawID6

 tv, bytarr(200, 100)

 lambda_min  = (1.+redshift)*3900
 lambda_max = (1.+redshift)*4000

; find region around OII 
oii_lambda_blue = where(lambda_blue gt lambda_min and lambda_blue lt lambda_max,nb)
oii_lambda_red = where(lambda_red gt lambda_min and lambda_red lt lambda_max,nr)

; figure out which side, red or blue, is the one with OII
 blue = 0
 red = 0
 if nb gt 0 then blue = 1
 if nr gt 0  then red = 1


if blue eq 1 then begin
  nrows = (size(lambda_blue))[2]
  ncols = n_elements(oii_lambda_blue)/float(nrows)
  oii_lambda_blue = oii_lambda_blue[0:ncols-1]
  flux = slit_blue.flux[oii_lambda_blue, *]
endif
if red eq 1 then begin
  nrows = (size(lambda_red))[2]
  ncols = n_elements(oii_lambda_red)/float(nrows)
  oii_lambda_red = oii_lambda_red[0:ncols-1]
  flux = slit_red.flux[oii_lambda_red, *]
endif
if red eq 0 and blue eq 0 then begin
  GOTO, NEXT2
endif


flux = flux[*,spatially]
nlam1 = FINDGEN(200) * N_ELEMENTS(flux[*,0])/200
nlam2 = FINDGEN(100)
flux = INTERPOLATE(flux,nlam1,nlam2,/grid)
fullrange = minmax(flux)


if fullrange[1] lt (-5) then minscl = fullrange[0] else minscl =  (-5)
if fullrange[0] gt (25) then maxscl = fullrange[1] else maxscl =  (25)

tv, bytscl(flux, min = minscl, max = maxscl, top = 255b-8b) +8b

NEXT2:



IF redshift GT 0.3 THEN BEGIN

; Hbeta
;*******
 WSET, drawID7

 tv, bytarr(200, 100)
 lambda_min  = (1.+redshift)*4840
 lambda_max = (1.+redshift)*4880

; find region around OII 
oii_lambda_blue = where(lambda_blue gt lambda_min and lambda_blue lt lambda_max,nb)
oii_lambda_red = where(lambda_red gt lambda_min and lambda_red lt lambda_max,nr)

; figure out which side, red or blue, is the one with OII
 blue = 0
 red = 0
 if nb gt 0 then blue = 1
 if nr gt 0  then red = 1


if blue eq 1 then begin
  nrows = (size(lambda_blue))[2]
  ncols = n_elements(oii_lambda_blue)/float(nrows)
  oii_lambda_blue = oii_lambda_blue[0:ncols-1]
  flux = slit_blue.flux[oii_lambda_blue, *]
endif
if red eq 1 then begin
  nrows = (size(lambda_red))[2]
  ncols = n_elements(oii_lambda_red)/float(nrows)
  oii_lambda_red = oii_lambda_red[0:ncols-1]
  flux = slit_red.flux[oii_lambda_red, *]
endif
if red eq 0 and blue eq 0 then begin
  GOTO, NEXT3
endif

flux = flux[*, spatially]
nlam1 = FINDGEN(200) * N_ELEMENTS(flux[*,0])/200
nlam2 = FINDGEN(100)
flux = INTERPOLATE(flux,nlam1,nlam2,/grid)


;IF MAX(flux) GT 0 THEN BEGIN
; flux = HIST_EQUAL(flux)
; flux = SMOOTH(flux,3)
;ENDIF

fullrange = minmax(flux)


if fullrange[1] lt (-5) then minscl = fullrange[0] else minscl =  (-5)
if fullrange[0] gt (25) then maxscl = fullrange[1] else maxscl =  (25)

tv, bytscl(flux, min = minscl, max = maxscl, top = 255b-8b) +8b

NEXT3:




; [OIII]
;*******
 WSET, drawID8

 tv, bytarr(200, 100)
 lambda_min  = (1.+redshift)*4940
 lambda_max = (1.+redshift)*5030

; find region around OII 
oii_lambda_blue = where(lambda_blue gt lambda_min and lambda_blue lt lambda_max,nb)
oii_lambda_red = where(lambda_red gt lambda_min and lambda_red lt lambda_max,nr)

; figure out which side, red or blue, is the one with OII
 blue = 0
 red = 0
 if nb gt 0 then blue = 1
 if nr gt 0  then red = 1


if blue eq 1 then begin
  nrows = (size(lambda_blue))[2]
  ncols = n_elements(oii_lambda_blue)/float(nrows)
  oii_lambda_blue = oii_lambda_blue[0:ncols-1]
  flux = slit_blue.flux[oii_lambda_blue, *]
endif
if red eq 1 then begin
  nrows = (size(lambda_red))[2]
  ncols = n_elements(oii_lambda_red)/float(nrows)
  oii_lambda_red = oii_lambda_red[0:ncols-1]
  flux = slit_red.flux[oii_lambda_red, *]
endif
if red eq 0 and blue eq 0 then begin
  RETURN
endif

flux = flux[*, spatially]
nlam1 = FINDGEN(200) * N_ELEMENTS(flux[*,0])/200
nlam2 = FINDGEN(100)
flux = INTERPOLATE(flux,nlam1,nlam2,/grid)

;IF MAX(flux) GT 0 THEN BEGIN
; flux = HIST_EQUAL(flux)
; flux = SMOOTH(flux,3)
;ENDIF

fullrange = minmax(flux)

if fullrange[1] lt (-5) then minscl = fullrange[0] else minscl =  (-5)
if fullrange[0] gt (25) then maxscl = fullrange[1] else maxscl =  (25)

tv, bytscl(flux, min = minscl, max = maxscl, top = 255b-8b) +8b






ENDIF  ELSE BEGIN
;
; Low z stuff
;************

 WSET, drawID7

 tv, bytarr(200, 100)
 lambda_min  = (1.+redshift)*6540
 lambda_max = (1.+redshift)*6600

; find region around Halpha
oii_lambda_blue = where(lambda_blue gt lambda_min and lambda_blue lt lambda_max,nb)
oii_lambda_red = where(lambda_red gt lambda_min and lambda_red lt lambda_max,nr)

; figure out which side, red or blue, is the one with OII
 blue = 0
 red = 0
 if nb gt 0 then blue = 1
 if nr gt 0  then red = 1


if blue eq 1 then begin
  nrows = (size(lambda_blue))[2]
  ncols = n_elements(oii_lambda_blue)/float(nrows)
  oii_lambda_blue = oii_lambda_blue[0:ncols-1]
  flux = slit_blue.flux[oii_lambda_blue, *]
endif
if red eq 1 then begin
  nrows = (size(lambda_red))[2]
  ncols = n_elements(oii_lambda_red)/float(nrows)
  oii_lambda_red = oii_lambda_red[0:ncols-1]
  flux = slit_red.flux[oii_lambda_red, *]
endif
if red eq 0 and blue eq 0 then begin
  GOTO, NEXT8
endif

flux = flux[*, spatially]
nlam1 = FINDGEN(200) * N_ELEMENTS(flux[*,0])/200
nlam2 = FINDGEN(100)
flux = INTERPOLATE(flux,nlam1,nlam2,/grid)


;IF MAX(flux) GT 0 THEN BEGIN
; flux = HIST_EQUAL(flux)
; flux = SMOOTH(flux,3)
;ENDIF

fullrange = minmax(flux)


if fullrange[1] lt (-5) then minscl = fullrange[0] else minscl =  (-5)
if fullrange[0] gt (25) then maxscl = fullrange[1] else maxscl =  (25)

tv, bytscl(flux, min = minscl, max = maxscl, top = 255b-8b) +8b

NEXT8:




; [SII]
;*******
 WSET, drawID8

 tv, bytarr(200, 100)
 lambda_min  = (1.+redshift)*6700
 lambda_max = (1.+redshift)*6745

; find region around OII 
oii_lambda_blue = where(lambda_blue gt lambda_min and lambda_blue lt lambda_max,nb)
oii_lambda_red = where(lambda_red gt lambda_min and lambda_red lt lambda_max,nr)

; figure out which side, red or blue, is the one with OII
 blue = 0
 red = 0
 if nb gt 0 then blue = 1
 if nr gt 0  then red = 1


if blue eq 1 then begin
  nrows = (size(lambda_blue))[2]
  ncols = n_elements(oii_lambda_blue)/float(nrows)
  oii_lambda_blue = oii_lambda_blue[0:ncols-1]
  flux = slit_blue.flux[oii_lambda_blue, *]
endif
if red eq 1 then begin
  nrows = (size(lambda_red))[2]
  ncols = n_elements(oii_lambda_red)/float(nrows)
  oii_lambda_red = oii_lambda_red[0:ncols-1]
  flux = slit_red.flux[oii_lambda_red, *]
endif
if red eq 0 and blue eq 0 then begin
  RETURN
endif

flux = flux[*, spatially]
nlam1 = FINDGEN(200) * N_ELEMENTS(flux[*,0])/200
nlam2 = FINDGEN(100)
flux = INTERPOLATE(flux,nlam1,nlam2,/grid)

;IF MAX(flux) GT 0 THEN BEGIN
; flux = HIST_EQUAL(flux)
; flux = SMOOTH(flux,3)
;ENDIF

fullrange = minmax(flux)

if fullrange[1] lt (-5) then minscl = fullrange[0] else minscl =  (-5)
if fullrange[0] gt (25) then maxscl = fullrange[1] else maxscl =  (25)

tv, bytscl(flux, min = minscl, max = maxscl, top = 255b-8b) +8b



ENDELSE


END







;*************************** MAIN ROUTINE ***********************

PRO pspec2, mname, curdir=curdir,  noremove=noremove, oldfile=oldfile, zquality=zquality
; 
; WIDGET initialisation.
;
; -  mname is the mask name e.g. 4206
;
;**********************************************


!QUIET = 1
COMMON SPEC_STUFF, maskname, slitnames, slitcount, nslit, objnames, zall
COMMON DRAW_STUFF, drawID1, drawID2, drawID3, drawID4
COMMON DRAW_STUFF_2D, drawID5, drawID6, drawID7, drawID8
COMMON TEXT_STUFF, text_widg, text_widgz, pgal
COMMON SPECTRUM,  lambda, spec, ivar, synspec, redshift, oldspec, oldsynspec, nsmooth
COMMON RESULTS, fin_redshifts, fin_quality, fin_comments, fin_canned, finprevious_data
common pspec2_dirs, topdir0, fulldir
COMMON ZJUNK, zresfile
COMMON BUTTON_LABELS, redshift_list, canned_list
zresfile = -1
nsmooth = 1
COMMON SKIP_STUFF, indx_pos, zqq, indx_cnt
indx_pos = 0
indx_cnt = 0
zqq = 0   ; Keeps track of the ZQUALITY keyword


 ; Default to mask 3206 if no other selection is made.
 ;+++++++++++++++++++++++++++++++++++++++++
 IF (N_PARAMS() LT 1) THEN BEGIN
   PRINT, 'Should supply maskname on command line!'
   PRINT, 'Assuming mask is 3206 instead...'
   maskname = '3206'
 ENDIF ELSE BEGIN
   maskname = STRING(mname,FORMAT='(I4.4)')
 ENDELSE





slitcount = 0   ; Keeps track of what slit we are looking at.



 ; Use zall.slitname to get slitnumbers
 ;++++++++++++++++++++++++++++++++++++++++++++++
   zall = READ_Z(maskname)

   if keyword_set(curdir) then fulldir='.' $
   else begin
       topdir = getenv('D2_RESULTS') ;either local or global
       filelist=findfile(topdir + '/' + $
                         maskname + '/*/obj_info*.fits', count=nobjfile) 
       if nobjfile eq 0 then $
         message, 'No obj_info.' + maskname + '.fits file found!!'
       if nobjfile gt 1 then begin
           print, 'More than one reduction for mask: '
           for ii=0,nobjfile-1 do begin
               pos = strpos(filelist[ii], 'obj_info')
               print, strmid(filelist[ii], 0, pos)
           endfor
           read, 'Enter file number you want (starting with 0)', num
           if num eq -1 then begin
               fulldir = ''
               read, 'Enter full path of data directory:', fulldir
               filelist = findfile(fulldir + '/obj_info*.fits', count=nobjfile)
               if nobjfile gt 0 then num = 0 $
               else message, 'No obj_info.' + maskname + '.fits file found!!'
           endif
           filelist = filelist[num]
       endif else filelist=filelist[0]
       dirlist = strsplit(filelist, '/', /extract)
       ndir = n_elements(dirlist)
       date = dirlist[ndir-2]
       if not(keyword_set(fulldir)) then $
              fulldir = concat_dir((topdir+'/'+maskname),date)
       topdir0 = topdir
   endelse


   DEIMOS_READSPEC, maskname, '001', zbest=zbest, fulldir=fulldir

   slitlist = zbest.slitname
   slitnames = STRING(slitlist, FORMAT='(I3.3)')
   nslit = N_ELEMENTS(zbest.slitname)

   objlist  = zbest.objname

   ; Change to object names (remove serendips)
   serendip_count = 0
   FOR i=0, nslit-1 DO BEGIN
     IF(STRMID(objlist[i],0,1) EQ 's') THEN $
      serendip_count = serendip_count + 1
   ENDFOR
   PRINT, 'No of serendips: ', serendip_count

   serendip = 0
   objnames = STRARR(nslit - serendip_count)
   slnames  = STRARR(nslit - serendip_count)
   FOR i=0, nslit-1 DO BEGIN
     IF(STRMID(objlist[i],0,1) NE 's') THEN BEGIN
      objnames[serendip] = objlist[i]
      slnames[serendip] = slitnames[i]
      serendip = serendip + 1
     ENDIF
   ENDFOR

   slitnames = slnames
   nslit = N_ELEMENTS(objnames)


   ; Remove allignment stars and sky slits (only works for v1_0)
 IF NOT(keyword_set(noremove)) THEN BEGIN
   info_file =  findfile(topdir+'/'+maskname+'/*/obj_info*.fits')
;   PRINT,  info_file[0]
   
   counter = 0
   aln_stars =  0
   sky_stuff =  0
   objnames_new =  STRARR(nslit)
   slitnames_new =  STRARR(nslit)
   info =  MRDFITS(info_file[0], 1, /silent)
   FOR i=0,  nslit-1 DO BEGIN
     cnt = WHERE(info.objno EQ objnames[i], c)
     IF c NE 2 THEN PRINT, 'Error: wrong number of entries in infofile for obj: ',  objnames[i]
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
 ;++++++++++++++++++++++++++++++++++++++++++++++++++++++

;TBD: read_z


  print, 'Using slitfiles and spec1dfiles from ' + fulldir







 ; Redshift flags and comments:
 ;*****************************

  fin_redshifts = FLTARR(nslit)
  fin_quality   = INTARR(nslit)
  fin_comments  = STRARR(nslit)
  fin_canned    = MAKE_ARRAY(nslit,/STRING,VALUE='-')   ; 'Canned' comments



 ; If oldfile is set then read in the previous values of z etc
 ;*************************************************************

 IF (keyword_set(oldfile)) THEN BEGIN
 PRINT,  'Reading in previous results from:', oldfile
  OPENR, 1, oldfile
  buffer = ''
  READF, 1, buffer
  READF, 1, buffer
  objn = 0L
  slitn = 0
  zz =  0.0
  q =  0
  can = ''
  com = ''
  WHILE NOT EOF(1) DO BEGIN
   READF, 1, buffer
   IF buffer EQ '' THEN CONTINUE  ; Needed if blank line at end of file
   READS, buffer, objn, slitn, zz, q, can, com, FORMAT='(I9,I6,F9.5,I4,A7,A0)'
;   PRINT, objn, slitn, zz, q, can, com
   cnt =  WHERE(objnames EQ objn, c)
   IF c EQ 1 THEN BEGIN
    fin_redshifts[cnt] =  zz
    fin_quality[cnt] =  q
    fin_canned[cnt] =  STRCOMPRESS(can, /REMOVE_ALL)
    fin_comments[cnt] =  STRTRIM(com, 1)
   ENDIF ELSE BEGIN
    PRINT, 'No match for ', objn, ' from file:  ', oldfile
   ENDELSE
  ENDWHILE
  CLOSE, 1
 ENDIF





 ; Grab pgal for each object
 ;***************************
 binfile = findfile(topdir+'/'+maskname+'/*/*bintab*.fits*', count=nbin)
 if nbin eq 0 then $
        message, 'ERROR: NO bin table file found!'
 binfile = binfile[0]
 fits_open, binfile, fcb
 ; find the pcat extension in the bin table file and then
 ; read it.
 pdex = where(fcb.extname eq 'PhotCat', pcnt)
 if pcnt eq 0 then $
         print, 'No pcat info found for mask!' $
 else pcat = mrdfits(binfile, pdex[0], /silent)

 pgal = FLTARR(nslit)
 FOR i=0, nslit-1 DO BEGIN
   cnt = WHERE(pcat.objno EQ LONG(objnames[i]), c)
   IF c EQ 1 THEN BEGIN
     pgal[i] = pcat[cnt].pgal
   ENDIF ELSE $
     PRINT, 'ERROR: No pgal for object ', objnames[i]
 ENDFOR




 ; If Q is set then only look at objects with the specified Q value
 ;****************************************************************

 nslit_tot = nslit
 IF (keyword_set(zquality)) THEN BEGIN
 PRINT, 'Extracting only spectra with zquality = ', zquality
  cnt = WHERE(fin_quality EQ zquality, c)

  IF c GT 0 THEN BEGIN
    PRINT, c, ' objects to be reviewed.'
    indx_pos = cnt
    nslit =  c
    zqq = 1
  ENDIF ELSE BEGIN
    PRINT, 'ERROR: ', c, ' objects found!'
    PRINT, 'Ignoring the ZQUALITY flag'
  ENDELSE
 ENDIF







;********************************
;     WIDGET STUFF
;********************************

; Main widget thing
base = WIDGET_BASE(/COLUMN, TITLE='pspec2')

; Graphical tool
draw_base1 = WIDGET_BASE(base,/ROW)
draw1 = WIDGET_DRAW(draw_base1, XSIZE=200, YSIZE=200)
draw2 = WIDGET_DRAW(draw_base1, XSIZE=200, YSIZE=200)
draw3 = WIDGET_DRAW(draw_base1, XSIZE=200, YSIZE=200)
draw4 = WIDGET_DRAW(draw_base1, XSIZE=200, YSIZE=200)
draw_base2 = WIDGET_BASE(base,/ROW)
draw5 = WIDGET_DRAW(draw_base2, XSIZE=200, YSIZE=50)
draw6 = WIDGET_DRAW(draw_base2, XSIZE=200, YSIZE=50)
draw7 = WIDGET_DRAW(draw_base2, XSIZE=200, YSIZE=50)
draw8 = WIDGET_DRAW(draw_base2, XSIZE=200, YSIZE=50)



basex = WIDGET_BASE(base, /ROW)

; Main text field
text1 = WIDGET_TEXT(basex, XSIZE=55, YSIZE=8, /SCROLL)


basex2 = WIDGET_BASE(basex, /COLUMN)


; Redshifts
;**********
base4 = WIDGET_BASE(basex2, /ROW, /FRAME)
values = ['Best', '2nd', '3rd', '4th','5th','6th','7th','8th','9th','10th']
redshift_list = WIDGET_DROPLIST(base4, VALUE=values, $ 
   TITLE='Redshift selection:', $ 
   EVENT_PRO='pspec2_selectz')
text2 = WIDGET_TEXT(base4, XSIZE=35, YSIZE=1, /scroll)

; Smoothing
;***********
base5 = WIDGET_BASE(basex2, /ROW)
values = ['None', '3 pix', '5 pix', '7 pix', '11 pix', '13 pix']
smooth_list = WIDGET_DROPLIST(base5, VALUE=values, $
    TITLE='Smoothing:', /FRAME, EVENT_PRO='pspec2_smooth')

; ATV
;******
atv_button = WIDGET_BUTTON(basex2, VALUE='ATV window', EVENT_PRO='pspec2_atv')

; Jump
;************
base_jump = WIDGET_BASE(base5, /ROW, /FRAME)
text_slit1 = WIDGET_LABEL(base_jump, VALUE='Jump to slit:')
text_slit2 = WIDGET_TEXT(base_jump, VALUE='', XSIZE=10, EVENT_PRO='pspec2_jump', /EDITABLE)

; Redshift quality flags
;***********************
  IF zqq EQ 1 THEN $
    slitcount = indx_pos[0]

base6 = WIDGET_BASE(base, /ROW)
values = ['Non-galactic', 'Junk', 'Q = 2', 'Q = 3', 'Q = 4'] 
redshift_buttons = CW_BGROUP(base6, values, /ROW, /FRAME, $
    LABEL_LEFT='Quality Selection', EVENT_FUNC='pspec2_setz', /RETURN_INDEX)
red_str = 'Current settings:  Q = ' + $
    STRING(fin_quality[slitcount], FORMAT='(I1)') + $
    '  z = ' + STRING(fin_redshifts[slitcount], FORMAT='(F5.3)') + $
    '   ' + fin_comments[slitcount]
finprevious_data = WIDGET_LABEL(base6, VALUE=red_str, /DYNAMIC_RESIZE)


; Comments about the spectrum
;*****************************
base7 = WIDGET_BASE(base, /ROW)
values = ['None', 'Strong OII', 'Strong abs features', 'Emission (not OII)', $
  'Rotation candidate']
canned_list = WIDGET_DROPLIST(base7, VALUE=values, TITLE='Canned comments', $
    EVENT_PRO='pspec2_canned')
text_label = WIDGET_LABEL(base7, VALUE='Other comments:')
text_comment  = WIDGET_TEXT(base7, VALUE='', /EDITABLE, $
  EVENT_PRO='pspec2_comment')

  


; Forwards and backwards
;************************
base_x3 = WIDGET_BASE(base, /ROW)
button1 = WIDGET_BUTTON(base_x3, VALUE='Next', uvalue='NEXT', XSIZE=230)
button2 = WIDGET_BUTTON(base_x3, VALUE='Back', uvalue='BACK', XSIZE=230)

; Finish session
;****************
button3 = WIDGET_BUTTON(base_x3, VALUE='Done', uvalue='DONE', XSIZE=230)





WIDGET_CONTROL, base, /REALIZE

WIDGET_CONTROL, draw1, GET_VALUE=drawID1
WIDGET_CONTROL, draw2, GET_VALUE=drawID2
WIDGET_CONTROL, draw3, GET_VALUE=drawID3
WIDGET_CONTROL, draw4, GET_VALUE=drawID4
WIDGET_CONTROL, draw5, GET_VALUE=drawID5
WIDGET_CONTROL, draw6, GET_VALUE=drawID6
WIDGET_CONTROL, draw7, GET_VALUE=drawID7
WIDGET_CONTROL, draw8, GET_VALUE=drawID8
text_widg = text1
text_widgz = text2
WIDGET_CONTROL, base, SET_UVALUE=drawID1


;***++++++++++++++++++++++++++++++++++++++++++++++++++
; Read and plot the first spectrum
;
;WSET, drawID1

  IF zqq EQ 1 THEN $
    slitcount = indx_pos[0]

  DEIMOS_READOBJ, maskname, objnames[slitcount], lambda=lambda, spec=spec, $
     ivar=ivar, zans=zans, synspec=synspec, fulldir=fulldir

  SELECT_Z, zall, objnames[slitcount], zres=zres

  redshift = zans.z  

  MY_PLOT_SPEC
  MY_SHOW_TEXT, zans, zres
  oldspec = spec
  oldsynspec = synspec
  postage_stamp, maskname, slitnames[slitcount], redshift
;
;***+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


; Start the GUI
;****************
XMANAGER, 'pspec2', base


cntn1 = WHERE(fin_quality EQ -1, cn1)
cnt1 = WHERE(fin_quality EQ 1, c1)
cnt2 = WHERE(fin_quality EQ 2, c2)
cnt3 = WHERE(fin_quality EQ 3, c3)
cnt4 = WHERE(fin_quality EQ 4, c4)
comp = 100.0*FLOAT(c3+c4)/FLOAT(cn1+c1+c2+c3+c4)
comp_pot =  100.0*FLOAT(c2+c3+c4)/FLOAT(cn1+c1+c2+c3+c4)

PRINT, '**********************************************'
PRINT, '************ PSPEC2 COMPLETE *****************'
PRINT, '**********************************************'
PRINT, ''
PRINT, ' Your selections are as follows:'
PRINT, '  Stars etc: ', cn1
PRINT, '  Q = 1:     ', c1
PRINT, '  Q = 2:     ', c2
PRINT, '  Q = 3:     ', c3
PRINT, '  Q = 4:     ', c4
PRINT, ' % completeness (Q > 2): ', comp
PRINT, ''
PRINT, '**********************************************'
PRINT, 'Would you like to save these results? (y/n):'

REPT:
  ans = STRUPCASE(GET_KBRD(1))
  PRINT, ans
  CASE ans OF
    'N': RETURN
    'Y': BEGIN
      filex = ''
      WHILE(1) DO BEGIN
       READ, filex, PROMPT='Enter filename to write to: '
       IF filex NE '' THEN BREAK
      ENDWHILE
      OPENW, 3, filex
      SPAWN, 'whoami', foo, /NOSHELL
      PRINTF, 3, 'Maskfile: ', maskname, '   user = ', foo 
      PRINTF, 3, 'Completeness = ', STRING(comp, FORMAT='(F6.2)'), '%   ',  '(potentially:', STRING(comp_pot, FORMAT='(F6.2)'), ' with Q>1)'
      FOR i=0, nslit_tot-1 DO $
        PRINTF, 3, FORMAT='(I9,I6,F9.5,I4,A7,A0)', $
           objnames[i], slitnames[i], fin_redshifts[i], $
           fin_quality[i], fin_canned[i], '  '+fin_comments[i]
      CLOSE, 3
    END
    ELSE: BEGIN
      PRINT, 'ERROR: Enter y or n'
      GOTO, REPT
    END
  ENDCASE



END






