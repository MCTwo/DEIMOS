; top level control script to reduce a mask of DEIMOS data
; intended to be run in batch mode
; formerly known as reduce_it.pro
; This is for slider 4 data only, taken prior to May of 2009
; or subsequent to August 19 of 2009

; you must CD to the directory you want to run this in. 

; SET NLSKY to run nonlocal sky subtraction in addition to local

; MD & DPF

pro domask_slider4_Marusa, planfile, slitcat=slitcat, chiplist=chiplist, nlsky=nlsky, skipcalib=skipcalib, badarcs=badarcs, arcshift=arcshift	;BL added slitcat 5/12, BCL added badarcs and arcshift 7/14

; slitcat is an acsii input list of the slits that are desired to be reduced on the mask, i.e., 
; if there is a problem with certain slits they can be omitted from the reduction. To create the 
; proper ascii file do in IDL after the reduction has failed:
; a = mrdfits('*bintabs*',1)
; openw, lun, 'XXXslits.cat', /get_lun
; for i=0, n_elements(a.slitname)-1 do prinf, lun, a[i].slitname
; free_lun, lun
; And remove the offending slits

; badarcs is for masks that had mechanical failures during the calibration process (warping or inproper 
; insertion), and need more wavelength tolerance w.r.t. the optical model guess for catching the arcs 
; when generating the initial (pre-skytweaked) wavelength solution. It sets the lagrange in 
; discrete_tweak_omodel.pro to double the initial value. Fully tested now. The wavelength solution is precise
; and accurate. This is for the initial difference between the first coefficient in the wavelength solution 
; and the optical model guess, if it's > 20A, need to use this keyword.

; arcshift is for shifting the wavelength solution coming from the arclamp frame for the science frame.
; It is only to be used if a large offset in wavelength space between the science frames and arc frame is 
; noticed. This phenomenon appears to come from some mechanical failure in the masks or from some FCS issue
; (maybe? I have no idea). The shift is applied at the end of arc tracing right before the creation of the 
; calibSlit files, which is then applied to the spSlit files and to the rest of the data. Only a translational
; shift is applied, i.e., only the 0th order coefficient is changed, there is no warping assumed. arcshift 
; is defined in Angstroms as a float 

  if n_params() eq 0 then begin
    planfile = findfile('*.plan')
    planfile = planfile[0]
  endif

  t1 = systime(1)
  if NOT keyword_set(planfile) then message, 'You must specify a plan file!'

; put into the log the location of the spec2d files
  findpro,'domask',dirlist=directories,/NOPRINT

  message,'Running spec2d code in: '+directories[0],/INFO
  print,'Running spec2d code in: '+directories[0]
; read bin tables from raw frame and store as .fits file.
  make_bintab_file, planfile

  deimos_isdeep,isdeep
  if n_elements(nlsky) eq 0 then nlsky = (isdeep eq 1)


; write calibSlit files
  if n_elements(slitcat) ne 0 then begin		;BL added 5/12, slitlist is a list of long integers passed to deimos_mask_calibrate
     readcol, slitcat, format='L', slitlist
     print, 'WARNING: Only Reducing Slits ', slitlist
  endif	
 
  calib_test = findfile('*calib*.fits', count=filecount)			;Added BL 5/13, skip making calibSlit files if they're already made, will skip also if it crashed at a slit later than 50
  if n_elements(arcshift) then print, 'narfnarf', n_elements(arcshift), arcshift
  if n_elements(skipcalib) gt 0 then begin
  	print, 'narf', n_elements(skipcalib)
	if filecount lt 20 then begin								
	if n_elements(badarcs) then begin
		
		if n_elements(arcshift) then begin
			deimos_mask_calibrate_slider4_Marusa, planfile, /noplot, chiplist=chiplist, slitlist=slitlist, badarcs=badarcs, arcshift=archshift
		endif else begin
			deimos_mask_calibrate_slider4_Marusa, planfile, /noplot, chiplist=chiplist, slitlist=slitlist, badarcs=badarcs
		endelse
	endif else begin
		if n_elements(arcshift) then begin	
			deimos_mask_calibrate_slider4_Marusa, planfile, /noplot, chiplist=chiplist, slitlist=slitlist, arcshift=arcshift
		endif else begin
			deimos_mask_calibrate_slider4_Marusa, planfile, /noplot, chiplist=chiplist, slitlist=slitlist
		endelse
	endelse
	endif else begin
	print, 'Skiping calibration frames, you have written enough to satisfy me...'
	endelse
  endif else begin
	;print, 'UHHHHHHH.... this should print.... #############################33'
	if n_elements(badarcs) then begin
		if n_elements(arcshift) then begin
                        ;print, 'GET AN EXTINGUISHER!!!!'
			deimos_mask_calibrate_slider4_Marusa, planfile, /noplot, chiplist=chiplist, slitlist=slitlist, badarcs=badarcs, arcshift=arcshift
		endif else begin
			;print, 'THERES THE GIANTS FIRST HIT!!!!'
                        deimos_mask_calibrate_slider4_Marusa, planfile, /noplot, chiplist=chiplist, slitlist=slitlist, badarcs=badarcs
                endelse
	endif else begin
		if n_elements(arcshift) then begin      
                        ;print, 'HES SCUFFLING!!!!!'
			deimos_mask_calibrate_slider4_Marusa, planfile, /noplot, chiplist=chiplist, slitlist=slitlist, arcshift=arcshift
                endif else begin
			;print, 'AGAINST THE RANGERS!!!!'
                        deimos_mask_calibrate_slider4_Marusa, planfile, /noplot, chiplist=chiplist, slitlist=slitlist
                endelse
	endelse
  endelse

; write spSlit files
  if n_elements(badarcs) then begin
    deimos_2dreduce, planfile, /badarcs
  endif else begin
    deimos_2dreduce, planfile
  endelse

; combine multiple exposures
  list=findfile('spSlit*.fits')
  spslit_combine,list, nlsky = nlsky

  read_planfile, planfile, maskname, rawdatadir, outdatadir, $
             flatnames

  head_flat = headfits(flatnames[0])
  deimos_grating, head_flat, grating, grangle, lambda_c

  minlambda=lambda_c - 1300*(1200./grating) 
  bluelim=5500.-1500.*(minlambda lt 5500)

  slitfiles = findfile('slit*.fits', count=nfiles)
  slitfiles = slitfiles[sort(slitfiles)]
  isred = (strpos(slitfiles, 'R.fits') ne -1)
  trans = max(where(isred and lindgen(nfiles) le nfiles/2))
  if trans eq -1 then trans = floor(nfiles/2)

; Make simple 2d images, 2 per mask
  epos = strpos(slitfiles[0], '.fits')
  masknumber = strmid(slitfiles[0], 4, epos-4-4) ;get mask name '.xxxx.'
  image = slit_merge_lambda(slitfiles[0:trans], hdr,blue=bluelim)
  writefits, 'Allslits0'+masknumber+'fits', image, hdr
  if trans gt 1 then begin
    image = slit_merge_lambda(slitfiles[trans+1:nfiles-1], $
                              hdr, blue = bluelim)
    writefits, 'Allslits1'+masknumber+'fits', image, hdr
  endif


; make non-local allslits.
  if keyword_set(nlsky) then begin
    img = slit_merge_lambda(slitfiles[0:trans], hdr, /nonlocal,blue=bluelim)
    writefits, 'Allslits0' + masknumber + 'nonlocal.fits', img, hdr
    if trans gt 1 then begin
      img = slit_merge_lambda(slitfiles[trans+1:nfiles-1], $
                              hdr, /nonlocal, blue = bluelim)
      writefits, 'Allslits1' + masknumber + 'nonlocal.fits', img, hdr
    endif
  endif

; do 1d extraction.
  deimos_isdeep, isdeep, maskname
  slitfiles = findfile('slit*.fits', count=nfiles)
  if isdeep then $
    do_extract, files=slitfiles, nsigma_optimal=1.75, nsigma_boxcar=1.1 $
  else begin
      if strpos(strupcase(maskname), 'KTRS') ge 0 then $
        do_extract, files=slitfiles, nsigma_optimal=1.75, nsigma_boxcar=1.5 $
      else do_extract, files=slitfiles, nsigma_optimal=1.75, nsigma_boxcar=1.1+0.4*(grating eq 600.)
  endelse

; do 1d extraction w/ non-local sky
  if keyword_set(nlsky) then $
    do_extract, files=slitfiles, /nonlocal

  print
  print, systime()+'  Done with slit processing!!!'
  print
  print, 'Total time elapsed: ', (systime(1)-t1)/3600., ' hours.'

; compress the files.
  spawn, 'gzip -1 -vf *lit*.*.fits ' ;operations suspend

  print
  print, systime()+'  Done with gzipping!!!'
  print
  print, 'Total time elapsed: ', (systime(1)-t1)/3600., ' hours.'


; do quality assurance
  if isdeep then qa_check, /doplot

; make a done-processing file to signify the completetion of the
; spec2d pipeline.
  openw, 2, 'doneprocessing.txt'
  printf, 2, spec2d_version()
  close, 2

  exit
end







