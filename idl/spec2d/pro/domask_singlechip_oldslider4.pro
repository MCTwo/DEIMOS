; top level control script to reduce a mask of DEIMOS data
; intended to be run in batch mode
; formerly known as reduce_it.pro

; you must CD to the directory you want to run this in. 

; SET NLSKY to run nonlocal sky subtraction in addition to local

; MD & DPF

pro domask_singlechip_oldslider4, planfile, slitcat=slitcat, chiplist=chiplist, nlsky=nlsky, skipcalib=skipcalib	;BL added slitcat 5/12



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
 
  calib_test = findfile('/Volumes/Data2/orelse/lemaux/deimos/1429/1429A/*/calib*.fits', count=filecount)			;Added BL 5/13, skip making calibSlit files if they're already made, will skip also if it crashed at a slit later than 50
  if n_elements(skipcalib) gt 0 then begin
  	print, 'narf', n_elements(skipcalib)
	if filecount lt 20 then begin								
  	deimos_mask_calibrate_oldslider4, planfile, /noplot, chiplist=chiplist, slitlist=slitlist
	endif else begin
	print, 'Skiping calibration frames, you have written enough to satisfy me...'
	endelse
  endif else begin
  	deimos_mask_calibrate_oldslider4, planfile, /noplot, chiplist=chiplist, slitlist=slitlist
  endelse

; write spSlit files
  deimos_2dreduce, planfile


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
  ;deimos_isdeep, isdeep, maskname
  ;slitfiles = findfile('slit*.fits', count=nfiles)
  ;if isdeep then $
    ;do_extract, files=slitfiles, nsigma_optimal=1.75, nsigma_boxcar=1.1 $
  ;else begin
      ;if strpos(strupcase(maskname), 'KTRS') ge 0 then $
        ;do_extract, files=slitfiles, nsigma_optimal=1.75, nsigma_boxcar=1.5 $
      ;else do_extract, files=slitfiles, nsigma_optimal=1.75, nsigma_boxcar=1.1+0.4*(grating eq 600.)
  ;endelse

; do 1d extraction w/ non-local sky
  ;if keyword_set(nlsky) then $
    ;do_extract, files=slitfiles, /nonlocal

  print
  print, systime()+'  Done with slit processing!!!'
  print
  print, 'Total time elapsed: ', (systime(1)-t1)/3600., ' hours.'

; compress the files.
  ;spawn, 'gzip -1 -vf *lit*.*.fits ' ;operations suspend

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







