;+
;
; NAME
;      reduce1d.pro
;
; PURPOSE
;      A wrapper routine to organize and run the 1d analysis
;      pipeline. 
;
; SYNTAX
;      reduce1d, [ result, files=files, /doplot, /debug, allresult= ]
;
; INPUTS
;      files = an optional parameter giving a vector of spec1d file
;              names. 
;
; KEYWORDS
;      /TODAY -- if set, the date in the filename will be today's
;                date, not the mask creation date
;      /nlsky - this keyword instructs the reduce1d to process the
;               non-local sky 1-d extractions. if this keyword is NOT
;               set or is set to zero, then the nlsky extractions are
;               ignored. if this keyword is set to 2, then ONLY the
;               nlsky extractions are processed. and for all others
;               values of nlsky, then both local and nlsky extractions
;               are processed. The nlsky results are written to
;               separate HDUs in the zresult output file.
;
; OUTPUTS
;     result -- output structure giving result of last object in call
;               (or single object)
;   
; PROCEDURES CALLED 
;      splog
;      zfind 
;      zrefind
;      vdispfit
;      fill_gap
;      cpbackup
;      dfpsplot
;      concat_dir
;      remove_telluric
;
; EXAMPLES
;
;
; COMMENTS
;
;
; HISTORY
;      Created August 7, 2002 by mcc.
;      various revisions 30sep02 MD
;      Revised March 10, 2002 by mcc - included ability to handle
;                                      nlsky reductions.
;
; TBD--
; 1. FITS keywords in output files
; 2. better PCA??
; 3. tools for better QA! 
;-

pro reduce1d, cresult, files=files, doplot=doplot, debug=debug, $
              allresult=allresult, today=today, nlsky=nlsky

; put into the spec1d log stderr file the location of the spec1d
; routines.
  findpro, 'reduce1d', dirlist=directories, /NOPRINT
  message, 'Running spec1d code in: '+ directories[0], /INFO
  print, 'Running spec1d code in: ' + directories[0]

; define a variable file_created which determines whether or not to
; initialize and create the zresult output FITS file.
  file_created = 0


;------------------------------------------------------------------------
; check how the nlsky keyword was set and accordingly call reduce1d
; recursively.
  if keyword_set(nlsky) then begin
      case nlsky[0] of
          0: nlsky = 0
          2: nlsky = 1
          else: begin
              reduce1d, fakevar, files=files, $
                doplot=doplot, debug=debug, allresult=allresult
              file_created = 1
              nlsky = 1
          end
      endcase
  endif else nlsky = 0

; fit skytweaks, once.
  if nlsky eq 0 then fittweaks

; print line to log file if doing nlsky work.
  if nlsky then print, 'Running on nlsky extractions....'

; check if 1-d spectrum files (spec1d files) were passed by user. if
; not, then search the current directory for them.
  if n_elements(files) eq 0 then files = findfile('spec1d.*.fits*')
  nfiles = n_elements(files)
; check that the file names are strings.
  if size(files, /tname) ne 'STRING' then $
    message, 'File names must be strings!'
; check that spec1d files were found.
  if nfiles eq 0 then message, 'No 1d files specified by user ' + $
    'or found in current directory!'

; setup housekeeping: set file names, directories, etc.
; read the header of the first spec1d file.
  hdr = headfits(files[0], ext=1, /silent)
; extract the mask name from the header.
  maskname = strcompress(sxpar(hdr, 'SLMSKNAM'), /rem)
; get mask name from filename, if necessary.
  if maskname eq 0 then maskname = (strsplit(files[0], '.', /extract))[1]
; remove suufix to mask (e.g. trim 4207.W to 4207).
  maskname = (strsplit(maskname, '.', /extract))[0]
; extract the observation date from the header.
  date = strcompress( sxpar(hdr, 'DATE-OBS'), /rem)
; if the today keyword is set, then define the date to be today's date
; rather than the observation date.
  if keyword_set(today) then get_date, date
; determine the output directory and name of the log file.
  outdir = getenv('D2_RESULTS') + '/zresult/'
  outname = maskname + '.' + date
  if nlsky then logfile = outdir + 'log/spRedux.' + outname + '.nlsky.log' $
  else logfile = outdir + 'log/spRedux.' + outname + '.log'
; make sure that the ps and log subdirectories exist.
  spawn, 'mkdir -p ' + concat_dir(outdir, 'log')
  spawn, 'mkdir -p ' + concat_dir(outdir, 'ps')

; open a post-script plotting file.
  if keyword_set(doplot) and not(keyword_set(debug)) then begin
      if nlsky then plotfile = outdir + 'ps/spPlot.' + outname + '.nlsky.ps' $
      else plotfile = outdir + 'ps/spPlot.' + outname + '.ps'
      dfpsplot, plotfile, /color, /landscape
  endif

; check the system time and begin analysis.
  stime0 = systime(1)
; open the log file for writting. only generate the log file if not in
; debug mode.
  if not(keyword_set(debug)) then begin 
      splog, filename=logfile
      splog, 'Log file ' + logfile + ' opened ' + SYSTIME()
  endif
; if a plot file will be generated, then specify the file name in the
; log file.
  if (keyword_set(plotfile)) then $
    splog, 'Plot file ' + plotfile
; write the IDL version and operating system to the log.
  splog, 'IDL version: ' + string(!version, format='(99(a," "))')
  spawn, 'uname -a', uname
  splog, 'UNAME: ' + uname[0]
; write the version numbers for the idlutils and spec1d packages to
; the log file.
;  splog, 'idlspec2d version ' + idlspec2d_version()
  splog, 'idlutils version ' + idlutils_version()

; specify the plot title for the post-script plots. WHY DO WE DO THIS HERE?
  plottitle = 'Galaxy Redshift'
; get the directory and file name for the galaxy and stellar templates.
  eigendir = concat_dir( GETENV('IDLSPEC1D_DIR'), 'templates' )
  eigenfile_gal = 'spDEEP.fits'
  eigenfile_star = 'spEigenStarDeep2.fits'
;  eigenfile = 'spEigenGal-*.fits'
;  eigenfile = 'emi_line.fits'

; select the stars eigen-file here to determine how many templates are in it.
  allfiles = findfile(djs_filepath(eigenfile_star, root_dir=eigendir), $
                      count=ct)
  if (ct eq 0) then $
    message, 'Unable to find EIGENFILE matching ' + eigenfile_star
  eigenfile_star = fileandpath(allfiles[ (reverse(sort(allfiles)))[0] ])
  shdr = headfits(djs_filepath(eigenfile_star, root_dir=eigendir))
  nstar = sxpar(shdr, 'NAXIS2') > 1
  subclass = strarr(nstar)      ;types of stars
  for istar=0, nstar-1 do $
    subclass[istar] = $
    strtrim( sxpar(shdr, 'NAME'+strtrim(string(istar),2)), 2)
  
; write info to log file.
  splog, 'Compute GALAXY redshifts: ', $
    ' ZMIN=', zmin, ' ZMAX=', zmax, ' PSPACE=', pspace
; determine airmass value for the header of the 0th spec1d file.
  airmass = sxpar(hdr, 'AIRMASS')
  splog, 'Mean airmass for this mask: ', airmass 

;--------------------------------------------------
; loop over all of the spec1d files.
  for ii=0,nfiles-1 do begin
; check the system time at start, so that we can determine the time
; needed for processing each individual spec1d file.
      t0 = systime(1)
; grab the mask name and slit number from the file
; header. alternatively, these values can be determined from the
; spec1d file name (if the header fails to contain the info).
      hdr = headfits(files[ii], ext=1, /silent)
      maskname = strcompress(sxpar(hdr, 'SLMSKNAM'), /rem)
      slitname = strcompress(sxpar(hdr, 'SLITNO'), /rem)
      objname = strcompress(sxpar(hdr, 'OBJNO'), /rem)
; grab the slitpa and maskpa values from the header of the spec1d
; file.
      slitpa = sxpar(hdr, 'SLITPA')
      maskpa = sxpar(hdr, 'MASKPA')

; get the maskno, slitno, objno from the filename.
;      pieces = strsplit(files[ii], '.', /extract)
;      if maskname eq 0 then maskname = pieces[1]
;      if slitname eq 0 then slitname = pieces[2]
;      objname = pieces[3]

; get the observation date from the header.
      date = strcompress(sxpar(hdr, 'DATE-OBS'), /rem)
      mjd = double(sxpar(hdr, 'MJD-OBS'))
; trim the mask name.
      deimos_isdeep, isdeep, maskname
; define the plottitle.
      plottitle = 'rchi2 for mask: ' + maskname +' slit: ' + slitname
      
; fill the gap between the blue and red portions of the 1-d spectra.
      if strpos(maskname, 'KTRS') ge 0 then boxsprof = 1 else boxsprof = 0
      if nlsky then ss1d = fill_gap(files[ii], /nlsky, boxsprof=boxsprof,/tweak) $
      else ss1d = fill_gap(files[ii], boxsprof=boxsprof,/tweak)
; check that the fill_gap routine succeeded in filling the gap.
      if size(ss1d, /tname) ne 'STRUCT' then begin
          print, 'Skipping slit ' + slitname + '!'
          result = replicate({zresult}, 10)
          result.objname = objname
          result.slitname = slitname
          result.maskname = maskname
      endif else begin
          if total(ss1d.spec) eq 0 then begin
              print, 'Skipping slit ' + slitname + ', no data!'
              result = replicate({zresult}, 10)
              result.objname = objname
              result.slitname = slitname
              result.maskname = maskname
          endif else begin
; correct for telluric absorption bands (do this before shifting to
; the vacuum wavelengths!).
          remove_telluric, ss1d, airmass
	  fix_response,ss1d
          
; convert lambda values from air to vacuum.
          lambda = ss1d.lambda 
          airtovac, lambda    
          ss1d.lambda = float(lambda) ;replace, reconvert from double

; estimate signal-to-noise in 1-d spectrum.
          s2n = mean(ss1d.spec*sqrt(ss1d.ivar))
          splog, ' '
          splog, 'mean S2N for frame: ', files[ii], ' : ', s2n
          
;;; SMOOTH THE 1d SPECTRUM - needed to ease interpolation.
;   should we do this before we fill_gap?
;  WE were getting rchi2<1 with this smoothing.  Eliminate it!!
;    ss1d.spec = ivarsmooth(ss1d.spec, ss1d.ivar, 3, ss1d.ivar)
          
; check for NaNs and Infs...interpolate over them.
          nfn = where(finite(ss1d.spec) eq 0, nfnum, compl=fn)
          if nfnum gt 0 then begin 
              ss1d.spec[nfn] = interpol(ss1d.spec[fn], ss1d.lambda[fn], $
                                        ss1d.lambda[nfn])
          endif

; ------------------------------
; determine galaxy redshifts.
; first set the parameters for call to zfind.
          zmin = -0.0001
          zmax = max(ss1d.lambda)/3727. -1.  ;allows for variable end
          pspace = 5            ; was 1
          width = 5*pspace
          nfind = 5
          npoly = 0             ; 3
          columns = [0,  1, 2]     ;include absorption and emission templates in all cases
          
; determine the redshift of object by comapring to the galaxy
; templates.
          result = zfind(ss1d, /linear_lambda, eigenfile=eigenfile_gal, $
                         eigendir=eigendir, npoly=npoly, $
                         zmin=zmin, zmax=zmax, nfind=nfind, pspace=pspace, $
                         width=width, plottitle=plottitle, doplot=doplot, $
                         debug=debug, objflux=objflux, objivar=objivar, $
                         loglam=loglam, columns=columns)
          
; redo the analysis at higher sampling, just around 5 minima chisqr regions
          pspace = 1
;            splog, 'Locally re-fitting GALAXY redshifts'
          res_gal = zrefind(ss1d, objflux, objivar, hdr=hdr, pwidth=15, $
                            pspace=pspace, width=3.*pspace, zold=result, $
                            loglam=loglam, plottitle=plottitle,  $
                            doplot=doplot, debug=debug, columns=columns)
; calculate the delta chi^2 values from minimum to minimum.
          result = res_gal
          delta_chisqr = (result[1].rchi2-result[0].rchi2)*.5* $
            (result[1].dof + result[0].dof)
          splog, 'Minimum chisqr is better than next min. by: ', delta_chisqr
          result[0].rchi2diff = result[1].rchi2-result[0].rchi2
          result.class = 'GALAXY' 
          result.subclass = ' '
          
; now find the velocity dispersions aacording to the galaxy models.
;            splog, 'Find velocity dispersions for galaxies'
          vdisp_wrapper, files[ii], zresult=result, $
            nfind=nfind, airmass=airmass, nlsky=nlsky
         
;            splog, 'CPU time to fit GALAXY velocity dispersions = ', systime(1)-t0
          
; -------------------
; Find STAR redshifts - do all stellar templates with one call.

          npoly = 0             ; With only 1 eigen-template, fit 3 poly terms as well.
          zmin = -0.004         ; -1200 km/sec
          zmax = 0.004          ; +1200 km/sec
          pspace = 1
          nfind = 1
; check the system time.
          ts0 = systime(1)
            
;            splog, 'Compute STAR (' + subclass + ') redshifts:', $
;              ' ZMIN=', zmin, ' ZMAX=', zmax, ' PSPACE=', pspace
          res_star = zfind_star(ss1d, /linear_lambda, $
                                eigenfile=eigenfile_star, npoly=npoly, $
                                zmin=zmin, zmax=zmax, pspace=1, $
                                nfind=nfind, width=5*pspace, $
                                subclass=subclass, doplot=0, debug=debug)
          
          res_star.class = 'STAR'
          res_star = res_star[sort(res_star.rchi2)] ;sort by rchi2 
          
          result = [result, res_star[0:2]] ; Append result of top 2 star

          splog, 'CPU time to compute STAR redshifts = ', systime(1)-ts0
 
   ;----------
   ; Find QSO redshifts


          npoly = 0
          zmin = 0.0033         ; +1000 km/sec
;          zmax = max(ss1d.lambda)/1215. -1  ; Max range we can expect to see
          zmax = 5.
          pspace = 10
          nfind = 2             ;find best QSO candidate z
          plottitle = 'QSO Redshift'

          eigenfile = 'spEigenQSOdeep.fits'

          splog, 'Compute QSO redshifts:', $
                 ' ZMIN=', zmin, ' ZMAX=', zmax, ' PSPACE=', pspace
          t0 = systime(1)
          
          res_qso =  zfind_qso(ss1d, /linear_lambda, $
                      eigenfile = eigenfile, npoly = npoly, zmin = zmin, $
                      zmax = zmax, pspace = pspace, loglam = loglam, $
                      nfind = nfind, width = 5*pspace, objflux = objflux, $
                      objivar=objivar, plottitle = plottitle, $
                       doplot = doplot, debug = debug )


          splog, 'CPU time to compute QSO redshifts = ', systime(1)-t0


;skip refitting QSO', as not necessary
;
;          splog, 'Locally re-fitting QSO redshifts'
;         t0 = systime(1)
;          res_qso = zrefind(ss1d, objflux, objivar, hdr=hdr, pwidth=91, $
;                            pspace=1, width=3.*pspace, zold=res_qso, $
;                            loglam=loglam, plottitle=plottitle, $
;                            doplot=doplot, debug=debug  )
;
;
;          splog, 'CPU time to re-fit QSO redshifts = ', systime(1)-t0

          res_qso.class = 'AGN'
          res_qso.subclass = ' '

; don't append QSO results here. we don't want them to be sorted by
; chi2 value. instead, we will sort the GALAXY and STAR z values and
; then append the 2 QSO redshifts.
;          result = [result, res_qso] ; Append results




          
; now decide which result gives best rchisqr
          
;----------
          

   ;----------
   ; Sort results for each object by ascending order in chi^2/DOF,
   ; but putting any results with zero degrees-of-freedom at the end.

          minvdiff = 1000.0     ; km/s
          cspeed = 2.99792458e5
          rchi2 = result.rchi2
          
          isort = sort(rchi2 + (result.dof EQ 0)*max(rchi2))
          result = result[isort]

; append the QSO results here! -- so that they aren't sorted by chi2!
          result = [result, res_qso]
          nper = (size(result,/dimens))[0]
; Find the difference in reduced chi^2 between each result and the next
          rchi2 = result.rchi2
          for ia=0,nper-2 do begin
              inext = (where(abs(result[ia+1:nper-1].z - result[ia].z) GT $
                        minvdiff/cspeed AND result[ia+1:nper-1].dof GT 0))[0]
              if (inext ne -1) then $
                result[ia].rchi2diff = rchi2[ia+1+inext] - rchi2[ia]
          endfor

; put the slitname, maskname, etc. into the zresult structures!         
          result.objname = objname
          result.maskname = maskname
          result.slitname = slitname
          result.date = date
          result.mjd = mjd
;          result.slitpa = slitpa
          endelse
      endelse                   ;end query on selected object

; define how many results to keep.
      nres = (n_elements(result) < 10) - 1          

; form result structure for output-- will evolve as Z's become real
      if not(keyword_set(cresult)) then begin
          cresult = result[0]   ;if first in list
          allresult = result[0:nres] ;save best nres fits
      endif else begin
          cresult = [[cresult], [result[0]]] ;if not first, append
          allresult = [[allresult], [result[0:nres]]]
      endelse
      
; write out how long it took to determine z values.
      splog, 'CPU time for this Object = ', SYSTIME(1)-t0

        
  endfor                        ;end loop over input files
    
;----------
; if not in debug mode, then write output file.
  if not(keyword_set(debug)) then begin
; make header for 0th extension.
      slitfiles = findfile('slit.*.fits*', count=nfiles)
      if nfiles eq 0 then mkhdr, hdr, 0, /extend $
      else begin
          mkhdr, hdr, 0, /extend
          slithdr = headfits(slitfiles[0], ext=1, /silent)
          keywords = ['DATE', 'OBJECT', 'ROTATVAL', 'DATE-OBS', $
                      'UT', 'AIRMASS', 'TARGNAME', 'EQUINOX', $
                      'RA', 'DEC', 'AZ', 'EL', 'HA', 'ST', 'UTC', $
                      'MJD-OBS', 'PARANG', 'SLMSKNAM', 'GRATEPOS', $
                      'EXPTIME', 'OBSERVER', 'SYNOPSIS', $
                      'HPLOGTIM', 'SP2DVERS', 'AUTHOR']
          for qq=0,n_elements(keywords)-1 do begin
              value = sxpar(slithdr, keywords[qq])
              sxaddpar, hdr, keywords[qq], value
          endfor
      endelse
      splog, 'Writing output file.....'
; put the spec1d version number in the header.
      dir1d = strsplit(getenv('IDLSPEC1D_DIR'), '/', /extract)
      num = n_elements(dir1d)
      if num gt 1 then $
        sxaddpar, hdr, 'SP1DVERS', dir1d[num-2], 'Version of spec1d' $
      else sxaddpar, hdr, 'SP1DVERS', '', 'Version of spec1d'
;      sxaddpar, hdr, 'VERS1D', spec2d_version(), $
;        'Version of idlspec2d for 1D reduction', after='VERSCOMB'
; add the maskpa to the zresult file header.
      spec1dhdr = headfits(files[0], ext=1, /silent)
      maskpa = sxpar(spec1dhdr, 'MASKPA')
      sxaddpar, hdr, 'MASKPA', maskpa, 'Slit Mask PA on sky w.r.t. N'

      spawn, 'uname -n', uname
      sxaddpar, hdr, 'UNAME', uname[0]

      zbestfile = outdir + 'zresult.' + outname + '.fits'

; put the header in the zeroth hdu.
      if not(file_created) then mwrfits, 0, zbestfile, hdr, /create 

; set the comment tag to none....if we leave it as an empty string,
; then mrdfits will not read the tag in!!!
      cresult.comment = 'none'
      allresult.comment = 'none'

      if nlsky then begin
          fxbhmake, hdr, 1, 'zbest-nl', /date, /initialize
          mwrfits, cresult, zbestfile, hdr
          fxbhmake, hdr, 1, 'zall-nl', /date, /initialize
          mwrfits, allresult, zbestfile, hdr
      endif else begin
          fxbhmake, hdr, 1, 'zbest', /date, /initialize
          mwrfits, cresult, zbestfile, hdr
          fxbhmake, hdr, 1, 'zall', /date, /initialize
          mwrfits, allresult, zbestfile, hdr
      endelse

;   mwrfits, zans, zbestfile
;   mwrfits, synflux, zbestfile
;   mwrfits, dispflux, zbestfile

;   sxaddpar, hdr, 'DIMS0', nper, ' Number of fits per objects'
;   sxaddpar, hdr, 'DIMS1', nobj, ' Number of objects'
;   mwrfits, 0, zallfile, hdr, /create ; Retain the original header in first HDU
;   mwrfits, res_all, zallfile

  endif


;   if (keyword_set(debugfile)) then dfpsclose

   ;----------
   ; Generate final QA plots

;   splog, 'Generating QA plots'

;   if (keyword_set(plotfile)) then begin
;      cpbackup, plotfile
;      dfpsplot, plotfile, /color
;   endif

;   plottitle = string(zans[0].plate, zans[0].mjd, $
;    format='("Flux-Calibration Errors Plate=", i4, " MJD=", i5)')
;   qaplot_fcalibvec, objloglam, objflux, objivar, synflux, plugmap, zans, $
;    plottitle=plottitle

;   if (keyword_set(plotfile)) then dfpsclose

   ;----------
   ; Close log file

    if keyword_set(doplot) then dfpsclose ;close the ps file
    
    splog, 'Total time for SPREDUCE1D = ', systime(1)-stime0, ' seconds', $
      format='(a,f6.0,a)'
    splog, 'Successful completion of SPREDUCE1D at ', systime()
    if (keyword_set(logfile)) then splog, /close


; make a done-processing file to signify the completetion of the
; spec1d pipeline.
  openw, 2, 'donespec1d.txt'
  printf, 2, spec1d_version()
  close, 2


; now run the rotation curve code. to do this, find the most recent
; zcat file and run with it.
;    d2dir = getenv('D2_RESULTS')
;    zpath = concat_dir(d2dir, 'zresult/zcat*.fits*')
;    zcat = findfile(zpath, count=ncats)
;    if ncats eq 0 then begin
;        print, 'ERROR: No zcat found!'
;        print, 'rcurve pipeline canNOT run w/o zcat!'
;    endif else begin
;        zcat = zcat[(reverse(sort(zcat)))[0]]
;        do_rcurve, zcat=zcat, /doplot
;    endelse


end








