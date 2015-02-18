;+
; NAME:
;   zfind
;
; PURPOSE:
;   Find possible redshift matches for a set of spectra using a set of
;   eigen-templates. 
;     (this version does NOT save a filtered version of template)
;
; CALLING SEQUENCE:
;   result = zfind( objflux, objivar, hdr=hdr, $
;    eigenfile=, [eigendir=, columns=, npoly=, linear_lambda=, $
;     zmin=, zmax=, zguess=, pwidth=, nfind=, width=, $
;     objflux=, objinv=,  _EXTRA= ]
;
; INPUTS:
;   objflux    - Object fluxes [NPIXOBJ,NOBJ]
;   objivar    - Object inverse variances [NPIXOBJ,NOBJ]
;
; REQUIRED KEYWORDS:
;   hdr        - FITS header for objects, used to construct the wavelengths
;                from the following keywords: COEFF0, COEFF1.
;   eigenfile  - Input FITS file with an [NPIXSTAR,NSTAR] image with
;                either templates or eigenspectra.  If a wildcard appears
;                in the file name, then the file that appears last in a sort
;                is used.
;                The header keywords COEFF0, COEFF1 are used to specify
;                the wavelength mapping in log-10 Angstroms.
;
; OPTIONAL KEYWORDS:
;   eigendir   - Directory for EIGENFILE; default to $IDLSPEC2D/templates.
;   columns    - Column numbers of the eigenspectra image to use in the
;                PCA fit; default to all columns.
;   npoly      - Number of polynomial terms to append to eigenspectra;
;                default to none.
;   zmin       - Minimum redshift to consider; default to no lower bound.
;   zmax       - Maximum redshift to consider; default to no upper bound.
;   zguess     - Initial guess for redshift; search for a solution about
;                this value.  If specified with PWIDTH, then ZMIN and ZMAX
;                are ignoreed.
;   pwidth     - Search width in pixels about the intial guess redshift ZGUESS.
;                If specified with ZGUESS, then ZMIN and ZMAX are
;                ignored.
;   linear_lambda - if set, input data is linear and will be
;                   transformed to log lambda
;   nfind      - Keyword for ZCOMPUTE().
;   width      - Keyword for ZCOMPUTE().
;   linear_lambda - flag set if input data is linear in
;                   lambda. triggers the mapping to loglambda
;   objflux      -log lambda data; if present, data is not read
;                  to signify data is passed, in correct units
;   obivar       -goes with objflux
;         
;   _EXTRA     - Keywords for ZCOMPUTE(), such as PSPACE, DOPLOT, DEBUG.
;
; OUTPUTS:
;   result     - Structure with redshift-fit information.  Structure
;                elements are left blank if fewer than NFIND peaks are found.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   One can specify a search domain for the redshift with ZMIN and ZMAX, or
;   with ZGUESS and PWIDTH.  If none of those parameters are set, then all
;   possible redshifts that overlap the object and star (template) are tested.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   concat_dir()
;   djs_filepath()
;   fileandpath()
;   readfits()
;   splog
;   sxpar()
;   zcompute()
;
; INTERNAL SUPPORT ROUTINES:
;   sp1d_struct()
;
; REVISION HISTORY:
;   28-Jun-2000  Written by D. Schlegel, Princeton
;   Revised June 27, 2002 by mcc:
;      Added additional parameter linear_lambda= so that a vector 
;   of wavelength values can be passed to routine. This way the 
;   wavelength info isn't passed via the header (hdr) parameter.
;   Commented out the portion where the hdr is referenced and 
;   added new lines which take input linear lambda vector and 
;   switch to log-lambda values. 
;------------------------------------------------------------------------------
function zfind_star, ss1d, eigenfile=eigenfile, eigendir=eigendir, $
                npoly=npoly, zmin=zmin, zmax=zmax, $
                zguess=zguess, pwidth=pwidth, nfind=nfind, width=width, $
                linear_lambda=linear_lambda,objflux=objflux, $
                objivar=objivar, loglam=loglam, subclass=subclass, $
                _EXTRA=EXTRA

common com_zfind_star, starflux_in, starloglam0, stardloglam, $
     nstars, npoints, thisfile,  starflux_corrected
       
  if NOT keyword_set(starflux_corrected) then begin 
;if first call, read template

;;; CHECK IF THE TEMPLATE DIRECTORY WAS SUPPLIED BY USER. IF NOT, USE
;;; DEFAULT DIRECTORY TO FIND TEMPLATE FILES.
    IF (N_ELEMENTS(eigendir) EQ 0) THEN $
      eigendir = concat_dir( GETENV('IDLSPEC1D_DIR'), 'templates' )

;;; GET THE MOST RECENT EIGENFILE.
    allfiles = FINDFILE(djs_filepath(eigenfile, root_dir=eigendir), COUNT=ct)
    IF ct EQ 0 THEN $
      message, 'Unable to find EIGENFILE matching '+eigenfile
    thisfile = allfiles[ (REVERSE(SORT(allfiles)))[0] ]
;;; WRITE TO LOG NOTING WHICH TEMPLATE WE ARE USING.
;   splog, 'Selecting EIGENFILE=' + thisfile

    starflux_in = readfits(thisfile, shdr, /SILENT)
;;; FROM THE HEADER, GRAB THE WAVELENGTH INFORMATION.
    starloglam0 = sxpar(shdr, 'COEFF0')
    stardloglam = sxpar(shdr, 'COEFF1')
    npoints = n_elements(starflux_in[*, 0])
    starloglam = findgen(npoints)*stardloglam + starloglam0

;    print, 'dloglam for stars: ', stardloglam

    nstars = (size(starflux_in, /dimen))[1] ;how many stars?
; remove continuum from template at this level

    starflux = starflux_in*0.

    for i=0, nstars-1 do begin ;correct all stars
      starcont = djs_median(starflux_in[*, i], width=2500, boundary='reflect')
      starflux[*, i] =  starflux_in[*, i] - float(starcont)

    endfor
    starflux_corrected = starflux ;save for next call

  endif else starflux = starflux_corrected ;if not 1st call, restore template.

; DETERMINE GRID SIZE (IN LOG LAMBDA) FOR THE OBJECT SPECTRUM AND
; DETERMINE THE WAVELENGTH RANGE FOR THE OBJECT. Force binning to 
; match that of template 
  IF KEYWORD_SET(linear_lambda) THEN $
    loglam = linear2log(ss1d, binsize=stardloglam, flux=objflux, ivar=objivar)
  objloglam0 = loglam[0]
  objdloglam = stardloglam; loglam[1] - loglam[0] 
;should be same as stardloglam, if code works
;  print, ' log lam range: ', loglam[0], loglam[n_elements(loglam)-1]

;remove continuum same as for templates
  objcont = djs_median(objflux, width=2500, boundary='reflect')
  objflux = objflux - float(objcont)


;;; CHECK IF THE zmin AND zmax ARGUMENTS WERE PASSED. IF SO, THEN
;;; CONVERT THE REDSHIFT VALUES INTO PIXEL VALUES pmin AND pmax.
;;; THIS IS ONLY TRUE IF objloglam0 = temploglam0?
  IF n_elements(zmin) NE 0 THEN $
    pmin = FLOOR( ALOG10(1.0 + zmin) / objdloglam )
  IF n_elements(zmax) NE 0 THEN $
    pmax = CEIL( ALOG10(1.0 + zmax) / objdloglam )


;;; CHECK IF A GUESS REDSHIFT zguess WAS PASSED ALONG WITH A PIXEL
;;; WINDOW pwidth. IF SO, THEN RESET pmin AND pmax ACCORDING TO THE
;;; GUESS VALUE AND THE WINDOW.
  IF N_ELEMENTS(zguess) GT 0 AND KEYWORD_SET(pwidth) THEN BEGIN
    IF KEYWORD_SET(width) THEN width1 = width $
    ELSE width1 = pwidth
      pmin = FLOOR( ALOG10(1.0 + zguess) / objdloglam - 0.5*(pwidth+1+width1))
      pmax = FLOOR( ALOG10(1.0 + zguess) / objdloglam + 0.5*(pwidth+1+width1))
  ENDIF

; if pmax is too large, reset it!
  maxp = fix((objloglam0 + objdloglam*.99*n_elements(objflux) - starloglam0)/ $
          objdloglam)  
;  print, 'test: ', pmax, maxp
  if maxp lt pmax then begin ;limit upper redshift range to have overlap
     pmax = maxp
     print, 'resetting pmax to: ', maxp 
  endif


  IF abs(objdloglam - stardloglam) GT 0.05*objdloglam THEN $
    MESSAGE, 'Template and object lambda resolution do NOT match!'

   ;----------
   ; Compute the redshift difference between the first pixel of the object
   ; spectra and the template.
   poffset = (objloglam0 - starloglam0) / objdloglam

;  print, 'poffset, pmin,pmax :', poffset, pmin, pmax


   for istar=0,  nstars-1 do begin ;loop over stellar templates

     star = starflux[*, istar]
;;; IF THE OPTIONAL PARAMETER npoly IS PASSED, THEN ADD THE
;;; APPROPRIATE NUMBER OF POLYNOMIAL TERMS TO THE TEMPLATE ARRAY.
     IF KEYWORD_SET(npoly) THEN $
       star = [ [star], [poly_array(npoints,npoly)] ]



;;; CALL zcompute.pro TO COMPUTE THE REDSHIFT(S).
       zans = zcompute(objflux, objivar, star, poffset=poffset, $
                  pmin=pmin, pmax=pmax, nfind=nfind, width=width, $
                  plottitle=plottitle, _EXTRA=EXTRA)

   ;----------
   ; Convert redshift (and error) from pixels to the conventional dimensionless
   ; value.  Do not modify any errors that are less than zero, since those
   ; can be used as just warning flags from the fit.

      indx = where(zans.dof GT 0, npeak)
      if (npeak GT 0) then $
        zans[indx].z = 10.^(objdloglam * zans[indx].z) - 1.

      jndx = where(zans.dof GT 0 and zans.z_err GE 0)
      if (jndx[0] NE -1) then $
         zans[jndx].z_err = $
         alog(10d) * objdloglam * zans[jndx].z_err * (1 + zans[jndx].z)

   ;----------
   ; Copy valid peaks into the output structure
      nobj = 1
      result = replicate({zresult}, nfind, nobj)
      if (npeak GT 0) then begin
         result[indx].z = zans[indx].z
         result[indx].z_err = zans[indx].z_err
         result[indx].rchi2 = zans[indx].chi2 / (zans[indx].dof > 1)
         result[indx].dof = zans[indx].dof
         ntheta = n_elements(zans[0].theta)
         result[indx].theta[0:ntheta-1] = reform(zans[indx].theta)
         result[indx].tfile = fileandpath(thisfile)
         result[indx].tcolumn[0] = istar ;get appropriate column from template
         result.npoly = npoly
         result.class = 'STAR'
         result.subclass = subclass[istar] ;type of sta
;print, 'finished with stellar type: ', subclass[istar]
       endif
       
       if istar eq 0 then tresult = result $
         else tresult = [tresult, result] ;concatentate outputs 
   endfor ;end loop over stellar types
   
   tresult = tresult[sort(tresult.rchi2)] ;sort by rchi2
   return, tresult[0:4] ;return 5 best entries
end
;------------------------------------------------------------------------------















