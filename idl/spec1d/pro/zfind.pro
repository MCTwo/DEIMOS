;+
; NAME:
;   zfind
;
; PURPOSE:
;   Find possible redshift matches for a set of spectra using a set of
;   eigen-templates.
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
;-----------------------------------------------------------------------------

function findpix, lambda, lambda0
  diff = abs(lambda - lambda0)
  minval = min(diff, minsub)
  minsub = (minsub - 1) > 0
  return, minsub
end

function zfind, ss1d, eigenfile=eigenfile, eigendir=eigendir, $
                columns=columns, npoly=npoly, zmin=zmin, zmax=zmax, $
                zguess=zguess, pwidth=pwidth, nfind=nfind, width=width, $
                linear_lambda=linear_lambda,objflux=objflux, $
                objivar=objivar, loglam=loglam, $
                nosubtract=nosubtract, wvmin=wvmin, wvmax=wvmax, $
                _EXTRA=EXTRA

  common com_zfind, starflux, starloglam0, stardloglam, $
    nstars, npoints, starflux_corrected, thisfile


; check if continuum subtraction is required.
  if n_elements(nosubtract) gt 0 then $
    nosubtract = nosubtract[0] ge 1 $
  else nosubtract = 0

  if (NOT keyword_set(starflux)) then begin ;if first call, read template

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
;  splog, 'Selecting EIGENFILE=' + thisfile

    starflux = readfits(thisfile, shdr, /SILENT)
;;; FROM THE HEADER, GRAB THE WAVELENGTH INFORMATION.
    starloglam0 = sxpar(shdr, 'COEFF0')
    stardloglam = double(sxpar(shdr, 'COEFF1'))
    npoints = n_elements(starflux[*, 0])
    starorder =  round(3.*npoints/7000.) ;3 terms/7000 points 
    starloglam = findgen(npoints)*stardloglam + starloglam0

; remove continuum from template at this level
    nstars = (size(starflux, /dimen))[1] ;how many templates?

;TEST!!  do subtraction!!
    if nosubtract eq 0 then begin
      for i=0, nstars-1 do begin 
;       starparams = svdfit(starloglam, starflux[*, i], starorder, $
;          /double,  /legendre,  yfit=starcont)
        if i eq 1 then continue ; skip the emission line template.
        starcont = djs_median(starflux[*, i], width=2500, boundary='reflect')
        starflux[*, i] = starflux[*, i] - float(starcont)

      endfor
    endif
    starflux_corrected = starflux ;save for next call
    
  endif else starflux = starflux_corrected ;if not 1st call, restore template.

;;; CHECK IF THE columns ARGUMENT WAS PASSED. IF SO, TRIM THE TEMPLATE
;;; ACCORDINGLY ONLY KEEPING THE DESIRED TEMPLATE SPECTRA. 
  IF n_elements(columns) NE 0 THEN BEGIN
    starflux = starflux[*,columns]
  ENDIF ELSE BEGIN
    columns = LINDGEN(nstars)
  ENDELSE

;;; IF THE OPTIONAL PARAMETER npoly IS PASSED, THEN ADD THE
;;; APPROPRIATE NUMBER OF POLYNOMIAL TERMS TO THE TEMPLATE ARRAY.
  IF KEYWORD_SET(npoly) THEN $
    starflux = [ [starflux], [poly_array(npoints,npoly)] ]

; DETERMINE GRID SIZE (IN LOG LAMBDA) FOR THE OBJECT SPECTRUM AND
; DETERMINE THE WAVELENGTH RANGE FOR THE OBJECT. Force binning to 
; match that of template 
  IF KEYWORD_SET(linear_lambda) THEN $
    loglam = linear2log(ss1d, binsize=stardloglam, flux=objflux, ivar=objivar)
  objloglam0 = loglam[0]
  objdloglam = stardloglam; loglam[1] - loglam[0] 
;should be same as stardloglam, if code works
;  print, ' log lam range: ', loglam[0], loglam[n_elements(loglam)-1]

;;; NOW SMOOTH THE OBJECT SPECTRUM AND REMOVE THE CONTINUUM.
  nobj = n_elements(objflux)
  objorder = 3 ;coeff. of polynomial to subtract

; TEST!! try do continuum, lower order. 
;  objparams = svdfit(loglam, objflux, objorder, /legendre, $
;       /double, yfit=objcont)
  if nosubtract eq 0 then begin
      objcont = djs_median(objflux, width=2500, boundary='reflect')
      objflux = objflux - float(objcont)
  endif

; check if we need to trim the input spectrum to just a subregion.
  if n_elements(wvmin) gt 0 then wvmin = wvmin[0] else wvmin = -1
  if n_elements(wvmax) gt 0 then wvmax = wvmax[0] else wvmax = -1
  if wvmax ge 0 and wvmin ge 0 then begin
      minpix = findpix(loglam, alog10(wvmin))
      maxpix = findpix(loglam, alog10(wvmax))
;      ss1d = trimss(ss1d, minpix, maxpix)
      loglam = loglam[minpix:maxpix]
      objflux = objflux[minpix:maxpix]
      objivar = objivar[minpix:maxpix]
      objloglam0 = loglam[0]
  endif

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
; old version - parentheses error?
;  maxp = fix(objloglam0 + objdloglam*.99*n_elements(objflux) - starloglam0/ $
;          objdloglam)  


;ALTERATION BY BJW, 8/21/03 
 maxp =long((objloglam0 + objdloglam*.99*n_elements(objflux) - starloglam0)/$
            objdloglam)  

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
;;; CALL zcompute.pro TO COMPUTE THE REDSHIFT(S).
   zans = zcompute(objflux, objivar, starflux, poffset=poffset, $
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
      for icol=0, n_elements(columns)-1 do $
       result[indx].tcolumn[icol] = columns[icol]
      result.npoly = npoly
   endif

   return, result
end
;------------------------------------------------------------------------------















