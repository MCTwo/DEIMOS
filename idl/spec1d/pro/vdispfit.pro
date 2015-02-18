;+
; NAME:
;   vdispfit
;
; PURPOSE:
;   Compute velocity dispersions for galaxy spectra.
;
; CALLING SEQUENCE:
;   vdispfit, objflux, objivar, [ objloglam, hdr=, zobj=, npoly=, $
;    sigma=, sigerr= ]
;
; INPUTS:
;   objflux    - Galaxy spectrum (spectra); array of [NPIX,NGALAXY].
;   objivar    - Galaxy inverse variance; array of [NPIX,NGALAXY].
;
; OPTIONAL INPUTS:
;   objloglam  - Log-10 wavelengths; this can be either an NPIX vector
;                if all the galaxy spectra have the same wavelength mapping,
;                or an array with the same dimensions as OBJFLUX.
;                Either OBJLOGLAM or HDR must be specified.
;   hdr        - FITS header from which to read COEFF0, COEFF1 for the
;                wavelength mapping.
;                Either OBJLOGLAM or HDR must be specified.
;   zobj       - Redshift for each galaxy; default to 0.
;   npoly      - Number of polynomial terms to append to eigenspectra;
;                default to 5.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   sigma      - Velocity dispersion in km/sec.
;   sigerr     - Error for SIGMA in km/sec.
;
; COMMENTS:
;   Note that the wavelength spacing in the galaxy and stellar template spectra
;   must be the same.
;   NOTE: this code uses the emission line template and is not intended for
;   abs. line dispersions.

; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC1D_DIR/templates/TEMPLATEFILES
;
; PROCEDURES CALLED:
;   airtovac
;   combine1fiber
;   computechi2()
;   djs_filepath()
;   findchi2min
;   mrdfits()
;   poly_array()
;   splog
;   sxpar()
;
; INTERNAL SUPPORT ROUTINES:
;   vdisp_gconv()
;
; REVISION HISTORY:
;   13-Mar-2001  Written by D. Schlegel, Princeton
;   01-Oct-2002  Revised by MD
;   05-Nov-2002  Revised by mcc
;------------------------------------------------------------------------------
function vdisp_gconv, x, sigma, _EXTRA=EXTRA

   ; Special case for no smoothing
   if (sigma EQ 0) then return, x

   ksize = round(4*sigma+1) * 2
   xx = findgen(ksize) - ksize/2

   kernel = exp(-xx^2 / (2*sigma^2))
   kernel = kernel / total(kernel)

   return, convol(x, kernel, _EXTRA=EXTRA)
end

;------------------------------------------------------------------------------
pro vdispfit, ss1d, zobj=zobj, npoly=npoly, sigma=sigma, $
              sigerr=sigerr, doplot=doplot, silent=silent

   common com_vdispfit, bigflux, bigloglam, bigmask, nsamp, bigsig, $
    nbigpix, nsig, dsig,  edloglam

;;; CHECK IF KEYWORDS npoly OR zobj WERE PASSED. IF NOT, SET TO
;;; DEFAULT VALUES (npoly = 0 AND zobj = 0.0).
   if (n_elements(npoly) EQ 0) then npoly = 0
   if keyword_set(zobj) then zobj = zobj[0] ELSE zobj = 0.0
   if n_elements(doplot) gt 0 then doplot = doplot[0] ge 1 $
   else doplot = 0
   if n_elements(silent) gt 0 then silent = silent[0] ge 1 $
   else silent = 0

   ;---------------------------------------------------------------------------
   ; Generate the over-sampled eigen-templates for the stellar spectra.
   ; This is saved in a common block between calls.

   if (NOT keyword_set(bigflux)) then begin
      nsamp = 1 ; (was 10)
      nsig = 40 ;was 25
      dsig = 10. ;was 25.
      bigsig = findgen(nsig) * dsig ; in km/sec

      ;----------
      ; Find the most recent template file matching EIGENFILE
      eigenfile = 'Vdisp*.fits'
      eigendir = concat_dir(getenv('IDLSPEC1D_DIR'), 'templates')
      allfiles = findfile(djs_filepath(eigenfile, root_dir=eigendir), count=ct)
      if (ct EQ 0) then $
       message, 'Unable to find EIGENFILE matching '+eigenfile
      thisfile = allfiles[ (reverse(sort(allfiles)))[0] ]
      if not(silent) then splog, 'Selecting EIGENFILE=' + thisfile

      ;----------
      ; Read the dispersion templates
      eflux = mrdfits(thisfile, 0, ehdr, /silent)
      naxis1 = sxpar(ehdr,'NAXIS1')
      nstar = sxpar(ehdr,'NAXIS2') > 1
      eloglam0 = sxpar(ehdr, 'COEFF0')
      edloglam = double(sxpar(ehdr, 'COEFF1'))
      eloglam = eloglam0 + dindgen(naxis1) * edloglam

      ; Pixel size in km/sec for these oversampled (smaller) pixels
      cspeed = 2.99792458e5
      pixsz = (10.^(edloglam)-1) * cspeed / nsamp

      ;----------
      ; Re-sample the templates to higher resolution by a factor of NSAMP
; don't do the resample for DEIMOS
      nbigpix = (naxis1 - 1) * nsamp + 1
      bigloglam = eloglam0 + dindgen(nbigpix) * edloglam / nsamp
      bigflux = fltarr(nbigpix, nstar, nsig)

;      for istar=0, nstar-1 do begin
;         combine1fiber, eloglam, eflux[*,istar], $
;          newloglam=bigloglam, newflux=tmpflux, maxiter=0
         bigflux[*,*,0] = eflux ;will work for multiple stars
;         if (istar EQ 0) then bigmask = tmpflux NE 0
;      endfor

      ;----------
      ; Generate array of broadened templates
      for isig=1, nsig-1 do begin
         for istar=0, nstar-1 do begin
            bigflux[*,istar,isig] = $
             vdisp_gconv(bigflux[*,istar,0], bigsig[isig]/pixsz, $
                   /edge_truncate)
         endfor
      endfor
      bigmask = bigflux[*, 0, nsig-1] ne 0 
;mask keeps only regions with signal in template 
;this is a simple means of cutting out non-emission regions for DEEP2
   endif

   ;----------
   ; Determine the wavelength mapping for the object spectra,
   ; which are the same for all of them.
;  TRANSFORM THE LINEAR LAMBDA VALUES IN THE STRUCTURE ss1d INTO
; LOG(LAMBDA) BINNED VALUES.
   objloglam = linear2log(ss1d, binsize=edloglam, flux=objflux, ivar=objivar)
   objloglam0 = objloglam[0]
   objdloglam = objloglam[1] - objloglam[0]
; DE-REDSHIFT THE OBJECT SPECTRUM.
   restloglam = objloglam - alog10(1 + zobj) 
; define npixobj, the number of pixels in the object spectrum.
   npixobj = N_ELEMENTS(objflux)


  ;----------
   ; Find the pixel numbers to use from the object and the templates

   ; Find the sub-pixel shifts in the object
;   subshift = round(((bigloglam[0]-restloglam[0]) / objdloglam MOD 1) * nsamp)
   ; Bug fix since the IDL MOD function does the wrong thing for negatives.
   subdum = (bigloglam[0]-restloglam[0]) / edloglam
   subshift = round((subdum-floor(subdum)) * nsamp)
   indx = subshift + nsamp * lindgen(nbigpix/nsamp)

   if (max(restloglam) LT min(bigloglam[indx]) $
    OR min(restloglam) GT max(bigloglam[indx])) then begin
;      if not(silent) then splog, 'No wavelength overlap with template'
      sigma = 0.0
      sigerr = 9999.
      return
   endif

   if (restloglam[0] LT bigloglam[indx[0]]) then begin
      ipixt0 = 0L
      junk = min(abs(restloglam - bigloglam[indx[0]]), ipixo0)
   endif else begin
      ipixo0 = 0L
      junk = min(abs(bigloglam[indx] - restloglam[0]), ipixt0)
   endelse

   npixcomp = (npixobj - ipixo0 ) < (n_elements(indx) - ipixt0)

   indxo = ipixo0 + lindgen(npixcomp) ; Indices for object spectrum
   indxt = indx[ipixt0 + lindgen(npixcomp)] ; Indices for template spectra

   ;----------
   ; Add more eigen-templates that represent polynomial terms.

   if (keyword_set(npoly)) then $
    polyflux = poly_array(npixcomp,npoly)

   ;----------
   ; Fit for chi^2 at each possible velocity dispersion

   chi2arr = fltarr(nsig)
   yfit = fltarr(nsig,n_elements(indxt))

   objsmall = objflux[indxo]
   sqivar = sqrt( objivar[indxo] > 0 ) * bigmask[indxt]

   for isig=0, nsig-1 do begin
      eigenflux = bigflux[indxt,*,isig]
      if (keyword_set(npoly)) then eigenflux = [[eigenflux], [polyflux]]
      chi2arr[isig] = computechi2(objsmall, sqivar, eigenflux, yfit=yfit1)
      yfit[isig,*] = yfit1
   endfor

   ;----------
   ; Fit for the dispersion value at the minimum in chi^2

   findchi2min, bigsig, chi2arr, minchi2, sigma, sigerr

   if doplot then begin
       plot, restloglam[indxo], objsmall, xr=[3.565, 3.58], $
         xtitle='log(lambda)', title='vdispfit', xthick=1.5, ythick=1.5
       minchi = min(chi2arr, minsub)
       oplot, bigloglam[indxt], yfit[minsub,*], thick=2, color=2, line=2
   endif
   ;----------
   ; If the best-fit value is at the maximum dispersion value tested,
   ; then we don't really know the answer and should set the error
   ; to a large value.

   if (sigma GE max(bigsig)) then sigerr = 9999.
; check that sigerr variable is defined. if not, then set it to big
; value. 
   if not(keyword_set(sigerr)) then sigerr = 9999.

;print, 'sigma,sigerr: ', sigma, sigerr

   return
end
;------------------------------------------------------------------------------










