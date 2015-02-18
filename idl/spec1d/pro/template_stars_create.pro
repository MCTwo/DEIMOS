;+
; NAME:
;    template_stars_create
;
; PURPOSE:
;    generate a version of SDSS stellar template at higher resolution
; 
; CALLING SEQUENCE:
;    template_stars_create,filename
;
; INPUTS:
;    filename  -- name of output file
;
; COMMENTS:
;    this procedure uses the high SN SDSS template of many stars,
;    simply changes the dloglam from .0001 to 2.e-5, to match the
;    other template
;
; REVISION HISTORY:
;   MD 18nov02
;----------------------------------------------------------------------
pro template_stars_create, fileout

  eigendir = concat_dir( GETENV('IDLSPEC1D_DIR'), 'templates' )
  eigenfile = 'spEigenStar.fits'
  sdss_stars = readfits(eigendir +'/' + eigenfile, hdr)

  c0 = sxpar(hdr, 'COEFF0')
  c1 = sxpar(hdr, 'COEFF1')
  
  nspec =  (size(sdss_stars, /dimen))[0]
  nstars = (size(sdss_stars, /dimen))[1]


  logl = findgen(nspec)*c1 +c0
  lambda0 = 10.^c0 ;initial wavelength  
  lambda1 = 10.^logl[nspec-1] ;final wavelength
  dloglambda = 2.e-5
  fxaddpar, hdr, 'COEFF1', dloglambda ;update FITS keyword

; CHANGED BY JAN 7/31/03
; Templates are sometimes missing outside 3788-9253 AA
;  and sometimes noisy outside 3820-9192 AA.  Before, we just
;  interpolated in those regions (from 0).  Now, we should only use
;  the template in the good regions.
  whtemplates=where(logl ge alog10(3820.) AND logl le alog10(9192.),ngood)
  newc0=logl[whtemplates[0]]
  fxaddpar, hdr, 'COEFF0', newc0 ;update FITS keyword
  specnum=(logl[whtemplates[ngood-1]]-logl[whtemplates[0]])/dloglambda

;  specnum = nspec*c1/dloglambda ;number of elements in stellar spectrum
  stars_deep = fltarr(specnum, nstars)

  biglogl = newc0 + findgen(specnum)*dloglambda

  for i=0, nstars-1 do begin
    good = where(sdss_stars[*, i] gt 0.) ;don't include dropouts

; the following caused trouble (as we interpolated from 0's:)
;    good = [0, good+1, nspec-1] ;append the first and last points to the good
; and offset pointer array by one!

    sspline = spl_init(logl[good], sdss_stars[good, i])
; spline interpolates over dropout points

    star_refined = spl_interp(logl[good], sdss_stars[good, i], sspline, $
         biglogl)
    stars_deep[*, i]= star_refined/mean(star_refined) ;normalize
  endfor

  writefits, fileout, stars_deep, hdr


; write trimmed version of file
  if findfile('spEigenStarDeep2.fits') eq 'spEigenStarDeep2.fits' then begin

    naddlow = 5000
    naddhigh = 4000
    newc0 = newc0-naddlow*dloglambda


      hdr2=headfits('spEigenStarDeep2.fits')
      fxaddpar, hdr2, 'COEFF0', newc0 ;update FITS keyword
      fxaddpar, hdr2, 'COEFF1', dloglambda ;update FITS keyword
  
      whnames=where(strpos(hdr,'NAME') eq 0,nnames)
      whnames2=where(strpos(hdr2,'NAME') eq 0,nnames2)

      startnames=min(whnames)
      startnames2=min(whnames2)
      value=hdr
      value2=hdr2

; need to match up stellar types in the two files
      for i=0,n_elements(hdr)-2 do $
        value[i]=(strsplit(hdr[i],'=',/extract))[1]
      for i=0,n_elements(hdr2)-2 do $
        value2[i]=(strsplit(hdr2[i],'=',/extract))[1]

      noldspec = n_elements(stars_deep(*, 0))
      nspecpix = noldspec+naddlow+naddhigh
      stars_deep2 = fltarr(nspecpix, n_elements(whnames2))

      for i=0,nnames2-1 do begin
         idx=where(value[whnames] eq value2[whnames2[i]])
; do constant extensions of the template at either end.  It is not
; clear that this is the optimum solution, but it is at least a
; solution to the problem of stellat templates not covering our full
; wavelength range.

         stars_deep2[0:naddlow-1,i]= $
           median(stars_deep[0:10, idx])
         stars_deep2[noldspec+naddlow:nspecpix-1,i]= $
           median(stars_deep[noldspec-11:noldspec-1, idx])
         stars_deep2[naddlow:noldspec+naddlow-1,i]=stars_deep[*,idx]

        
      endfor


      fitspos=strpos(fileout,'.fits')
      fileout2=strmid(fileout,0,fitspos)+'2.fits'

      writefits, fileout2, stars_deep2, hdr2

  endif

end


