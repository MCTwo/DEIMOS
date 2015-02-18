;+
; NAME:
;   makewindowf
;
; PURPOSE:
;   create a window function of regions observed, based on mask coords.
;
; CATEGORY:
;   correlation analysis
;
; CALLING SEQUENCE:
;   windowf =makewindowf,  nra,ndec, masks=masks, ralim=ralim, $
;               declim=declim, nozcorr=nozcorr
;
; INPUTS:
;   nra -- number of pixels in RA direction (default =3000)
;   ndec -- number of pixels in dec direction (default=1000)
;   masks--     list of masks data structures
;
; KEYWORDS:
;   ralim  -- if set, limits on RA, otherwise gleaned roughly from central mask
;   declim -- if set, limits on dec
;   nozcorr -- if set, do NOT apply the z-completeness correction per mask 
;
; OUTPUTS:
;   windowf -- window function
;
; MODIFICATION HISTORY:
;  md 09dec02
;-


pro mask_ccr,  x, y,  ra,  dec, rac, decc, pa, bow,  height
; rotate from mask coords. into ra,dec (degrees)
  sinpa = sind(pa) &  cospa= cosd(pa)
  x = x -bow*((2.*y/height)^2-.5)  ;maximum bow at y=height/2. and middle 
  ra =  (x*cospa + y*sinpa)/3600./cosd(decc) + rac
  dec = (-x*sinpa + y*cospa)/3600. + decc

  return
end

function makewindowf,  nra, ndec, masks=masks,  ralim=ralim, declim=declim, nozcorr = nozcorr

   maskx = [0., .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, 0.0]
   masky = [-.5, -.5, -.4, -.3, -.2, -.1, 0., .1,  .2,  .3, .4, .5, .5]
   maskx = [maskx, -maskx] ;concatenate array
   masky = [masky, reverse(masky)] ;concatenate again
   maskxf = [-0.5,  0.5, 0.5,-0.5]
   maskyf = [-0.5, -0.5, 0.5, 0.5]
   sfraction = 0.65

   if n_params() eq 0 then begin
     nra = 3000
     ndec =  1000
   endif
   if n_params() eq 1 then ndec = 1000 ;default array sizes

   if NOT keyword_set(ralim) then  ralim = $ ;degrees
            [min(masks.ra)-.1, max(masks.ra)+.1]
   if NOT keyword_set(declim) then declim = $ ;degrees
            [min(masks.dec)-.14, max(masks.dec)+.14]
   ra_inc = (ralim[1]-ralim[0])/nra
   dec_inc = (declim[1]-declim[0])/ndec

   
   tempwindowf = dblarr(nra, ndec) ;declare window array

; read in redshift completeness
   zfile='$DEEP2PRODUCTS/zspec_archive/zcomp.dat'
;   zfile='$DEEP2PRODUCTS/zcomp.dat'
   readcol, zfile, maskname, zcomp, zcompblue, zcompred, format = 'I,F,F,F'


   Y_ang_offset = ((187. + masks[0].width +187.)/2 -4.5*60.)/3600.
 ;offset of ra1 from pointing center, in arcsec, converted to degrees
;start first pass on selection
   nmasks = n_elements(masks)

   windowf = 1.d0

   for i=0,nmasks-1  do begin ;plot boundaries of masks

     tempwindowf = tempwindowf*0.d0

     rac  = masks[i].ra - Y_ang_offset/cosd(masks[i].dec)*cosd(masks[i].pa) 
     decc = masks[i].dec+ Y_ang_offset*sind(masks[i].pa)

; get redshift completeness for this mask
    thismask = where(maskname eq masks[i].masknumber, ct)

; if mask wasn't observed use z completeness of 100% (this allows us
     ; to make a cut in chi-squared instead of only using galaxies which
     ; have been pspec'd
    if ct gt 0 then zfactor = zcompred[thismask] else zfactor = 1.

; get mask coords core region
    mask_ccr,  maskx*masks[i].cwidth, masky*masks[i].height , ra_mask, $
        dec_mask, rac, decc, masks[i].pa, masks[i].bow, $ ;include bow
        masks[i].height 

    xwindow = (ra_mask-ralim[0])/ra_inc ;get mask x,y
    ywindow = (dec_mask-declim[0])/dec_inc
    jcore = polyfillv(xwindow, ywindow, nra, ndec)
    if keyword_set(nozcorr) then tempwindowf[jcore] = tempwindowf[jcore] + $
        0.75*sfraction else tempwindowf[jcore] = tempwindowf[jcore] + $
        0.75*sfraction*zfactor[0]    ;set to 3/4 full value

; get mask coords, outer region
    mask_ccr,  maskxf*masks[i].width, maskyf*masks[i].height , ra_mask, $
        dec_mask, rac, decc, masks[i].pa, 0. , masks[i].height ;no bow
; get mask coords outer region
    xwindow = (ra_mask-ralim[0])/ra_inc
    ywindow = (dec_mask-declim[0])/dec_inc
    j2 = polyfillv(xwindow, ywindow, nra, ndec)
    if keyword_set(nozcorr) then tempwindowf[j2] = tempwindowf[j2] + $
     .25*sfraction else tempwindowf[j2] = tempwindowf[j2] + .25*sfraction*zfactor[0]  
;given another chance- adds to core as well giving the core 4* the
     ;sel. probability as the edge, as observed in real masks

    windowf = windowf*(1.d0- (tempwindowf < 1)) ;deal with probabilities properly

   endfor

   windowf =  1.d0 -windowf;limit overlapping zones

   whlow = where(abs(windowf) lt 1.e-6, lowct) ; just to make sure we don't have numerical problems
   if lowct gt 0 then windowf[whlow] = 0.
return, float(windowf)
end





