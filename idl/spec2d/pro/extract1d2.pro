;+
;
; NAME
;      extract1d2.pro
;
; PURPOSE
;      To extract a 1d spectrum from a sky-subtracted 2-d slit file.  
;      The function extract1d takes a slitfile,position, a full width
;      half max,a weight, inverse variance weight,a shift type,
;      and extraction type as its required arguments and returns a 
;      structure which contains the 1-dimensional spectrum for the given slit. 
;      The structure contains the flux, wavelength, and inverse 
;      variance along the slit. The 1-d spectrum is extracted according to a
;      user-specified algorithm (either a weighted boxsprof 
;      extraction or a horne extraction
;      algorithm) with presently the boxsprof extraction set
;      as the default extraction setting. 
;
;      This routine differs from the old DEEP2/DEIMOS EXTRACT1D
;      routine in that it:
;
;      1) by default interpolates, rather than shifts, pixels in the X
;      (wavelength) direction
;      2) Accounts for the effects of atmospheric dispersion in
;      defining the weight/extraction region: e.g. boxcar extraction
;      regions are now curved, rather than the same at all wavelengths
;      3) is written in a more general-purpose fashion.  Other
;      weighting algorithms can be added straightforwardly



;
; SYNTAX
;      ss1d = extract1d2(slitfile, pos, fwhm [hdr=hdr,
;                        weight = weight /optimal, /boxcar, /horne, $
;                        /boxsprof,/nonlocal,/ivarwt, /shift,$
;                        flux = flux  nsigma=nsigma])
;
; INPUTS
;      slitfile = a structure containing the sky-subtracted,
;                 cosmic-ray rejected 2-d flux values, the wavelength
;                 solution, and the inverse variance
;                 information. Also, recall that the 2-d spectrum from
;                 a given slit is broken into two slit files: red
;                 chip/blue chip. That is, a single slitfile only
;                 contains the blue portion of the spectrum or the red
;                 end of the spectrum. Can also be a file containing
;                 the structure (a slit file).
;      pos = the position of the object in the slit. pos gives the
;            pixel column about which to extract the object.
;      fwhm = a parameter specifying the width in the spatial
;             direction over which to extract the 1d spectrum. In the
;             optimal extraction algorithm this specifies the
;             full-width at half-maximum (fwhm) for the Gaussian
;             profile.
;      hdr = this optional parameter is provided so that if the user
;            passes a slit structure (rather than file), the header
;            for the slit FITS file can be passed to the routine. This
;            header is employed to locate the slitno and then the
;            corresponding spSlit file.
;
; KEYWORDS
;      optimal = if this keyword is set, then the spectrum is
;                extracted using the horne extraction
;                algorithm. Recall that if the /optimal keyword is not
;                set, then the spectrum is extracted according to a boxsprof
;                inverse variance weighted technique meaning that
;                each pixel within a given spatial width is weighted
;                according to its inverse variance.
;
;      horne = if this keyword is set, then the spectrum is extracted
;              using horne extraction algorithm as detailed by
;              K. Horne (1986, PASP, 98, 609). Note that in this 
;              approach we assume that the spatial profile is
;              described by a Gaussian, but weight by inverse variance. 
;      boxsprof = this is the default setting, then the spectrum is
;                 extracted using an un-weighted boxcar
;                 extraction algorithm. That is, outside a
;                 given width the weight is 0 and inside the defined
;                 spatial range the weight is given by the inverse
;                 variance. This algorithm differs from
;                 the /boxcar extraction in how it handles
;                 bad-pixels. Here bad-pixels are masked and the
;                 fraction of flux
;                 missed within these pixels is determined using the
;                 object's spatial profile. The final spectrum is
;                 scaled so as to account for this missed flux.
;      boxcar =  If this keyword is used it will be set as a boxsprof 
;      nonlocal = if this keyword is set, then the extracted 1-d
;                 spectrum is extracted from the non-local
;                 sky-subtracted 2-d slit spectrum.
;
; OUTPUTS
;      ss1d = a structure containing the 1-d extracted spectrum:
;             spec.spec = a vector containing the flux values as a
;                         function of wavelength. Note that the
;                         spec.flux array has curious dimensions. 
;             spec.ivar = a vector containing the inverse variance of
;                         the extracted 1-d spectrum. 
;             spec.lambda = a vector containfng the wavelength values
;                           corresponding to the the flux and inverse
;                           variance values given in spec.flux and
;                           spec.ivar.
;             spec.flux = a 2-D array of the intensity at each pixel 
;             spec.weight = a 2-D array giving the weight assigned to
;                           each pixel 
;             spec.ivarfudge = an empirical factor by which to
;             multiply spec.ivar to make it match observed
;             fluctuations (in a median sense).  It is saved for
;             diagnostic purposes. 
;
; PROCEDURES CALLED 
;      lambda_eval
;      find_object
;      mcc_polyfit

; estimate how far off our ivar estimate is, by comparing deviations
; from median-smooth version of spectrum to expected deviations
function ivarfudge, spectrum, ivar

  nspec=n_elements(spectrum)
  nbins=8

  fudgearr=fltarr(nbins)
  binarr=fudgearr

  for i=0,nbins-1 do begin
      whbin = where(lindgen(nspec) ge (nspec/8.*i>50) and $
                    lindgen(nspec) le  (nspec/8.*(i+1)<(nspec-50)) and $
                    ivar ne 0. and spectrum ne 0.,binct)
      binarr[i]=binct
      if binct gt 100 then begin
          medlevel = median(spectrum,50)
          fudgearr[i] = djsig((spectrum[whbin]-medlevel[whbin]) * $
                              sqrt(ivar[whbin]))
      endif else fudgearr[i] = 0.
  endfor

  whok=where(binarr gt 200,okct)

  if okct gt 0 then $
    fudgefactor=1/(median(fudgearr[whok],/even))^2 $
  else fudgefactor=1.
;  if okct gt 0 then print,median(fudgearr[whok],/even)

  return,fudgefactor
end



function extract1d2, slitfile, pos, fwhm, shift=shift, horne = horne, $
   flux=flux,weight = weight,  nsigma=nsigma, optimal=optimal,boxcar = boxcar,$
    boxsprof=boxsprof, nonlocal=nonlocal, hdr=hdr,ivarwt=ivarwt


; check the various keyword settings.
;  if n_elements(optimal) gt 0 then $
;    optimal = optimal[0] ge 1 else optimal = 0
   if n_elements(optimal) gt 0 then if optimal gt 0 then horne = 1 $
    else begin
     horne=optimal
     optimal=0
    endelse


 ; if n_elements(horne) gt 0 then $
   ; horne = horne[0] ge 1 else horne = 0
  if n_elements(horne) eq 0 then horne=0 $
     else horne = horne[0] gt 0
 ; boxcar is always boxsprof
  if n_elements(boxcar) gt 0 then if boxcar gt 0 then boxsprof = 1 

  if n_elements(boxsprof) gt 0 then $
    boxsprof = boxsprof[0] ge 1 else boxsprof = 1-horne
  if n_elements(nonlocal) gt 0 then $
    nonlocal = nonlocal[0] ge 1 else nonlocal = 0

  if n_elements(shift) eq 0 then shift = 2
  if n_elements(ivarwt) eq 0 then ivarwt=0
 
if horne gt 0 and boxsprof gt 0 then horne=0


if keyword_set(nsigma) then nsigma=nsigma[0] else nsigma=1.5

; check if the argument slitfile is a FITS file or an IDL structure
; and then define slit accordingly.
  if size(slitfile, /tname) eq 'STRING' then begin
     if keyword_set(nonlocal) then begin
         fits_info, slitfile, /silent, n_ext=n_ext
         if n_ext gt 2 then begin
             slit = mrdfits(slitfile, 3, hdr, /silent)
          endif else return, 0
      endif else slit = mrdfits(slitfile, 1, hdr, /silent)
      slitno = strcompress(sxpar(hdr, 'SLITNO'), /rem)
  endif else begin
      slit = slitfile
      if keyword_set(hdr) then $
          slitno = strcompress(sxpar(hdr, 'SLITNO'), /rem) $
      else slitno = ''
  endelse

; construct an empty spectrum: filled with zeros. this will be
; returned as the result in case of various errors in the extraction.
  n = n_elements(slit.flux[*,0])
  nrows=n_elements(slit.flux[0,*])

  skyspec=fltarr(n)
  usesky=0
  nsky=0
  readnoise = 2.32


; JAN first sky code here
  if slitno ne '' then begin
; find the spSlit file corresponding to the given slit file/structure.
      endpos = strpos(slitfile, '.fits')
      spname = 'spSlit' + strmid(slitfile, 4, endpos-4) + '.fits*'
      spslitfile = findfile(spname, count=numsp)
      if numsp gt 0 then begin
          spslitfile = spslitfile[0]
          fits_info, spslitfile, /silent, n_ext=n_ext
          nbsplines=(n_ext / 2)
          usesky=1
          skyarr=fltarr(n,nbsplines)
          spslit=mrdfits(spslitfile, 1, /silent)
                                ; keep track of # of good sky pixels
                                ; in each column.  this should roughly
                                ; correspond to the extracted
                                ; pixels... INCORPORATE CRMASK!
          skydex = where(spslit.skyrow, skydexcnt)
          if skydexcnt gt 0 then $
            nsky=total((spslit.mask[*,where(spslit.skyrow)] AND 22b) eq 0b,2) $
          else begin
              print, 'No skyrows in slit!'
              nsky = 0
          endelse
      endif
      exptime=0.
  endif


; define a 1-d spectrum structure to return in case of error.
  zero1d = {spec:fltarr(n), lambda:fltarr(n), ivar:fltarr(n), $
            crmask:intarr(n), bitmask:intarr(n), ormask:intarr(n), $
            nbadpix:intarr(n), infomask:intarr(n), $
            objpos:float(pos), fwhm:float(fwhm), $
            nsigma:0.0, r1:0, r2:0, skyspec:fltarr(n), ivarfudge:1.} 

; get the 2-d lambda array, takes into account the three diff types of lambda
  lambda2d = lambda_eval(slit)
; extract the wavelength values along the center of the object.
  cwave = lambda2d[*,pos]
  npix = n_elements(cwave)
; the dispersion level (angstrom/pixel) (wavelength/pixel), dlambda/dx
  dldx = (cwave[npix-1] - cwave[0]) / npix
  ;this is dlambda/dy
  tiltx = (lambda2d[npix/2,nrows-1] - lambda2d[npix/2,0]) / nrows
  dxdp = tiltx/dldx ;slope of constant lambda/vert pix, dx/dy



; build up sky spectrum, for use in determining its inverse variance
; JAN - sky ivar code here
  if usesky gt 0 then begin
      for i=0,nbsplines-1 do begin
          bspline=mrdfits(spslitfile,i*2+2, /silent)
          skyarr[*,i]=bspline_valu(cwave,bspline)
          exptime=exptime + $
            sxpar(headfits(spslitfile,ext=i*2+1, /silent),'EXPTIME')
  endfor

; assume average spec. is ok for statistics, and that don't have to
; worry about bad columns, etc. in sky region.  
      avgexptime=exptime/nbsplines
      if nbsplines gt 1 then skyspec=total(skyarr,2)/nbsplines $
      else skyspec = skyarr
; note the bspline spectrum is an average, not a sum, of the values at
; each pixel that went into it.
      skyctsper2dpixperspec=skyspec*(avgexptime/3600.)
      skyivar= nsky*nbsplines/(skyctsper2dpixperspec+readnoise^2)
      skyivar=skyivar*(avgexptime/3600.)^2

; maybe add...
;      skyivar=skyivar*(2./3.) ; as typically are using ~nsky*2/3 pixel
        ; to define the sky value at a breakpoint. only slightly kludgey.

; this is the inverse-variance in the sky spectrum at each pixel of
; the extraction.  Note that this contributes, wholly covariantly, to
; EVERY pixel at a given wavelength.
  endif else begin
      print, 'No sky info found!!!'
      skyivar = fltarr(n) + 1.0
  endelse





if horne eq 1 then begin

   ; define range in the spatial direction [r1:r2] over which to do
   ; extraction. exclude the first and last few pixels due to slit edge
   ; effects.
     r1 = 4
     r2 = nrows-5

     width = fwhm * 2
     ext = (ceil(width/2.)) < pos < (nrows-1-pos)

     r1_record = (round(pos-ext)) > 4
     r2_record = (round(pos+ext)) < (nrows-5)
 

   ; estimate the spatial profile map (P) by modeling the spatial profile
   ; as a gaussian. let's define the width (sigma) of the gaussian
   ; according to the fwhm of the spatial profile as measure by the
   ; find_object routine. and take the center of the gaussian (x0) to be
                                ; at the object position pos (as
                                ; determined using find_object +
                                ; peakinfo).

   ; do these lines for any case, now for horne, before determine weights
     sigma = fwhm / 2.35482
     x0 = pos

     x = findgen(nrows)
     xdata = exp( -(x - x0)^2 / (2.*sigma^2) )*(x ge r1)*(x le r2) 
     P = (xdata) ## (fltarr(n)+1.)
     P = abs(P)
          
     Ptot = total(P, 2)
     P = P / (Ptot # (fltarr(nrows)+1.) )
      
     weight = P*(P gt 1E-3)   

 endif else begin
   
   ; default boxsprof
   ; figure out size to extract for object from fwhm
     width = fwhm * nsigma

   ; determine width for boxcar extraction
   ; half the width can not exceed the position of the object
     ext = (ceil(width/2.)) < pos < (nrows-1-pos)
   ; must be greater than four pixels
     r1 = (round(pos-ext)) > 4
     r2 = (round(pos+ext)) < (nrows-5)

     r1_record=r1
     r2_record=r2

  
 
     if r2 lt r1 then begin
         print, '(extract1d.pro) ERROR: invalid extraction extrema defined!'
         print, '(r1, r2, pos):', r1, r2, pos
         return, zero1d
      endif

     weight = slit.flux*0 

   ; for every y index that is between r1 and r2, set its value to 1 
     weight[*,r1:r2] = 1

endelse



; do these lines for any case
  sigma = fwhm / 2.35482
  x0 = pos

  x = findgen(nrows)
  xdata = exp( -(x - x0)^2 / (2.*sigma^2) )*(x ge r1)*(x le r2)
  P = (xdata) ## (fltarr(n)+1.)
  P = abs(P)

  M = (slit.mask EQ 0b)
; now define some arrays of the proper length which we will fill with
; the extracted 1-d spectrum, ivar, and bad pixel masks. 
  spec_num = fltarr(n)
  spec_denom = fltarr(n)
  var_num = fltarr(n)
  sky_num=var_num
  var_denom = fltarr(n)
  crmask = intarr(n)
  bitmask = intarr(n) + 255
  ormask = intarr(n)
  infomask = intarr(n)
  nbadpix = intarr(n)
  var2d = slit.ivar*0. 
  ivar = fltarr(n)
  spec = ivar


  shiftarr = fltarr(nrows)


if shift eq 2 then begin
    for i=0,nrows-1 do begin
    idx = interpol(findgen(n), lambda2d[*,i], cwave)
          idx = round(idx)
          slit.flux[*,i] = slit.flux(idx,i)
          slit.ivar[*,i] = slit.ivar(idx,i)
          weight[*,i] = weight(idx,i)

          whbad = where((idx lt 0) OR (idx gt n-1), ct)
          if ct gt 0 then begin
              slit.ivar[whbad,i] = 0
              weight[whbad,i]=0
          endif

      endfor 
endif 


; where shift equals 1 or 0
  if((shift eq 1) OR (shift eq 0)) then begin

      for i=0,nrows-1 do begin
      ; JAN do this from 0 to nrows-1 instead
      ; determine the shift needed to rectify the rows in wavelength.
  
        if shift eq 0 then cshift = round((i-pos)*dxdp) $
            else cshift=round(mean(lambda2d[*,i]-lambda2d[*,pos])/dldx) 

      ; at a particular value of x, as you move along in y
      ; shift value changes

        shiftarr[i] = cshift

      ; find all of the points in the ith row which will be improperly
      ; wrapped when vector shift operations are applied.
        badpts = where(cwave[0] gt lambda2d[*,i] or $
                       cwave[n-1] lt lambda2d[*,i], bcnt)
      ; at all points where we improperly wrap, set the slit.ivar to zero.
        if bcnt gt 0 then slit.ivar[badpts,i] = 0.0
   
      ; fill the slit.var array with ivar at the good points....exclude the
      ; points where slit.ivar=0 since this will blow up upon inversion. the
      ; variance values at the bad pixels will be set at zero, but this
      ; isn't a concern since these variance values are not used anywhere.
        good = where(slit.ivar gt 0, goodcnt)
         if goodcnt gt 0 then slit.ivar[good] = 1./slit.ivar[good]


         slit.flux[*,i] = shift(slit.flux[*,i],cshift)

         weight[*,i] = shift(weight[*,i], cshift)

         slit.ivar[*,i] = shift(slit.ivar[*,i], cshift)

       ; sky_num = sky_num + shift(sky_num[i], cshift)

     endfor
 endif

; eval is a function that combines the different lambda information
  lambda2d=lambda_eval(slit)

  lambda = lambda2d[*,pos]
; convert lambda to angstroms
  lambda0 = lambda / 1E4


  slitpa=sxpar(hdr,'SLITPA')

  parang = sxpar(hdr,'PARANG')
 
  geom = cos(!dtor*(parang - slitpa))

  chipno = sxpar(hdr,'CHIPNO')

  el = sxpar(hdr,'EL')

; 90 degrees minus elevation
  z = 90. - el
; z to radians
  z = !dtor*z  
   
; (ie 610 millibar); press in mm(Hg)
  p = 465.
; temp in Celcius 
  t = 0.
  q = (1./lambda0) ^ 2 

; if chipno lt 5 then q0 = q[n-1] else q0 = q[0] 

  if chipno lt 5 then q0 = q[2730] else q0 = q[2048] 

  delndx_std = (64.328 + 29498.1/(146.-q) + 255.4/(41.-q)) * (1E-6)

  delndx_std0 = (64.328 + 29498.1/(146.-q0) + 255.4/(41.-q0)) * (1E-6)

  adjust =  p * (1. + p * (1.049*(1E-6)) - t * 0.0157*(1E-6)) $ 
     / (720.883 * (1.+ t * 0.003661))

  delndx = adjust * delndx_std

  delndx0 = adjust * delndx_std0
 
  r = delndx / (1. + 2.*delndx)
  r0 = delndx0 / (1. + 2.*delndx0)
 
; refraction = r * tan(z)

  delta_pix = geom * tan(z) * (r - r0) * 206265. / 0.1192


  for i = 0, n -1 do begin
    ; weight[i,*] = shift(weight[i,*],round(delta_pix[i]))
      weight[i,*] = interpolate(weight[i,*],$
         findgen(nrows)-delta_pix[i],missing=0)
  endfor




  if shift eq 2 then begin
      for i=0,nrows-1 do begin
          idx = interpol(findgen(n), lambda2d[*,i], cwave)
          idx = round(idx)
          crmask = crmask + (weight[*,i] gt 0)*(slit.crmask[idx,i] gt 0b)
         
          bitmask  =  ((weight[*,i] eq 0)* bitmask) $
                          or ((weight[*,i] gt 0) * ( bitmask and slit.mask[idx,i]))
         
          ormask = ormask or  (weight[*,i] gt 0)* slit.mask[idx,i]

          nbadpix = nbadpix + (weight[*,i] gt 0)*(slit.mask[idx,i] ne 0b)

          infomask = infomask or(weight[*,i] gt 0) * slit.infomask[idx,i]
       endfor
  endif


  if((shift eq 1) OR (shift eq 0)) then begin
       for i=0,nrows-1 do begin
        ; determine the number of pixels at each wavelength which have a CR
        ; hit in them.
           cshift=shiftarr[i]

          crmask = crmask + weight[*,i]*shift(slit.crmask[*,i] gt 0, cshift)
        ; take the minimum bitmask value (along the row) at each wavelength.
          bitmask =  ((weight[*,i] eq 0)* bitmask) $
             or ((weight[*,i] gt 0)* (bitmask and shift(slit.mask[*,i], cshift)))
        ; also take the OR of the bitmask.
          ormask = ormask or (weight[*,i] gt 0) * shift(slit.mask[*,i], cshift)
        ; determine the number of badpixels at each wavelength.
          nbadpix = nbadpix +(weight[*,i] gt 0) * shift(slit.mask[*,i] ne 0b, cshift)
        ; take the OR of the infomask.
          infomask = infomask or (weight[*,i] gt 0) * shift(slit.infomask[*,i], cshift)
       endfor  
  endif





; ADAM: do this after the new for loops
  weight = weight *(slit.ivar ne 0)  

  if ivarwt eq 0 then weight_ivar = 1 $
      else weight_ivar=slit.ivar

  if ivarwt eq 0 then variance= 1/slit.ivar $
      else variance=1

  if n_elements(variance) gt 1 then variance[where(slit.ivar eq 0)]=0
 
  totweight = total(weight, 2)


; do these lines for boxsprof, after have shifted & determined weight 
; locate bad pixels in the flux array and interpolate over them.
; interpolate across CTE problems, bad spots in pixmap, in spatial
; direction. do not interpolate across vignetted regions.
  isbadpix = (slit.mask and 22b) gt 0
   
; if the boxsprof keyword is set, then do NOT interpolate. just set
; the bad pixels to zero and we will compensate for them using the
; spatial profile. for the horne and optimal extractions, set the ivar
; for interpolated (isbadpix)  pixels to zero.
   
  P = slit.flux * (1 - isbadpix) 
; P = djs_maskinterp(slit.flux, interpolate, iaxis=1, /const)
; interpolate the ivar values as well.
; ivarout = slit.ivar * (1 - isbadpix) 
; else ivar2d = djs_maskinterp(slit.ivar, interpolate, iaxis=1, /const)
  interpix = where(isbadpix, inter_cnt)
; increase the variance in the interpolated pixels, but exact factor
; depends on number of rows interpolated over.
; if boxcar then begin
;   if inter_cnt gt 0 then ivar2d[interpix] = ivar2d[interpix]/4.
; endif else begin
; if inter_cnt gt 0 then ivar2d[interpix] = 0.
; endelse

; figure out fraction of light lost to bad pixels
  if r2-r1 ge 1 then begin
       intprof = total(P * isbadpix*weight, 2) 
  endif else intprof = P * isbadpix
    

; boxsprof, figure out correction factor
; remember to do this and figure out isbadpix after have shifted weight, flux, etc.
; intscl is use to find the correction
  intscl = intprof / total(P*weight,2)
  intscl = 1.0 / (1.0 - intscl)

; handle division by 0 in intscl
  whbad=where(finite(intscl) eq 0,badct)
  if badct gt 0 then intscl[whbad]=1.


  if horne gt 0 then correction=1. else correction=median(totweight,/even)*intscl

; denominator=(weight)/slit.ivar
; denominator=(weight^2)/slit.ivar

  denominator = weight^2*weight_ivar

  denominator[where(slit.ivar eq 0)]=0

; denominator=total(denominator * weight_ivar,2)
  denominator=total(denominator,2)

  numerator=total(weight^2*variance*weight_ivar,2)

; calculating the output inverse variance, now weighted by ivar_weight
; ivarout = total(weight^2 * weight_ivar, 2 ) /denominator

  sky_num = total(weight * slit.ivar, 2)

  var_denom = total(weight^2 * slit.ivar, 2)

  count = total(weight *(1- isbadpix), 2)

  if horne eq 1 then f = (sky_num / var_denom) else f = count


; correction, dividing by the weight reduces it by a factor
; ivarout = ivarout* (1/correction^2)

  ivarout=1/(numerator/denominator^2*correction^2 + f^2/skyivar)

; updated to match horne code in extract1d
  totweight=total(weight^2*weight_ivar,2)

; summing up all of the y values of the weights and normalizing
  spec = total(slit.flux*weight*weight_ivar,2)  / totweight*correction

  whinfinite1=where(finite(spec) eq 0,ct)

  whinfinite2=where(finite(ivarout) eq 0 OR finite(spec) eq 0,ct2)

  if ct gt 0 then spec[whinfinite1]=0

  if ct2 gt 0 then ivarout[whinfinite2]=0       

  flux=slit.flux

RETURN,  {SPEC:spec,LAMBDA:lambda, IVAR:float(ivarout), CRMASK:crmask, BITMASK:bitmask, $
         ORMASK:ormask, NBADPIX:nbadpix,INFOMASK:infomask, OBJPOS:float(pos), $
          FWHM:float(fwhm), NSIGMA:nsigma, R1:r1_record, R2:r2_record, SKYSPEC:skyspec,$
           IVARFUDGE:ivarfudge(spec,ivarout)}
END         


