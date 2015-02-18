;+
; NAME:
;  coadd_spectra
;
; PURPOSE:
;  sums 1d spectra with ivar weighting to get composite having
;  improved SNR
;
; CALLING SEQUENCE:
;  result = coadd_spectra(list,normalize=normalize)
; 
; INPUTS:
;  list -- list of zresult structures tabulating files to coadd
;  
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;  normalize -- if set, spectrum is normalized to unity prior to
;               coaddition
;  zeroweight -- if set, Aband and Bband regions are given 0 ivar
;  ivartrim   -- if set, throw out points with worst ivar, as used to
;                do (can be problematic)
;  inorm   -- normalize all ivar to have mean 1 (~equally weight
;                each spectrum)
;  dloglam -- change sample spacing (default: 2E-5)
; OUTPUTS:
;   result -- a structure of the summed spectra with tags for spectra,
;             lambda, ivar
; 
; COMMENTS:  
;   this routine is useful for improving SNR by coadding many files.
;   A-band and B-band regions given 0 weight
;
; REVISION HISTORY:
;  md, rs 12dec02
;  md  29jan03 ivar=0. in Aband, Bband regions
;----------------------------------------------------------------------

function coadd_spectra, list, normalize=normalize, zeroweight=zeroweight, $
                        ivartrim=ivartrim,inorm=inorm,dloglam=dloglam

   nspec = n_elements(list)
   rlamlimit = [2600., 9300.]
   llimit = alog10(rlamlimit)

   if n_elements(dloglam) eq 0 then dloglam = 2.e-5
   npoints = long((llimit[1]-llimit[0])/dloglam)
   specsum = fltarr(npoints)
   norm=fltarr(npoints)
   llambda = findgen(npoints)*dloglam +llimit[0]
   lambda = 10.^llambda
   ivar = lambda*0.
   if n_elements(ivartrim) eq 0 then ivartrim = 0
   if n_elements(inorm) eq 0 then inorm = 0

   rescale=fltarr(nspec)

   for i=0, nspec-1 do begin ;begin loop over input list

     datestring=filenamefromdate(list[i].date)
     path=concat_dir(getenv('D2_RESULTS'),strcompress(list[i].maskname+'/'+datestring+'/',/remove))
     spec=(findfile(concat_dir(path,strcompress('spec1d*'+list[i].objname+'*',/rem))))[0]

     print, spec
     z_in = list[i].z

     ss1d=fill_gap(spec,/tweak,/telluric,/silent)

; correct for DEIMOS response

     fix_response,ss1d
     wave=ss1d.lambda
     airtovac,wave
     ss1d.lambda=wave

     spec_in=ss1d.spec
     lambda_in=ss1d.lambda
     ivar_in=ss1d.ivar

     badpix = ss1d.nbadpix gt 2 OR (ss1d.infomask AND 21b) gt 0 $
       OR (ss1d.ormask AND 30b) gt 0 

     badpix=dilate(badpix,intarr(5)+1)

; sometimes want to give spectra ~equal weight, but deweight bad
; pixels - set IVAR to have median 1 in that case
     if inorm then begin
         rescale[i]=1./median(ivar_in[where(ivar_in NE 0)])
         ivar_in=ivar_in*rescale[i]
     endif

; the data has noise in excess of Poisson due to skysub problems. 
; Get rid of the worst 5% of the pixels, plus their immediate
; neighbors

     if ivartrim then begin

         awful = where(ivar_in lt 1.e-8)
         ivar_in[awful] = 0.
         bad = (ivar_in*0.) 
         good = bad             ; how is that?
         nlam = n_elements(lambda_in)
         cut = sort(ivar_in)
         junk = min(ivar_in[cut] gt 0, cs) ;find location of worst non-zero point
         cut = cut[cs:cs+nlam/20]
         bad[cut] = 1.
         bad = smooth(bad, 3)
         good = bad eq 0.
         ivar_in = ivar_in*good
     endif

     ivar_in=ivar_in*(badpix eq 0)

; account for interpolation issues
     ivar_in = ivar_in < smooth(ivar_in,3)


;before coadding, set regions of telluric absorption to have ivar=0
     if keyword_set(zeroweight) then begin
        aband = where(lambda_in gt 7592. and lambda_in lt 7675., naband)
        bband = where(lambda_in gt 6866. and lambda_in lt 6881., nbband)
        if (naband gt 0) then ivar_in[aband] = 0.
        if (nbband gt 0) then ivar_in[bband] = 0.
     endif

; set to unit average flux if desired
     if keyword_set(normalize) then begin
       sumspec = total(spec_in*ivar_in > 0.) ;only sum positive zones
       sumivar = total(ivar_in*(spec_in gt 0.));only weight positive
       norm = sumspec/sumivar
       spec_in = spec_in/norm   ;normalize
       ivar_in = ivar_in*norm^2 ;keep snr constant

       ss1d.spec = spec_in
     endif

     ss1d.ivar=ivar_in

     lambdarange=minmax(lambda_in/(1.+z_in))
     defined=lambda gt (lambdarange[0]+1) AND lambda lt (lambdarange[1]-1)
     objflux=defined*interpol(spec_in,lambda_in/(1.+z_in),lambda)
     objivar=defined*interpol(ivar_in,lambda_in/(1.+z_in),lambda)
     if inorm then objnorm = defined*rescale[i]

     ivar =  ivar + $
           objivar  ;accumulate ivar
     specsum =  specsum + $
           objflux*objivar

     if inorm then norm=norm + objnorm
     
   endfor
;
; now normalize 
     good = where(ivar gt 0)
     specsum[good] = specsum[good]/ivar[good]
     if inorm then ivar[good]=ivar[good]/norm[good]*nspec

   result = { spec: specsum, lambda: lambda, ivar: ivar}

  return, result
end
