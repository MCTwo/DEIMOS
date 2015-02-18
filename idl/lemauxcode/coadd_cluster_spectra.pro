
; NAME:
;  coadd_spectra
;
; PURPOSE:
;  sums 1d spectra with ivar weighting to get composite having
;  improved SNR
;
; CALLING SEQUENCE:
;  result = coadd_cluster_spectra(cat,normalize=normalize,zero=zero,inorm=inorm,ivartrim=ivartrim,weight=weight)
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
;                (makes little difference)
;  ivartrim   -- if set, throw out points with worst ivar, as used to
;                do (can be problematic - DEPRECATED)
;  inorm   -- normalize all ivar to have mean 1, and include only
;             variations in ivar due to sky intensity, not
;             obj. intensity (~equally weight
;                each spectrum)
;  dloglam -- change sample spacing (default: 2E-5)
;  weight  -- a weight to apply to each element of list in building
;             the coadded spectra.  If weights are used, NORMALIZE and
;             INORM are also set.
;
; OUTPUTS:
;   result -- a structure of the summed spectra with tags for spectra,
;             lambda, ivar
; 
; COMMENTS:  
;   this routine is useful for improving SNR by coadding many files.
;
; REVISION HISTORY:
;  md, rs 12dec02
;  md  29jan03 ivar=0. in Aband, Bband regions
; jan, ry 2004, 2005 major changes
; feb 2007 bl, cat input, minor revisions
;----------------------------------------------------------------------

function Ez, x

   distance = 1/sqrt(0.27*(1+x)^3+0.73)
   return, distance 

end

function coadd_cluster_spectra, cat, normalize=normalize, zeroweight=zeroweight, $
                        ivartrim=ivartrim,inorm=inorm,dloglam=dloglam,weight=weight

   numread = 207

   readcol, cat, format='A,A,F,F,A,A', mask, slit, Iband, z, q, file, skipline=0,  numline = numread
   
   list = replicate({list}, numread)
   for i=0, numread-1 do begin
   
   list[i] = {list, maskname:mask[i], slitname:slit[i], Iband:Iband[i], z:z[i], zquality:q[i], file:file[i], absIflux:0.0}
   endfor

   ;absolute i' magnitude calculation using omega_M = 0.27 omega_L = 0.73 and h = 0.71
   absIflux = fltarr(numread)
   for i=0, numread-1 do begin

   Dcnorm = qromb('Ez', 0.0, list[i].z)
   Dc = 4.225e9*Dcnorm
   lumdist = (1+list[i].z)*Dc
   M = list[i].iband - 5*alog10(lumdist/10)
   absIflux[i] = 10^(.4*(-20-M))
   list[i].absIflux = absIflux[i]
   endfor
   
   nspec = n_elements(list)
   rlamlimit = [3200., 4900.]
   llimit = alog10(rlamlimit)

   if n_elements(dloglam) eq 0 then dloglam = 5e-5
   npoints = long((llimit[1]-llimit[0])/dloglam)
   specsum = dblarr(npoints)
   norm=fltarr(npoints)
   llambda = findgen(npoints)*dloglam +llimit[0]
   lambda = 10.^llambda
   ivar = lambda*0.d0
   if n_elements(ivartrim) eq 0 then ivartrim = 0
   if n_elements(inorm) eq 0 then inorm = 0
   if n_elements(normalize) eq 0 then normalize=0
   useweight=n_elements(weight) gt 0
   if n_elements(weight) lt n_elements(list) then weight=fltarr(n_elements(list))+1.
   
   weight = 100*list.absIflux
   ;print, weight
   ;weight=weight/total(weight)
   totweight=specsum*1.d0
   totn=totweight
   inorm=inorm OR useweight
   normalize = normalize OR useweight

   rescale=fltarr(nspec)

   for i=0, nspec-1 do begin ;begin loop over input list

    spec = ('/Volumes/Data2/orelse/lemaux/deimos/sc1604/' + list[i].file)[0]
    ;datestring=filenamefromdate(list[i].date)
    ;path=concat_dir(getenv('D2_RESULTS'),strcompress(list[i].maskname+'/'+datestring+'/',/remove))
    ;spec=(findfile(concat_dir(path,strcompress('spec1d*'+list[i].maskname+'*'+list[i].objname+'*',/rem))))[0]

; some slosh in dates
    ;if strlen(spec) eq 0 then begin
    ;        datestring=strmid(datestring,0,strlen(datestring)-1)+'*'            
    ;        path=concat_dir(getenv('D2_RESULTS'),strcompress(list[i].maskname+'/'+datestring+'/',/remove))
    ;        spec=(findfile(concat_dir(path,strcompress('spec1d*'+list[i].maskname+'*'+list[i].objname+'*',/rem))))[0]
    ;endif

    z_in = list[i].z

    if strlen(spec) gt 1 then begin

     if inorm eq 0 then begin 
       ss1d=fill_gap(spec,/tweak,/telluric,/silent,header=header) 
       wave=ss1d.lambda
       airtovac,wave
       ss1d.lambda=wave
     endif else ss1d=fill_gap(spec,/tweak,/silent,header=header)

; corrected for DEIMOS response in fill_gap or inorm section

     spec_in=ss1d.spec
     lambda_in=ss1d.lambda
     ivar_in=ss1d.ivar

     badpix = ss1d.nbadpix gt 5 $
       OR (ss1d.ormask AND 30b) gt 0 

     ;badpix=dilate(badpix,intarr(5)+1)

     ivar_in=ivar_in*(badpix eq 0)     
       ss1d.ivar=ivar_in


; sometimes want to give spectra ~equal weight, but deweight bad
; pixels - set IVAR to have median 1 in that case, and have it scale
;          ONLY with sky spectrum
     if inorm OR useweight then begin
; do two sides separately



; fit for IVAR vs. sky counts
         wh=where(ivar_in NE 0 AND ss1d.skyspec gt 2*median(ss1d.skyspec) $
                  AND badpix eq 0 AND $
                  lindgen(n_elements(ss1d.skyspec)) $
                          lt n_elements(ss1d.skyspec)/2 ,ct)
         if ct gt 10 then begin
             fit=svdfit(ss1d.skyspec[wh],1./ivar_in[wh],2)
             wh1=where(lindgen(n_elements(ss1d.skyspec)) $
                          lt n_elements(ss1d.skyspec)/2 ,ct)
; adjust ivar
             if ct gt 0 then ivar_in[wh1]= $
               1./(fit[0]+fit[1]*(ss1d.skyspec[wh1]>0.))*(ivar_in[wh1] ne 0)
         endif    



; now do it on the other side
         wh=where(ivar_in NE 0 AND ss1d.skyspec gt 2*median(ss1d.skyspec) $
                  AND badpix eq 0 AND $
                  lindgen(n_elements(ss1d.skyspec)) $
                          ge n_elements(ss1d.skyspec)/2 ,ct)
         if ct gt 10 then begin
             fit=svdfit(ss1d.skyspec[wh],1./ivar_in[wh],2)
             wh2=where(lindgen(n_elements(ss1d.skyspec)) $
                          ge n_elements(ss1d.skyspec)/2 ,ct)
             if ct gt 0 then ivar_in[wh2]= $
               1./(fit[0]+fit[1]*(ss1d.skyspec[wh2]>0.))*(ivar_in[wh2] ne 0)
         endif    



         ivar_in=ivar_in*(ss1d.skyspec ne 0) 

         if total(ivar_in lt 0) gt 0 then ivar_in=ivar_in*0.
         whok=where(ivar_in ne 0,ct)
         if ct gt 0 then rescale[i]=1./median(ivar_in[whok]) else rescale[i]=1.
         ivar_in=ivar_in*rescale[i]

; need to do telluric corrections after adjusting ivar for best results
         airmass = sxpar(header, 'AIRMASS')
         ss1d.ivar=ivar_in
         ss1dtmp=ss1d


         remove_telluric, ss1d, airmass,silent=silent
         fix_response,ss1d

         ivar_in=ss1d.ivar
         wave=ss1d.lambda

         airtovac,wave
         ss1d.lambda=wave
         lambda_in=ss1d.lambda
         spec_in=ss1d.spec
     endif



; set to unit average flux if desired
     if normalize OR useweight then begin
         whgood=where(ivar_in gt 0,ct)
         if ct gt 1 then norm=median(spec_in[whgood]) else norm=-1.E10
         if norm lt 0 then ivar_in=ivar_in*0.
       spec_in = spec_in/norm   ;normalize
       ivar_in = ivar_in*norm^2 ;keep snr constant

       ss1d.spec = spec_in
       ss1d.ivar=ivar_in
     endif

; the data has noise in excess of Poisson due to skysub problems. 
; Get rid of the worst 5% of the pixels, plus their immediate
; neighbors
; THIS IS OLD CODE AND SHOULD BE A LAST RESORT!

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



; account for interpolation issues
     ivar_in = ivar_in < smooth(double(ivar_in),3)



;before coadding, set regions of telluric absorption to have ivar=0
     if keyword_set(zeroweight) then begin
        aband = where(lambda_in gt 7592. and lambda_in lt 7675., naband)
        bband = where(lambda_in gt 6866. and lambda_in lt 6881., nbband)
        if (naband gt 0) then ivar_in[aband] = 0.
        if (nbband gt 0) then ivar_in[bband] = 0.
    endif





     lambdarange=minmax(lambda_in/(1.+z_in))
     defined=lambda gt (lambdarange[0]+1) AND lambda lt (lambdarange[1]-1)
     lambda_in_rest = lambda_in/(1.+z_in)
     objflux=defined*interpol(spec_in,lambda_in_rest,lambda);,/spline)
     objivar=defined*interpol(ivar_in,lambda_in_rest,lambda);,/spline)

     if inorm then objnorm = defined*rescale[i]
     totn=totn+(objivar gt 0)

     whdefined=where(defined gt 0)

     if useweight then ivar =  ivar + $
           objivar*weight[i]  $
     else ivar =  ivar + $
           objivar  ;accumulate ivar

     if useweight then specsum =  specsum + $
           objflux*objivar*weight[i] $
       else specsum =  specsum + $
           objflux*objivar

     if useweight then totweight=totweight+weight[i]*(objivar ne 0)
     
     if inorm then norm=norm + objnorm
     
    endif else print,'missing '+string(list[i].maskname)+' '+string(list[i].slitname)

   endfor
;
; now normalize 
     good = where(ivar gt 0)
     
     specsum[good] = specsum[good]/ivar[good]

     if inorm then ivar[good]=ivar[good]/norm[good]*totn[good]

     if useweight then ivar[good]=ivar[good]/totweight[good]

   result = { spec: specsum, lambda: lambda, ivar: ivar, totn:totn}
   mwrfits, result, 'fullcoadd_208galaxies.fits'
   device, retain=2
   splot, result.lambda, result.spec
  return, result
end
