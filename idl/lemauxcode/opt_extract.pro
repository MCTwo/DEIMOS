Function opt_extract,counts,lam,invvar,weight,speclam=speclam,$
               specvar=specvar,modelcounts=modelcounts,lam0=lam0,$
               bin=bin,loglin=loglin,clip=clip,noreject=noreject,$
               maxcount=maxcount,smooth=smooth,medfilt=medfilt,$
               boxcar=boxcar,disp=disp

;+--------------------------------------------------------------------
;
; OPT_EXTRACT   4/2003
;
; Perform optimal extraction using Horne 1986 algorithm.
;
; Note on weights: The weight given to a pixel is the fraction of the
; total flux from lam1 to lam2 expected to fall on that pixel, where
; lam1 and lam2 are the minimum and maximum wavelengths, respectively,
; spanned by the pixel.  In the case where the spectrum is being
; extracted in bins that are larger than the wavelength range of a
; single pixel (bin size > the dispersion at that pixel), the program
; automatically adjusts the weights to reflect that the fact that
; the flux in a bin is divided among more pixels.
;
; Note on binning: If the disp keyword is set (and it almost always
; should be), then the program is currently geared to work only if the
; bin size is everywhere at least 0.5 times the dispersion.  One
; normally wants to extract with bin sizes > or = the dispersion, so
; this should not generally a problem.  Limits on binning are so that
; the routine can be optimized while properly accounting for input
; pixels that fall partially in ouput bins.
;
; Warning: Boxcar and optimal extraction should give the same results
; for bright objects.  However, it can be difficult to compute the
; profile for a bright object accurately enough to avoid systematic errors,
; which are usually over- or underestimating the contributions from the
; profile wings.  A good check is to run both optimal and boxcar 
; extraction and compare.
;
; INPUTS
;       The following are 1-D or 2-D arrays of uniform size
;       giving information for each pixel:
;          counts   - Sky-subtracted counts
;          lam      - Wavelength (lambda)
;          invvar   - Inverse variance of counts.  Pixels whose
;                     inverse variance = 0 will be masked.
;          weight   - Weights determined from the object profile.
;                     See above.
;
; OUTPUTS
;          <result> - Optimally extracted spectrum, binned
;                     by wavelength
;
; KEYWORDS
;          speclam  - Wavelength array corresponding to binned
;                     spectrum.  Wavelengths are given for the
;                     center of each bin.
;          specvar  - Array of formal variance values for binned
;                     spectrum. 
;          modelcounts - Array with same size as input arrays,
;                     giving expected number of counts based
;                     the best spectrum estimate
;          bin      - Bin size, in wavelength units unless the loglin
;                     keyword is set (see below).  Default is 1.
;          loglin   - If set, binning is log-linear, with 'bin' in
;                     km/s.
;          lam0     - If set, bins will be placed such that lam0
;                     occurs at the center of a bin.  Default
;                     is 1000.
;          clip     - Threshold fpr rejecting outlying pixels from
;                     estimate, in numbers of standard deviations.
;                     Default is 5.
;          noreject - Do not perform iterative rejection of outlying
;                     pixels in each bin.
;          maxcount - Scalar, pixels having an absolute value greater
;                     than this are rejected.  Default is no threshold.
;          smooth   - Scalar, FWHM in wavelength units of Gaussian
;                     kernel by which to smooth extracted spectrum
;                     before generating modelcounts.  Useful, e.g.,
;                     for smoothing out poorly constrained absorption
;                     feautres or artifacts from skyline subtraction.
;                     Note: will not affect output extracted spectrum.
;          medfilt  - Scalar, size of median filter by which to
;                     smooth extracted spectrum before generating model 
;                     counts.  Will not affect output spectrum.
;                     Should not set both this keyword and the smooth
;                     keyword.
;          boxcar   - Use boxcar extraction instead of optimal
;                     extraction (weighting is ignored).  This is
;                     the default if only three input parameters are 
;                     provided.  Note: Boxcar extraction may correlate
;                     adjacent spectral bins, as partial input pixels
;                     are used and a single input pixel may contribute
;                     to more than one output bin.
;          disp     - RECOMMENDED. Dispersion, in wavelength units per
;                     pixel (unless the 'loglin' keyword is set - see
;                     below).  May either be a scalar (constant
;                     dispersion) or an array of the same size as
;                     counts giving the dispersion at each pixel).  If
;                     supplied, pixels that fall on the boundary
;                     between binned pixels will be sub-sampled.
;                     Default is no sub-sampling (assume infinitely
;                     narrow input pixels).  Should be smaller than
;                     the output bin size.
;
; HISTORY
;          Written 4/15/2003 GDB
;          BOXCAR and DISP keywords added 9/6/2004 GDB
;          LOGLIN keyword added 5/4/2007 GDB
;          MEDFILT keyword added 8/30/2007 GDB
;----------------------------------------------------------------------

if (n_params() lt 3) then begin
   print,'CALLING SEQUECE:'
   print,'   spectrum = opt_extract(counts,lam,invvar,weight,speclam=speclam,'
   print,'               specvar=specvar,modelcounts=modelcounts,lam0=lam0,'
   print,'               bin=bin,/loglin,clip=clip,/noreject,'
   print,'               maxcount=maxcount,smooth=smooth,medfilt=medfilt,'
   print,'               /boxcar,disp=disp)'
   return,-1
endif

; Speed of light (km/s)
c = 2.99792458d5

; If fewer than four inputs provided, assume boxcar extraction
if (n_params() eq 3) then begin
   print,'No weighting provided - doing boxcar extraction.'
   boxcar = 1
   weight = 1 + 0.*counts
endif

; Check that input arrays all have the same number of elements
n_counts = n_elements(counts)
n_lam    = n_elements(lam)
n_invvar = n_elements(invvar)
n_weight = n_elements(weight)
if (n_lam ne n_counts or  n_invvar ne n_counts or n_weight ne n_counts) $
   then begin
   print,'OPT_EXTRACT: Error - Input arrays must all be the same number'
   print,'             of elements.'
   return,-1
endif
  
; Check that variance array is everywhere non-negative
ls = where(invvar lt 0)
if (total(ls) ne -1) then begin
   print,'OPT_EXTRACT: Error - Inverse variance array should contain only'
   print,'                     non-negative values.'
   return,-1
endif

; Dispersion (optional)
if keyword_set(disp) then begin
   if (n_elements(disp) eq 1) then begin
      disper = 0.*counts
      disper(*) = disp
   endif else disper = disp
endif else disper = 0.*counts
; Calculate the wavelength at the edges of each pixel
lam_edge1 = lam - abs(disper) / 2.
lam_edge2 = lam + abs(disper) / 2.

; Reject pixels that have zero weight (optimal extraction only),or
; whose inverse variance is zero, or (optionally) whose counts exceed
; the specified limit
if keyword_set(boxcar) then begin
   if keyword_set(maxcount) then workpix = where(invvar ne 0 and $
                                                 abs(counts) lt maxcount) $
   else workpix = where(invvar ne 0)
endif else begin
   if keyword_set(maxcount) then workpix = where(weight ne 0 and $
                                                 invvar ne 0 and $
                                                 abs(counts) lt maxcount) $
   else workpix = where(weight ne 0 and invvar ne 0)
endelse
workcounts    = counts(workpix)
worklam       = lam(workpix)
worklam_edge1 = lam_edge1(workpix)
worklam_edge2 = lam_edge2(workpix)
workinvvar    = invvar(workpix)
workweight    = weight(workpix)
workdisper    = disper(workpix)

; Set bin size.  Avoid floating point errors - assume that bin size is
; > 1e-6.
if keyword_set(bin) then binsz=double(bin) else binsz=1d0 
binsz = round(1d6*binsz)/(1d6)

; Reference bin center
if keyword_set(lam0) then lamcent=double(lam0) else lamcent=1000d0

; Construct wavelength array.  Index j would be the array index
; for wavelength speclam(j) if speclam(0) = lamcent.  Also calculate
; wavelengths at the edge of each bin.
lam_min = min(worklam_edge1)
lam_max = max(worklam_edge2)
if keyword_set(loglin) then begin
   logbinsz   = alog10(binsz/c + 1)
   loglamcent = alog10(lamcent)
   j_min      = floor((alog10(lam_min/lamcent)/logbinsz) + 0.5)
   j_max      =  ceil((alog10(lam_max/lamcent)/logbinsz) - 0.5)
   nbins      = j_max - j_min + 1
   j_list     = j_min + dindgen(nbins)
   speclam    = 10^(loglamcent + j_list*logbinsz)
   speclam_lo = 10^(loglamcent + (j_list-0.5)*logbinsz)
   speclam_hi = 10^(loglamcent + (j_list+0.5)*logbinsz)
endif else begin
   j_min      = floor((lam_min-lamcent)/binsz + 0.5)
   j_max      =  ceil((lam_max-lamcent)/binsz - 0.5)
   nbins      = j_max - j_min + 1
   j_list     = j_min + dindgen(nbins)
   speclam    = lamcent + j_list*binsz
   speclam_lo = lamcent + (j_list-0.5)*binsz
   speclam_hi = lamcent + (j_list+0.5)*binsz
endelse

; Wavelength range per spectrum pixel
wavperpix = speclam_hi - speclam_lo

; Check that the bin size is everywhere > 0.5 * the dispersion
if keyword_set(disp) then begin
   nearest_binsz = interpol(wavperpix,speclam,worklam)
   toosmall_ls   = where(nearest_binsz lt 0.5*workdisper,n_toosmall)
   if (n_toosmall ne 0) then begin
      print,'OPT_EXTRACT: Bin size must everywhere be > 0.5 * dispersion'
      return,-1
   endif
endif

; Adjust the weighting to reflect the fact that the bin
; size may not be equal to the dispersion, and therefore multiple
; pixels may fall in the same bin at the same position along
; the object profile (may occupy the same region in the
; probability distribution).
workweight = workweight * workdisper / interpol(wavperpix,speclam,worklam)

; Step through wavelengths and perform weighted (or boxcar) sums.

spec    = dblarr(nbins)
specvar = dblarr(nbins)
specpix = fltarr(nbins)  ; Count of the total # of pixels in each bin

modelcounts     = 0*counts
workmodelcounts = 0*workcounts

; Convert input wavelengths to the index values of the output spectrum
; bins in which they would fall.
if keyword_set(loglin) then begin
   w0 = alog10(speclam(0))
   dw = logbinsz
   lamindex = (alog10(worklam) - w0) / dw
endif else begin
   w0 = speclam(0)
   dw = binsz
   lamindex = (worklam - w0) / dw
endelse

lam_hist = histogram(lamindex,min=-0.5,nbins=nbins,rev=lam_ri)

for i=0L,long(nbins)-1 do begin

   if (lam_ri(i+1) eq lam_ri(i)) then continue

   ; All pixels that contribute to the current bin.  Pixel centroids
   ; are guaranteed to be in either the current bin or one of the
   ; adjacent bins so long as the bin size is everywhere > 0.5 *
   ; dispersion (checked above).
   ilo     = (i-1) > 0
   ihi     = (i+1) < (nbins-1)
   temp_ls = lam_ri(lam_ri(ilo):lam_ri(ihi+1)-1)

   ; Compute the fraction of each pixel that falls in the bin.
   ; Note: Pixels that fall entirely outside the bin will have
   ; negative fractions.
   temp_lam_edge1 = worklam_edge1(temp_ls)
   temp_lam_edge2 = worklam_edge2(temp_ls)
   temp_disper    = workdisper(temp_ls)
   lamlo     = speclam_lo(i)
   lamhi     = speclam_hi(i)
   temp_frac = ((temp_lam_edge2 < lamhi) - (temp_lam_edge1 > lamlo)) / $
                 temp_disper

   ; If running optimal extraction, include only pixels that are at
   ; least half covered by the current bin.  (The requirement above
   ; that the binning be everywhere at least 0.5 times the dispersion
   ; guarantees that no pixels are missed on account of being less
   ; than half covered by a single bin.)  If running boxcar
   ; extraction, then include all pixels that overlap the current bin,
   ; and scale the counts in each pixel by the fraction that is in the
   ; bin.  Note that this method of boxcar extraction may correlate
   ; adjacent pixels.
   if keyword_set(boxcar) then begin
      include_ls = where(temp_frac gt 0, n_include)
   endif else begin
      include_ls = where(temp_frac gt 0.5, n_include)
   endelse
   if (n_include eq 0) then continue
   bin_ls        = temp_ls(include_ls)
   bin_frac      = temp_frac(include_ls)
   bin_counts    = workcounts(bin_ls)
   bin_weight    = workweight(bin_ls)
   bin_invvar    = workinvvar(bin_ls)
   bin_lam_edge1 = worklam_edge1(bin_ls)
   bin_lam_edge2 = worklam_edge2(bin_ls)
   bin_disper    = workdisper(bin_ls)

   ; Legacy naming conventions
   slicecounts = bin_counts
   sliceweight = bin_weight
   sliceinvvar = bin_invvar

   ; Total number of pixels in this bin
   specpix(i) = n_include

   ;
   ; Boxcar extraction (optional)
   ;
   
   if keyword_set(boxcar) then begin

      spec(i)    = total(slicecounts*bin_frac)
      slicevar   = sliceinvvar^(-1)
      specvar(i) = total(slicevar*(bin_frac^2))

   endif else begin

   ;
   ; Optimal extraction (default)
   ;

      ; Compute initial spec estimate   
      fluxsum   = total(sliceweight*slicecounts*sliceinvvar)
      weightsum = total((sliceweight^2)*sliceinvvar)
      estimate  = fluxsum/weightsum
      variance  = 1./weightsum

      ; Iteratively check for and reject outliers
      if keyword_set(noreject) then checkoutliers=0 $
         else checkoutliers = 1
      if keyword_set(clip) then nsig=clip else nsig=5
      
      while checkoutliers do begin
         excess = slicecounts - sliceweight*estimate
         badpix = where((excess^2)*sliceinvvar gt nsig^2,n_bad)
         if (n_bad ne 0) then begin
            ; Reject worst outlier and recompute spec estimate
            maxout  = max((excess^2)*sliceinvvar)
            goodpix = where((excess^2)*sliceinvvar ne maxout)
            slicecounts = slicecounts(goodpix)
            sliceweight = sliceweight(goodpix)
            sliceinvvar = sliceinvvar(goodpix)
            fluxsum     = total(sliceweight*slicecounts*sliceinvvar)
            weightsum = total((sliceweight^2)*sliceinvvar)
            estimate  = fluxsum/weightsum
            variance  = 1./weightsum
         endif else begin
            ; End loop if no more outliers
            checkoutliers = 0
         endelse
      endwhile
      
      spec(i)    = estimate
      specvar(i) = variance

   endelse

endfor

; Compute expected counts.  

; Optional: Smooth extracted spectrum before calculating counts.
if keyword_set(smooth) then begin
   kerfwhm = double(smooth) / binsz
   ker = gaussker(fwhm=kerfwhm)
   smspec = convol(spec,ker,/edge_truncate)
endif else smspec = spec

; Optional: Median filter extracted spectrum before calculating
; counts.
if keyword_set(medfilt) then smspec = median(smspec,medfilt)

; Fill in model array with expected counts
for i=0L,long(nbins)-1 do begin
   if (lam_ri(i+1) gt lam_ri(i)) then begin
      bin_ls = lam_ri(lam_ri(i):lam_ri(i+1)-1)
      workmodelcounts(bin_ls) = smspec(i)*workweight(bin_ls)
   endif
endfor
modelcounts(workpix) = workmodelcounts

; Divide spectrum by wavelength range per pixel to get flux / unit
; wavelength
spec      = spec / wavperpix
specvar   = specvar / (wavperpix^2)

; Return spectrum
return,spec

end
