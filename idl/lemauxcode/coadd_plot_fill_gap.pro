;+
; ROUTINE:  findex
;
; PURPOSE:  Compute "floating point index" into a table using binary
;           search.  The resulting output may be used with INTERPOLATE.
;
; USEAGE:   result = findex(u,v)
;
; INPUT:    
;   u       a monitically increasing or decreasing 1-D grid
;   v       a scalor, or array of values
;
; OUTPUT:
;   result  Floating point index. Integer part of RESULT(i) gives
;           the index into to U such that V(i) is between
;           U(RESULT(i)) and U(RESULT(i)+1).  The fractional part
;           is the weighting factor
;
;                          V(i)-U(RESULT(i))
;                          ---------------------
;                     U(RESULT(i)+1)-U(RESULT(i))
;
;
; DISCUSSION: 
;           This routine is used to expedite one dimensional
;           interpolation on irregular 1-d grids.  Using this routine
;           with INTERPOLATE is much faster then IDL's INTERPOL
;           procedure because it uses a binary instead of linear
;           search algorithm.  The speedup is even more dramatic when
;           the same independent variable (V) and grid (U) are used
;           for several dependent variable interpolations.
;
;  
; EXAMPLE:  
;
;; In this example I found the FINDEX + INTERPOLATE combination
;; to be about 60 times faster then INTERPOL.
;
;  u=randomu(iseed,20000) & u=u(sort(u))
;  v=randomu(iseed,10)     & v=v(sort(v))
;  y=randomu(iseed,20000) & y=y(sort(y))
;
;  t=systime(1) & y1=interpolate(y,findex(u,v)) & print,systime(1)-t
;  t=systime(1) & y2=interpol(y,u,v)            & print,systime(1)-t
;  print,f='(3(a,10f7.4/))','findex:   ',y1,'interpol: ',y2,'diff:     ',y1-y2
;
; AUTHOR:   Paul Ricchiazzi                        21 Feb 97
;           Institute for Computational Earth System Science
;           University of California, Santa Barbara
;           paul@icess.ucsb.edu
;
; REVISIONS:
;
;-
;
function findex,u,v
nu=n_elements(u)
nv=n_elements(v)

us=u-shift(u,+1)
us=us(1:*)
umx=max(us,min=umn)
if umx gt 0 and umn lt 0 then message,'u must be monotonic'
if umx gt 0 then inc=1 else inc=0

maxcomp=fix(alog(float(nu))/alog(2.)+.5) 

; maxcomp = maximum number of binary search iteratios

jlim=lonarr(2,nv)
jlim(0,*)=0          ; array of lower limits
jlim(1,*)=nu-1       ; array of upper limits

iter=0
repeat begin
  jj=(jlim(0,*)+jlim(1,*))/2
  ii=where(v ge u(jj),n) & if n gt 0 then jlim(1-inc,ii)=jj(ii)
  ii=where(v lt u(jj),n) & if n gt 0 then jlim(inc,ii)=jj(ii)
  jdif=max(jlim(1,*)-jlim(0,*))
  if iter gt maxcomp then begin
    print,maxcomp,iter, jdif
    message,'binary search failed'
  endif
  iter=iter+1
endrep until jdif eq 1 

w=v-v
w(*)=(v-u(jlim(0,*)))/(u(jlim(0,*)+1)-u(jlim(0,*)))+jlim(0,*)

return,w
end

pro smooth_ivar, arrayin, ivarin, sigma, arrayout, ivarout
     array = arrayin
     ivar = ivarin
     if n_elements(array) ne n_elements(ivar) then begin
         print, 'Wrong array dimensions in SMOOTH_IVAR.'
     endif
     asize = n_elements(array)
     n = floor(8. * sigma) + 1
     if n gt asize then begin
         message, 'Array too small to smooth.'
     endif
     half_n = (n-1)/2
     arraysm = dblarr(asize)
     ivarsm = dblarr(asize)
     ;window = dblarr(n) + 1d                                                                                   ;boxcar
     window = (1. / (double(sigma)) * sqrt(2.*!PI)) * exp(-(indgen(n)-double(half_n))^2./(2.*double(sigma)^2.))   ;Gaussian
     for i=0,half_n-1 do begin
         if total(ivar[0:i+half_n]*window[n-1-i-half_n:n-1]) ne 0.0 then begin
             arraysm[i] = total(array[0:i+half_n]*ivar[0:i+half_n]*window[n-1-i-half_n:n-1]) / total(ivar[0:i+half_n]*window[n-1-i-half_n:n-1])
             ivarsm[i] = (total(ivar[0:i+half_n]*window[n-1-i-half_n:n-1]))^2. / total(ivar[0:i+half_n] * (window[n-1-i-half_n:n-1])^2.)
         endif
     endfor
     for i=half_n,asize-half_n-1 do begin
         if total(ivar[i-half_n:i+half_n]*window) ne 0.0 then begin
             arraysm[i] = total(array[i-half_n:i+half_n]*ivar[i-half_n:i+half_n]*window) / total(ivar[i-half_n:i+half_n]*window)
             ivarsm[i] = (total(ivar[i-half_n:i+half_n]*window))^2. / total(ivar[i-half_n:i+half_n] * window^2.)
         endif
     endfor
     for i=asize-half_n,asize-1 do begin
         if total(ivar[i-half_n:asize-1]*window[0:asize-i-1+half_n]) ne 0.0 then begin
             arraysm[i] = total(array[i-half_n:asize-1]*ivar[i-half_n:asize-1]*window[0:asize-i-1+half_n]) / total(ivar[i-half_n:asize-1]*window[0:asize-i-1+half_n])
             ivarsm[i] = (total(ivar[i-half_n:asize-1]*window[0:asize-i-1+half_n]))^2. / total(ivar[i-half_n:asize-1] * (window[0:asize-i-1+half_n])^2.)
         endif
     endfor
     arrayout = arraysm
     ivarout = ivarsm
end

; This function blueshifts an array of (nspec) spectra by assuming
; that the serendip line is Lya.  It finds the maximum flux in a 10
; Angstrom window, and calls that rest wavelength 1215.68 Angstroms.
; Then, it interpolates the counts and inverse variances onto a
; standard 20,000 element grid with wavelength spacing d_lambdat.  The
; interpolation weights the flux in each bin by its inverse variance.
function Standardize_Spectra, $
        lambdas, $	; wavelength solution of each spectrum in Angstroms [nspec, *], restframe
        specs, $ ; flux, arbitrary units [nspec, *]
        ivars, $ 	; inverse variance, in ivar units of spectra, [nspec, *]
        z	; guess at the wavelength in Angstroms of the Lya line
    n = n_elements(z)                       ; number of individual spectra to blueshift and rebin

    d_lambdat = (5100.-3200.)/12000. ; wavelength spacing in Angstroms of coadded spectrum
    lambdat = d_lambdat * dindgen(12000) + 3200. ; wavelength array of coadded spectrum
    rspec = dblarr(n,12000) ; individual blueshifted and rebinned spectra [nspec, 12000]
    rivar = dblarr(n,12000) ; associated inverse variances [nspec, 12000]
    Lya_strength = dblarr(n)

    for i=0,n-1 do begin        ; for each individual spectrum
        
	; correct for the telluric absorption and the ccd gap  
        z_in = z[i]
	
        ss1d=fill_gap_nofitstable(specs,/telluric,/silent)

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

	;smoothspecs = smooth(specs[i,*], 3)
        sig = 1.17741 * (1.95)/(1+z[i])
	;sig = 2
	smooth_ivar, spec_in[i,*], ivar_in[i,*], sig, smoothspecs, smoothivars
        lambda = lambda_in[i,*] / (1.+z[i])
            
        spec = specs[i,*] ; individual spectrum yet to be blueshifted and rebinned
        ivar = ivars[i,*]       ; associated inverse variance
        ;print, ivar

        w = findex(lambda, lambdat) ; find the index mapping lambda -> lambdat
        ix = floor(w)           ; take the integer part of w
        leftover = w - ix       ; weighting factor from findex
        d_ix = long(-1*ts_diff(ix,1)) ; see where wavelength bin changes

	d_ix_0 = where(d_ix eq 0) ; indices where lambdat falls into the same lambda as previous lambdat
        d_ix_1 = where(d_ix eq 1, count1) ; indices where lambdat falls into the next lambda bi
        d_ix_2 = where(d_ix eq 2, count2) ; indices where lambdat falls into two lambda bins
        d_ix_3 = where(d_ix ge 3, count3) ; indices where lambdat falls into two lambda bins
        if count3 ne 0 then begin
            message, 'STANDARIZE_SPECTRA -- Rebinning can not handle input spectrum!!'
        endif
        
        ; indices of the rebinned coadded spectrum for those lambdat which fall uniquely into one lambda
        k = d_ix_0[where(ix[d_ix_0] ge 0 and ix[d_ix_0] le n_elements(spec)-1)]
        j = ix[k]               ; indices of the individual spectra

        rspec[i,k] = spec[j] ; for lambdat which falls into only one lambda, just copy the counts
        rivar[i,k] = ivar[j]    ; and the inverse variance

        if count1 gt 0 then begin
            ; indices of the rebinned coadded spectrum for those lambdat which straddle more than one lambda
            k = d_ix_1[where(ix[d_ix_1] ge 0 and ix[d_ix_1] le n_elements(spec)-2)]
            j = ix[k]           ; indices of the individual spectra
            w = where((ivar[j] + ivar[j+1]) ne 0d, countw)
            if countw gt 0 then begin
                k = k[where((ivar[j] + ivar[j+1]) ne 0d)] ; make sure denominators do not equal zero
                j = j[where((ivar[j] + ivar[j+1]) ne 0d)]

                ; weight the counts by ivar and lambda spacing
                rspec[i,k] = ( spec[j] * ivar[j] * (lambda[j+1] - lambdat[k]) + $
                             spec[j+1] * ivar[j+1] * (lambdat[k+1] - lambda[j+1]) ) $
                             / ( ivar[j] * (lambda[j+1] - lambdat[k]) + $
                                 ivar[j+1] * (lambdat[k+1] - lambda[j+1]) )
                ; do the same for inverse variance (weight ivar by itself)
                rivar[i,k] = ( ivar[j] * (lambda[j+1] - lambdat[k]) + $
                             ivar[j+1] * (lambdat[k+1] - lambda[j+1]))^2. $
                             / (ivar[j] * (lambda[j+1] - lambdat[k])^2. + $
                                ivar[j+1] * (lambdat[k+1] - lambda[j+1])^2.)
            endif
        endif

        if count2 gt 0 then begin
            k = d_ix_2[where(ix[d_ix_2] ge 0 and ix[d_ix_2] le n_elements(spec)-3)]
            j = ix[k]           ; indices of the individual spectra
            w = where((ivar[j] + ivar[j+1] + ivar[j+2]) ne 0d, countw)
            if countw gt 0 then begin
                k = k[w]    ; make sure denominators do not equal zero
                j = j[w]
            
                ; weight the counts by ivar and lambda spacing
                rspec[i,k] = ( spec[j] * ivar[j] * (lambda[j+1] - lambdat[k]) + $
                             spec[j+1] * ivar[j+1] * (lambda[j+2] - lambda[j+1]) + $
                             spec[j+2] * ivar[j+2] * (lambdat[k+1] - lambda[j+2]) ) $
                             / ( ivar[j] * (lambda[j+1] - lambdat[k]) + $
                                 ivar[j+1] * (lambda[j+2] - lambda[j+1]) + $
                                 ivar[j+2] * (lambdat[k+1] - lambda[j+2]))
                ; do the same for inverse variance (weight ivar by itself)
                rivar[i,k] = ( ivar[j] * (lambda[j+1] - lambdat[k]) + $
                             ivar[j+1] * (lambda[j+2] - lambda[j+1]) + $
                             ivar[j+2] * (lambdat[k+1] - lambda[j+2]))^2. $
                             / (ivar[j] * (lambda[j+1] - lambdat[k])^2. + $
                                ivar[j+1] * (lambda[j+2] - lambda[j+1])^2. + $
                                ivar[j+2] * (lambdat[k+1] - lambda[j+2])^2.)
            endif
        endif

        w = where(lambdat ge 4200 and lambdat le 4600)
        C = median(rspec[i,w])
        rspec[i,*] = rspec[i,*] / C
        rivar[i,*] = rivar[i,*] * C^2.
    endfor

    ; return a structure with wavelength, flux, and inverse variance
    raw = {lambda:lambdat, spec:rspec, ivar:rivar}
    return, raw

end


; This function coadds the spectra in the structure passed to it, but
; it only coadds the spectra referenced by the index array ix.
function Coadd, $
        raw, $		; structure: {lambda (rest Angstroms), spec (flux), ivar} all length 12000
        ix, $		; indices of objects to be coadded
        Iband=Iband

    lambda = raw.lambda     ; remove the wavelength from the structure
    spec = raw.spec(ix,*) ; weed out the spectra that haven't been selected for coaddition
    ivar = raw.ivar(ix,*)       ; do the same for the inverse variance

    spect = dblarr(12000)       ; total flux
    ivart = dblarr(12000)       ; total inverse variance

    if keyword_set(Iband) then begin
        Iband = 10.^(-Iband[ix] / 2.5)
        for k=0,n_elements(spect)-1 do begin
            if total(ivar[*,k], 1) ne 0.0 then begin
                spect[k] = total(spec[*,k] * Iband * ivar[*,k], 1) / total(Iband * ivar[*,k], 1)
                ivart[k] = (total(Iband * ivar[*,k], 1))^2. / total(Iband^2. * ivar[*,k], 1)
            endif
        endfor
    endif else begin
        ivart = total(ivar, 1) ; add up the inverse variance from each spectrum
        w = where(ivart ne 0.0) ; ignore places with zero total inverse variance
        spect[w] = total(spec[*,w] * ivar[*,w], 1) / ivart[w] ; weight each flux by its inverse variance
    endelse

    ; return a structure with wavelength, flux, and inverse variance
    ;w = uniq(lambda, sort(lambda))
    ;lambda = lambda[w]
    ;spect = spect[w]
    ;ivart = ivart[w]
    sst = {lambda:lambda, spec:spect, ivar:ivart}
    ;print, n_elements(spect), n_elements(where(spect ne 0))

    return, sst
end


; This function gets an entire spectrum (red-side and blue-side) based
; on an input filename.
function ReadSpec, file
    bs = mrdfits(file,1, /silent)
    rs = mrdfits(file,2, /silent)
    ss = {ss, lambda:[bs.lambda, rs.lambda], spec:[bs.spec, rs.spec], ivar:[bs.ivar, rs.ivar]}
    return, ss
end

; This program coadds and plots various subsets of serendip spectra.
pro Perform_Coadd
    Lya = 1215.68
    Ha = 6562.8
    Hb = 4861.3
    OIII = 5006.9

    sigma = 1.17741 * 2

    numread = 2             ; number of columns to read

    ; read in the specifications for each serendip
    readcol, 'sngl_all_cluster_for_coadd.dat',  format='A,I,D,D,A,A', mask, slit, Iband, z, q, file, skipline=2, numline=numread
    file = '/Volumes/Data2/orelse/lemaux/deimos/sc1604/'      + file

    lambdas = dblarr(numread,8192) ; 8192-element wavelength grid for all serendips
    specs = dblarr(numread,8192) ; 8192-element flux grid for all serendips 
    ivars = dblarr(numread,8192) ; 8192-element inverse variance grid for all serendips
	
    for i=0,numread-1 do begin  ; for each serendip
        ss = ReadSpec(file(i))  ; read in its
        lambdas(i,*) = ss.lambda ; wavelength grid,
        specs(i,*) = ss.spec    ; flux grid,
        ivars(i,*) = ss.ivar    ; and inverse variance grid
    endfor
	
    ; blueshift and rebin the spectrum for each serendip
    rawLya = Standardize_Spectra(lambdas, specs, ivars, z)
    ;rawHa = Standardize_Spectra(Ha, lambdas, specs, ivars, serlambda)
    ;rawHb = Standardize_Spectra(Hb, lambdas, specs, ivars, serlambda)
    ;rawOIII = Standardize_Spectra(OIII, lambdas, specs, ivars, serlambda)

    ;use subsets 1, 3, and 4 plus the nonLya

    subset1 = indgen(numread)

    finalspec = Coadd(rawLya, subset1, Iband=Iband) ; call the coaddition routine
    smooth_ivar, finalspec.spec, finalspec.ivar, sigma, spec, ivar
    finalspec.spec = spec
    finalspec.ivar = ivar
    mwrfits, finalspec, 'fill_gap_2galaxy_test.fits'
    device, retain=2
    splot, finalspec.lambda, finalspec.spec
end
