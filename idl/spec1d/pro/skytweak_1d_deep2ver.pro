;+
; NAME:
;   skytweak_1d
;
; PURPOSE: 
;   Tweaks the wavelength solution of DEEP2 DEIMOS 1d spectra by
;   cross-correlating LOCALLY with a template
;   sky emission spectrum, and fits the local shifts to a
;   polynomial of order degree-1 (decomposed into Legendre polynomials).
;   Based on deimos_skytweak_var in the spec2d pipeline.
;    This routine does NOT alter the original structure.
;
; CALLING SEQUENCE:
;   wave_fit = skytweak_1d(slit, [degree=,new_wave_grid=,$
;                                 old_wave_grid= , flux1d= , $
;                                 template= , temp_wave=,fitlam= , $
;                                 shift= , serr=,chisq=, $
;                                 errflag=,sigma=,/PLOT])
;
; INPUTS:
;   slit       - untweaked DEEP2 spec1d structure 
;
; OUTPUTS:
;   wave_fit   - Coefficients of Legendre polynomial fit.
;
; KEYWORDS:
;   /PLOT      - set to output diagnostic numbers/plots
;   degree     - # of coefficients of polynomial to fit
;
; OPTIONAL OUTPUTS:
;   new_wave_grid - 1-d tweaked wavelength solution in a regular 0.2A
;                grid.
;   old_wave_grid - 1-d untweaked wavelength solution, as above.
;   flux1d     - 1-d sky flux spectrum, corresponding to new_wave1d and,
;                originally, old_wave1d
;   template   - full smoothed template spectrum (passable for efficiency)
;   temp_wave  - wavelength array for template spectrum (ditto)
;   fitlam     - wavelengths at which the new wavelength solution was
;                fit.
;   shift      - amount by which old wavelength solution was shifted
;                at the fitlam, before refitting.
;   serr       - estimated errors in shift
;   chisq      - chi-square of new fit for 2-d wavelength array
;
;   errflag    - 1 if routine encountered an error; zero otherwise.
;   sigma      - DJSIG of residuals from fit
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
;
; REVISION HISTORY:
;   deimos_skytweak_var Written by BFG Sept. 02
;   This routine written around the skeleton of that one by BFG in Dec
;   04.
;   and modified in small ways by JAN, Jan 04
;----------------------------------------------------------------------


function skytweak_1d,  slit, new_wave_grid=new_wave_grid, $
               old_wave_grid=fullwave, flux1d=fullflux, $
               template=temp_flux, temp_wave=temp_wave,$
                       fitlam=fitlam, shift=shift, $
               serr=shifterr, chisq = chisq2,  errflag=errflag,$
               sigma=sigma,plot=plot,degree=degree

if n_elements(plot) eq 0 then plot = 0
  if n_elements(degree) eq 0 then degree = 3
  errflag =  0
;restore,  file='../etc/template.sav'

tmp=minmax(slit.lambda)

if n_elements(temp_flux) eq 0 OR n_elements(temp_wave) eq 0 then begin


  dir = getenv('DEEP_DIR')+'/'
  if dir eq '/' then message, 'you must set $DEEP_DIR'
  restore,  file=dir+'spec2d/etc/template.sav'

;----restores template sky spectrum from HiRes: 
;----wavelength array in temp_wave, flux in temp_flux

;----smooth the HiRes spectrum with a gaussian to make similar to
;----DEIMOS resolution

  sigma=9
  halfwidth=15
  kernel=findgen(2*halfwidth+1)-halfwidth
  kernel=exp(-kernel^2/2/sigma^2)
  kernel=kernel/total(kernel)

  temp_flux =  convol(temp_flux,  kernel,  /center,/edge_wrap)
endif

  sizex = n_elements(slit.lambda)

  whgood=where((slit.infomask AND 1b) eq 0b,goodct)


;---set up 0.2A wavelength grid
 if goodct gt 2 then begin 
     minlambda =  min(slit.lambda[whgood], max=max1) > min(temp_wave, max=max2)
     maxlambda =  max1 < max2
 endif else begin
     minlambda=0
     maxlambda=0
 endelse



  dlam =  0.2  ; size of shifts in cross correlation
  oversample=0.2/dlam

  nlag = 4*oversample*2+1   ; number of shifts in cross correlation

  minlambda = (minlambda+1) < 1.5E4
  maxlambda = maxlambda-1 >500      ; kill off spurious end effects


  npts =  floor((maxlambda-minlambda)/dlam)
  fullwave = findgen(npts>1)*dlam + minlambda

;---interpolate observed sky spectrum onto grid
  iwave = where((slit.lambda ge minlambda) and (slit.lambda le maxlambda),ct)
  if ct gt 0 then fullflux = interpol(slit.skyspec[iwave], slit.lambda[iwave],  fullwave) else fullflux=0.

  badspec=median(slit.nbadpix) gt 0.3*(slit.r2-slit.r1)

;---test whether bspline is useful:
  if total(fullflux) gt 0. then begin

;---interpolate template onto wavelength grid
  iwave = where((temp_wave ge minlambda) and (temp_wave le maxlambda))
  stemp_flux = interpol(temp_flux[iwave], temp_wave[iwave],  fullwave)


  window =  100.                ; 100A windows
  nregions =  fix(2*(maxlambda-minlambda)/window)

 
  if maxlambda gt 8500. then color='R' else color='B'

  if color eq 'B' then begin 
     sig_thresh =  400.         ; minimum signal required in stemp_flux over window
     peakthresh =  1.5          ; minimum height for a significant peak
  endif
  if color eq 'R' then begin 
     sig_thresh= 400.           ; minimum signal required in stemp_wave over window
     peakthresh = 10.           ;minimum signal required for a significant peak
  endif
  lag =  indgen(nlag) - nlag/2  ;shift values for x-correlation
  wavlag = lag*dlam
  npoly = 3
  cc =  fltarr(nlag, nregions)     ; cross-correlation values for lag
  cc_err =  cc                  ;  estimated error in cross-correlation

  shiftfit =  fltarr(npoly, nregions) ;polynomial coefficients for fit to cc
  fiterr = shiftfit             ; errors on coefficients

  shift =  fltarr(nregions)     ;peak of fit
  shifterr =  shift             ;error in peak finding

  dofit = shift                 ;which regions' shifts have enough signal to trust
  fitlam = shift                ;central wavelength of fit regions
  fitpix = shift                ; central pixel of fit regions
  tmed = shift                  ;median of each window

  for i=0, nregions-1 do begin
;---march across spectrum in 100A windows, overlapping by 50A
     minwave =  minlambda+i*window/2.+nlag/2*dlam+dlam
     maxwave =  minwave + window
     if maxwave ge maxlambda then begin
        minwave = maxlambda-window-nlag/2*dlam-dlam
        maxwave = maxlambda
     endif


     iwin = where((fullwave ge minwave) and (fullwave le maxwave))
     fitpix[i] = min(iwin)+(max(iwin) - min(iwin))/2.
     wave = fullwave[iwin]
     fitlam[i] = fullwave[fitpix[i]]
     flux = fullflux[iwin]
     tflux = stemp_flux[iwin]
     djs_iterstat,  tflux,  median=tmp
     tmed[i] = tmp
     
;----does the input bspline exist here?
     if total(flux) ne 0. then begin    
;----is there enough signal to do cross-correlation?
        if total(tflux-tmed[i]) ge sig_thresh then begin
           

;----do cross-correlation, fit with a polynomial, find max
           cc[*, i] = c_correlate(flux,  tflux,  lag)
     
           wpeak =  where(tflux gt peakthresh) ;where there's a significant peak
           if wpeak[0] ne -1 then error = .2/(sqrt(total(tflux[wpeak])/2.)) $
              else error = 1.

           ;if wpeak[0] ne -1 then print,total(tflux-tmed[i]),((total(tflux[wpeak])/2.))
           cc_err[*, i] = replicate(error,  nlag)

; ----only fit near the xcorr peak
           maxcorr = max(cc[*, i], maxpix)
           ; number of pix away from peak to go
           nout=2
           mintofit = (maxpix-nout*oversample) > 0
           if maxpix lt n_elements(wavlag) - nout*oversample -1 then $
             maxtofit = mintofit+2*nout*oversample else begin
              maxtofit = (maxpix+nout*oversample) < (n_elements(wavlag)-1)
              mintofit = maxtofit-2*nout*oversample
           endelse

           ;print, minwave,  maxwave,  cc_err[0]
           shiftfit[*, i] =  poly_fit(wavlag[mintofit:maxtofit], $
                               cc[mintofit:maxtofit, i], npoly-1, $
                               measure_errors=cc_err[mintofit:maxtofit, i], $
                               sigma=err)
           fiterr[*, i] =  err
           shift[i] = - 0.5*shiftfit[1, i]/shiftfit[2, i]
   
           shifterr[i] = (abs(shift[i])>0.005)*sqrt((fiterr[1, i]/shiftfit[1, i])^2 +$
                                            (fiterr[2, i]/shiftfit[2, i])^2)
;---deweight if we're in the noisy part of the template.
           if minwave ge 8950. then shifterr[i] = sqrt(5.)*shifterr[i]        
           
           dofit[i] = 1 
        endif else begin
           shifterr[i] =  1d10
           dofit[i] = 0
        endelse
     endif else begin
        shifterr[i] =  1d10
        dofit[i] = 0
     endelse
 endfor

  whfit =  where(dofit)
  nfit = n_elements(whfit)
;whfit =  whfit[1:nfit-2];drop first and last points.
  fitlam = fitlam[whfit]
  shift = shift[whfit]
  shifterr = shifterr[whfit]
  fitpix = fitpix[whfit]
;  cent_row = floor(sizey/2.)
;  cent_wave = wave2d[*, cent_row]
  pix = findgen(sizex)

; find interpolated pixel numbers corresponding 
; to fit wavelengths
  fitpix2 = fitlam
  fitpix2 = interpol(pix, slit.lambda, fitlam) 


; the new wavelengths for the fit pixels:
  shift_wave = shift;fitlam+shift


  npix = n_elements(fullwave)
  npix2 = sizex

shifterr=sqrt(shifterr^2+(median(shifterr) > 0.0025)^2)

; do new fits for regular grid and input wavelength array
  xx = fitpix2/(npix2/2.) -1

  pxx2 =  fitpix2

; fit vs. wavelength
;for degree = 2,3 do begin
  wave_fit2 = svdfit(xx, shift_wave, degree, /double, /legendre, $
             yfit=pxx2, measure_errors = shifterr, chisq=chisq,sigma=sigma)
;  print,degree,chisq,djsig(shift_wave-pxx2)
;endfor


  chisq2 = chisq
  pxx =  fitpix



; fit vs. pixel in wavelength array
; this fit generally looks better, but makes less sense for global fits
;for degree = 2,4 do begin
  wave_fit = svdfit(fitpix/(npix/2.)-1, shift_wave, degree,/double,/legendre, $
          yfit=pxx, measure_errors = shifterr, chisq=chisq,sigma=sigma)  
;endfor
;print,chisq
  

; evaluate new 1-d fit
  pixels = indgen(npix2)
  new_wave =  polyleg(pixels/(npix2/2.) -1 ,  wave_fit)
new_wave2 =  polyleg(pixels/(npix2/2.) -1 ,  wave_fit2)

  if plot then begin
      loadct,12
      plot,fitlam,shift,psym=4,yrange=[min(shift)-0.02,max(shift)+0.02], $
        xtit="Wavelength ("+ STRING("305B)+")", $
        ytit="Cross-correlation shift ("+ STRING("305B)+")",xthick=3,ythick=3
      oplot,fitlam,shift,psym=4,color=16*12,thick=3
!p.thick=3
      oploterr,fitlam,shift,shifterr
!p.thick=4
      oplot,slit.lambda,new_wave,line=0
;      oplot,slit.lambda,new_wave2,line=1

      print,wave_fit
      print,djsig(pxx-shift),n_elements(shift)
  endif

  sigma=djsig(pxx-shift)

  pixels = indgen(npix)
  new_wave_grid =  polyleg(pixels/(npix/2.) -1 ,  wave_fit)

  if badspec then begin
      print,'>1/3 of spectrum is bad - beware!'
      print,'ERRFLAG set to 2'
      errflag=2
  endif
  
  return, wave_fit2
  endif else begin
    print,  'WARNING: no useful bspline found in deimos_skytweak_var.  Wavelength solution not tweaked'
    errflag =  1
    return, fltarr(degree) 
  endelse
end


