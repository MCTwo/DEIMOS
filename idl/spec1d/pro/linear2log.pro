;+
;-----------------------------------------------------------------------------
; NAME
;      linear2log.pro
;-----------------------------------------------------------------------------
; PURPOSE
;      Converts a spectrum in linear lambda to log lambda. 
;
;-----------------------------------------------------------------------------
; SYNTAX
;      loglam = linear2log(spec1d, [binsize=binsize, flux=flux, ivar=ivar])
;-----------------------------------------------------------------------------
; INPUTS
;      spec1d = a structure from a spec1d output file. 
;-----------------------------------------------------------------------------
; KEYWORDS
;      
;-----------------------------------------------------------------------------
; OUTPUTS
;
;-----------------------------------------------------------------------------
; PROCEDURES CALLED 
;      makearr
;-----------------------------------------------------------------------------
; EXAMPLES
;
;-----------------------------------------------------------------------------
; COMMENTS
;
;-----------------------------------------------------------------------------
; HISTORY
;      Created August 7, 2002 by mcc.
;-----------------------------------------------------------------------------
;-

function linear2log, spec1d, binsize=binsize, flux=flux, ivar=ivar

;;; EXTRACT WAVELENGTH FROM spec1d STRUCTURE.
  lambda = double(spec1d.lambda)

;;; DETERMINE THE LENGTH OF THE WAVELENGTH VECTOR.
  nlam = n_elements(lambda)

;;; CHECK IF binsize ARGUMENT WAS PASSED. IF NOT, USE DEFAULT binsize
;;; OF 0.0001 IN LOG LAMBDA.
  if N_ELEMENTS(binsize) gt 0 then binsize = binsize[0] $
    else binsize = 0.0001

;;; DETERMINE THE LIMITS OF THE LOG LAMBDA VECTOR. 
  logwv = ALOG10(lambda)
  lims = MINMAX(logwv)
;;; DETERMINE THE NUMBER OF BINS NEEDED. HERE WE ARE ASSUMING THAT THE
;;; lambda ARGUMENT WAS CONTINUOUS IN WAVELENGTH.
  nloglam = floor((lims[1] - lims[0]) / binsize)

  loglam = lims[0]+binsize*dindgen(nloglam)

; center loglam amidst original wavelength array
loglam=loglam+mean(lims)-(loglam[0]+loglam[nloglam-1])/2.

;makearr(nloglam, lims[0], lims[1],/doub)

;;; INTERPOLATE TO GET OBJECT FLUX AND IVERSE VARIANCE.
  objflux = spec1d.spec
  objivar = spec1d.ivar
  ratio = n_elements(objflux)/nloglam
; smooth data if it is to be compressed
  if ratio gt 1.5 then begin
      print, 'reducing resolution of data by factor: ', ratio
      smoothfactor = round(ratio)
      objflux_out = ivarsmooth(objflux, objivar, smoothfactor, objivar_out)
      objflux = objflux_out
      objivar = objivar_out
  endif

  flux = INTERPOL(objflux, logwv, loglam, /SPLINE)
  ivar = INTERPOL(objivar, logwv, loglam, /SPLINE)

return, loglam

end











