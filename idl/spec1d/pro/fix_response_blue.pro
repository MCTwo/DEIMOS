;+ 
; NAME:
;      fix_response
;
; PURPOSE:
;      Make the throughput correction to the 1-D spectra.
; 
; CALLING SEQUENCE:
;      fix_response,spec1d
; INPUTS:
;      spec1d - ss1d structure containing spec,lambda,ivar tags. 
;
; OUTPUTS:
;      spec1d - spec and ivar tags are changed. 
;
; ROUTINES CALLED:
;      deimos_correction 
; 
; MODIFICATION HISTORY:
;  ry 16jun2004
;-

pro fix_response_blue,spec1d

  pixelwidth = (shift(spec1d.lambda,-1) - shift(spec1d.lambda,1))/2.0
  pixelwidth[0] = pixelwidth[1]
  num = n_elements(pixelwidth)
  pixelwidth[num-1] = pixelwidth[num-2]
  if min(pixelwidth) gt 0 and min(pixelwidth) lt 0.4 then begin
    inv_throughput = deimos_correction_masterblue(spec1d.lambda)*pixelwidth
    spec1d.spec = spec1d.spec*inv_throughput
    spec1d.ivar = spec1d.ivar/(inv_throughput*inv_throughput)
  endif
end
