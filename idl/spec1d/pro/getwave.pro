;**************************
pro getwave,h,wave
;+
; PURPOSE:
;	Create a wavelength vector from a FITS header of a spectrum
;	The FITS header must contain a CRVAL1 keyword and a CRPIX1 or
;	CD1_1 keyword.
;
; CALLING SEQUENCE:
;	GETWAVE, hdr, wave
;
; INPUTS:
;	hdr - FITS header for a 1-d spectrum, string array
;
; OUTPUTS:
;	wave - Wavelength vector.    
;-

 zparcheck, 'GETWAVE', h, 1, 7, 1, 'FITS header'

 xsize = sxpar(h,'NAXIS1')
 if xsize LE 0 then message, $
	'ERROR - FITS header does not contain positive value of NAXIS1'
 w0 = sxpar(h,'CRVAL1')
 if !ERR EQ -1 then message,'FITS header does not contain CRVAL1 value'
 wdelt = sxpar(h,'CDELT1')
 if !ERR EQ -1 then begin
	wdelt = sxpar(h,'CD1_1')
        if !ERR EQ -1 then message, $
             'FITS header does not contain CDELT1 or CD1_1 value'
 endif
 
 wave = w0 + wdelt*findgen(xsize)

 return
 end














