;+
; NAME:
;	ytoradec
;
; PURPOSE:
;	To find the RA, Dec of a position given the mask number,
;       slit number, and y offset in pixels from the bottom of the
;       slit as seen in zspec.
;	
;       
; CALLING SEQUENCE:
;	ytoradec, mask, slit, y, ra, dec
;
; INPUT:
;	mask    four-digit mask number
;       slit    one- to three-digit slit number
;       y       y spatial offset from slit bottom in pixels
;
; OUTPUTS:
;       ra      RA of that position in degrees
;       dec     Dec of that position in degrees
;
; PROCEDURE:
;	
;       
; REVISION HISTORY:
;       enk 2006-05-28
;-

pro ytoradec, mask, objno, y, ra, dec, delta=delta
    mn = strcompress(string(mask),/rem)
    filename = getenv('D2_RESULTS')+'/../1HSmasks/' + mn + '/1HSmask.' + mn + '.fits'
    if file_test(filename) then begin
        ip = mrdfits(filename, 1, /silent)
        mm = mrdfits(filename, 2, /silent)
        mm = mm[where(mm.sg eq 'G' and long(mm.objno) eq long(objno))]
    endif else begin
        filename = getenv('D2_RESULTS')+'/../1HSmasks/' + mn + '/1HSmask.' + mn + 'SN.fits'
        if file_test(filename) then begin
            ip = mrdfits(filename, 1, /silent)
            mm = mrdfits(filename, 2, /silent)
            mm = mm[where(mm.sg eq 'G' and long(mm.objno) eq long(objno))]
            print, objno
        endif else message, 'Mask file not found.'
    endelse

    slit_pa = (ip.pa + mm.slitpa) * !DTOR                                                      ;slit PA (rad)
    slit_pa_sin = sin(slit_pa)
    slit_pa_cos = cos(slit_pa)

    obj_ra = mm.ra                                                                             ;target RA (deg)
    obj_dec = mm.dec                                                                           ;target Dec (deg)
    obj_dec_cos = cos(obj_dec * !DTOR)

    slit_ra = obj_ra + 0.5d * (mm.topdist-mm.botdist)/3600d * slit_pa_sin/obj_dec_cos          ;slit center RA (deg)
    slit_dec = obj_dec + 0.5d * (mm.topdist-mm.botdist)/3600d * slit_pa_cos                    ;slit center Dec (deg)
    slit_dec_cos = cos(slit_dec * !DTOR)

    slit_length = mm.topdist+mm.botdist                                                        ;arcsec
    slit_width = mm.slitwidth                                                                  ;arcsec
    half_slit_length = 0.5d * slit_length / 3600d                                              ;deg
    half_slit_width = 0.5d * slit_width / 3600d                                                ;deg 

    deimosplatescale = 0.119d                                                                  ;arcsec/pixel
    y_deg = y * deimosplatescale / 3600d                                                       ;y offset from slit bottom (deg)

    if keyword_set(delta) then begin
        ra = obj_ra + (y_deg * slit_pa_sin / slit_dec_cos) ;RA of y offset (deg)
        dec = obj_dec + (y_deg * slit_pa_cos)    ;Dec of y offset (deg)        
    endif else begin
        y0_ra = slit_ra - (half_slit_length * slit_pa_sin / slit_dec_cos) ;RA of bottom center of slit (deg)
        y0_dec = slit_dec - (half_slit_length * slit_pa_cos) ;Dec of bottom center of slit (deg)
        
        ra = y0_ra + (y_deg * slit_pa_sin / slit_dec_cos) ;RA of y offset (deg)
        dec = y0_dec + (y_deg * slit_pa_cos)    ;Dec of y offset (deg)
    endelse
end
