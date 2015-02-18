function applytweaks,lambda,header,path 

;+
; NAME:
;
;  APPLYTWEAKS
;
; PURPOSE:
;
;  Apply wavelength tweaks from SKYTWEAK_1D skyline fits over mask
;
; CATEGORY:
;
;  spec1d
;
; CALLING SEQUENCE:
;
;   new_lambda=APPLYTWEAKS(oldlambda,header,path [FITBX=,FITBY=, $
;               FITRX=,FITRY=, MINX=,MINY=,MAXX=,MAXY=])
;
; INPUTS:
;   oldlambda -- original wavelength array
;   header    -- a FITS header from this spec1d file (any will do)
;   path      -- path to the mask's directory (assumed to be the
;                current directory if not set), containing the
;                skytweaks file
;
; OUTPUTS:
;   newlambda -- tweaked wavelength array.
;
; RESTRICTIONS:
;
;  ONLY TESTED FOR DEEP2 (1200/7800) DATA!
;
; MODIFICATION HISTORY:
;  JAN 5jan05
;-

if n_elements(path) eq 0 then path=''

skytweak_file=(findfile(concat_dir(path,'skytweaks*.fits' ) ))[0]

if strlen(skytweak_file) gt 2 then begin

    fitbx=mrdfits(skytweak_file,0,hdr,/sil)
    fitby=mrdfits(skytweak_file,1,hdr,/sil)
    fitrx=mrdfits(skytweak_file,2,hdr,/sil)
    fitry=mrdfits(skytweak_file,3,hdr,/sil)

    minx=sxpar(hdr,'MINX')
    miny=sxpar(hdr,'MINY')
    maxx=sxpar(hdr,'MAXX')
    maxy=sxpar(hdr,'MAXY')

    if max(lambda) lt 8500. then begin
        fit_x=fitbx
        fit_y=fitby
    endif else begin
        fit_x=fitrx
        fit_y=fitry
    endelse

    xmm=sxpar(header,'XMM')
    ymm=sxpar(header,'YMM')

; legendre values
    xl=2*(xmm-minx)/(maxx-minx)-1.
    yl=2*(ymm-miny)/(maxy-miny)-1.

    xpix=2.*dindgen(n_elements(lambda))/n_elements(lambda)-1.

    s=size(fit_x,/dim)
    ncoeff=s[0]

    coeffs=dblarr(ncoeff)

    for i=0,ncoeff-1 do coeffs[i]=polyleg(xl,fit_x[i,*]) + $
                     polyleg(yl,fit_y[i,*]) 


    return,lambda+polyleg(xpix,coeffs)
endif else begin
    print,'No skytweak file found!'
    return,lambda
endelse


end
