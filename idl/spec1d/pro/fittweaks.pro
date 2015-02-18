pro fittweaks,maskname
;+
; NAME:
;
;  FITTWEAKS
;
; PURPOSE:
;
;  Fit wavelength offsets vs. pixel over whole mask, using SKYTWEAK_1D
;
; CATEGORY:
;
;  spec1d
;
; CALLING SEQUENCE:
;
;  fittweaks,maskname
;
; INPUTS:
;  masknum -- mask name (e.g. 1251, not assumed to be string but can be)
;
; OUTPUTS:
;   file named skytweaks.MASKNUM.fits, where MASKNUM is the string
;   version of masknum.  This file includes all the factors used by
;   TWEAK_SPEC.
;
; RESTRICTIONS
; 
; ONLY TESTED FOR DEEP2 (1200/7800) DATA!
;
; MODIFICATION HISTORY:
;  JAN 5jan05
;-


if n_elements(masknum) eq 0 then begin
    cd,'.',current=currdir
    broken=str_sep(currdir,'/')
    broken=broken[where(broken ne '')]
    masknum=broken[n_elements(broken)-2]

    print,'I am assuming that the mask label is: ',masknum

endif



files=findfile('spec1d*.fits')
files=files[where(strpos(files,'.s') lt 0)]

nfiles=n_elements(files)



if n_elements(degree) eq 0 then degree=3




fitb=fltarr(nfiles,degree)
fitr=fitb
xmm=fltarr(nfiles)
ymm=fltarr(nfiles)
sigmab=xmm
sigmar=sigmab
slitno=intarr(nfiles)
simple_tables,'*bin*',slitnames=names,slitwid=widths
goodslit=slitno
numb=slitno
numr=slitno
errarr=slitno
dlamb=xmm
dlamr=xmm

for i=0,nfiles-1 do begin 
	s=mrdfits(files[i],1,hdr,/sil) 
	xmm[i]=sxpar(hdr,'XMM') 
	ymm[i]=sxpar(hdr,'YMM') 
	slitno[i]=sxpar(hdr,'SLITNO') 
	fitb[i,*]=skytweak_1d(s,degree=degree,template=stemp_flux,$
		temp_wave=temp_wave,sigma=sig,shift=shift,err=errb)

        if slitno[i] lt 10 then slitnum='00'+string(slitno[i],format='(i1)') $
          else if slitno[i] lt 100 then $
             slitnum='0'+string(slitno[i],format='(i2)') $
          else slitnum=string(slitno[i],format='(i3)') 
   
        calfile=findfile('calibSlit.*.'+slitnum+'B*')	;Changed syntax on 1/11/10 BL, cuz crashing on slit 005 for mask XL005B, not sure why it didn't for XL005A
        calstruct=mrdfits(calfile[0],1,/sil)
        slitfile=findfile('slit.*.'+slitnum+'B*')
        slitstruct=mrdfits(slitfile[0],1,/sil)

        if strlen(calfile) gt 0 then begin
            nobs=n_elements(slitstruct.dlam)/n_elements(calstruct.dlam)
            whgood=where(calstruct.dlam NE 0,ct)
            if ct gt 0 then dlamb[i]=mean((total(slitstruct.dlam,2)/nobs-calstruct.dlam)[whgood]) else dlamb[i]=0.
        endif

; get dlam here
	numb[i]=n_elements(shift) 
	sigmab[i]=sig 
	s=mrdfits(files[i],2,hdr,/sil) 
	fitr[i,*]=skytweak_1d(s,degree=degree,template=stemp_flux,$
		temp_wave=temp_wave,sigma=sig,shift=shift,err=errr) 

        calfile=findfile('calibSlit.*.'+slitnum+'R*')
        calstruct=mrdfits(calfile[0],1,/sil)
        slitfile=findfile('slit.*.'+slitnum+'R*')
        slitstruct=mrdfits(slitfile[0],1,/sil)

        if strlen(calfile) gt 0 then begin
            nobs=n_elements(slitstruct.dlam)/n_elements(calstruct.dlam)
            whgood=where(calstruct.dlam NE 0,ct)
            if ct gt 0 then dlamr[i]=mean((total(slitstruct.dlam,2)/nobs-calstruct.dlam)[whgood]) else dlamr[i]=0.
        endif





	numr[i]=n_elements(shift) 
	sigmar[i]=sig 
        errarr[i] = errb>errr 
	wh=where(fix(names) eq slitno[i],ct) 
	if ct eq 0 then goodslit[i]=0 else $
		goodslit[i]=widths[wh[0]] gt 0.745 AND $
			widths[wh[0]] lt 1.2 
    endfor

    sigdlamb=djsig(dlamb-shift(dlamb,1))/sqrt(2.)
    sigdlamr=djsig(dlamr-shift(dlamr,1))/sqrt(2.)

    baddlamb=abs(dlamb-djs_median(dlamb,width=5,boundary='reflect')) $
      gt 3*sigdlamb

    baddlamr=abs(dlamr-djs_median(dlamr,width=5,boundary='reflect')) $
      gt 3*sigdlamr



goodslit=goodslit AND fitb[*,0] NE 0. AND fitr[*,0] NE 0.AND $
  (baddlamr eq 0) AND (baddlamb eq 0)

;filename=strcompress('outputs'+string(degree)+'.sav',/remove)

;save,f=filename


; do fits as a function of deimos x/y
deimos_degree=5

xl=2*(xmm-min(xmm))/(max(xmm)-min(xmm))-1
yl=2*(ymm-min(ymm))/(max(ymm)-min(ymm))-1
minx=min(xmm)
miny=min(ymm)
maxx=max(xmm)
maxy=max(ymm)


infostruct={slitno:slitno,xmm:xmm,ymm:ymm,dlamb:dlamb,dlamr:dlamr,sigmar:sigmar,sigmab:sigmab,errarr:errarr, goodslit:goodslit}

   wh=where(goodslit and sigmar lt 3*median(sigmar) and sigmab lt 3*median(sigmab) and errarr eq 0,comp=comp, goodct)

fitcoeff_xr=dblarr(degree, deimos_degree)
fitcoeff_yr=fitcoeff_xr
fitcoeff_xb=fitcoeff_xr
fitcoeff_yb=fitcoeff_xr
   sig=sigmar

if goodct gt 0 then begin
    for element=0,degree-1 do begin

        fit=fitr[*,element] 
        fitpars=svdfit(xl[wh],fit[wh],deimos_degree,meas=3*sigmar[wh]+.01, $
                  yfit=yfit_x,/legen,/doub)  
        fitvsy=svdfit(yl[wh],fit[wh]-yfit_x,deimos_degree,meas=3*sigmar[wh]+.01, $
                 yfit=yfit_y,/legen,/doub) 
        fitcoeff_xr[element,*]=fitpars 
        fitcoeff_yr[element,*]=fitvsy 

        fit=fitb[*,element] 
        fitpars=svdfit(xl[wh],fit[wh],deimos_degree,meas=3*sigmab[wh]+.01, $
                  yfit=yfit_x,/legen,/doub)  
        fitvsy=svdfit(yl[wh],fit[wh]-yfit_x,deimos_degree,meas=3*sigmab[wh]+.01, $
                 yfit=yfit_y,/legen,/doub) 
        fitcoeff_xb[element,*]=fitpars 
        fitcoeff_yb[element,*]=fitvsy 
    endfor

mkhdr,hdr,fitcoeff_xb
sxaddpar,hdr,'MINX',minx
sxaddpar,hdr,'MAXX',maxx
sxaddpar,hdr,'MINY',miny
sxaddpar,hdr,'MAXY',maxy
sxaddpar,hdr,'DEGFIT',degree,'Degree used to fit delta-lambda vs. pix, +1'
sxaddpar,hdr,'DEGXY',deimos_degree,'Degree of fit for coeffs vs. DEIMOS X/Y, +1'
sxaddhist,'All fits are for Legendre polynomials',hdr,/COMMENT
sxaddhist,'mapping [0,4096] and [XMIN,XMAX] and [YMIN,YMAX] to [-1,1]',hdr,/COMMENT
sxaddhist,'Add the polynomial to original wavelength solution',hdr,/COMMENT
sxaddhist,'to get true, air wavelength',hdr,/COMMENT

mwrfits,fitcoeff_xb,strcompress('skytweaks.'+string(masknum)+'.fits',/remove),hdr,/create
mwrfits,fitcoeff_yb,strcompress('skytweaks.'+string(masknum)+'.fits',/remove),hdr
mwrfits,fitcoeff_xr,strcompress('skytweaks.'+string(masknum)+'.fits',/remove),hdr
mwrfits,fitcoeff_yr,strcompress('skytweaks.'+string(masknum)+'.fits',/remove),hdr

mwrfits,infostruct,strcompress('skytweaks.'+string(masknum)+'.fits',/remove)

endif else message,'No good slits!',/info

return
end
