;+
;  The purpose of this routine is to calculate the
;  redshift-completeness of each mask in a consistent manner for all
;  the masks observed so far.  The output is an ascii file containing
;  the maskname and the redshift completeness.
;
;  Note:  the output file 'zcomp.dat' will be written to the current directory
;-

pro zcomp, maskname, success, bluesuccess, redsuccess

  zcat = MRDFITS('/deep0/marc/deep/deep2products/zcat.latest.fits',1) ;new catalog
  tmp = zcat.maskname

  masks = tmp[sort(tmp)]
  masks = masks[uniq(masks)]
  
  bluelim = 0.4
  redlim = 0.5

  nmasks = n_elements(masks)
  nblue = fltarr(nmasks)

  nobj = nblue
  nbluegood = nblue
  nbluebad = nblue
  nredgood = nblue
  nredbad = nblue
  nred = nblue
  ngood = nblue
  nbad = nblue

; count blue objects:R-I<bluelim, PGAL>0.2
  for i = 0, n_elements(masks) -1 do $
    nblue[i] = total(zcat.magr-zcat.magi lt bluelim $
                     AND zcat.pgal gt 0.2 $
                     AND zcat.maskname eq masks[i])

; count good blue objects: Z>~0, ZQ=3 or 4
 for i = 0, n_elements(masks) -1 do $
    nbluegood[i] = total(zcat.magr-zcat.magi lt bluelim $
                     AND zcat.pgal gt 0.2 $
                     AND zcat.maskname eq masks[i] $
                        AND zcat.zquality ge 3 AND zcat.zquality le 4 $
                        AND zcat.z gt 0.005)

; count bad blue objects: ZQ=1 or 2
 for i = 0, n_elements(masks) -1 do $
    nbluebad[i] = total(zcat.magr-zcat.magi lt bluelim $
                     AND zcat.pgal gt 0.2 $
                     AND zcat.maskname eq masks[i] $
                        AND zcat.zquality ge 1 AND zcat.zquality le 2)

for i = 0, n_elements(masks) -1 do $
    nredgood[i] = total(zcat.magr-zcat.magi gt redlim $
                     AND zcat.pgal gt 0.2 $
                     AND zcat.maskname eq masks[i] $
                        AND zcat.zquality ge 3 AND zcat.zquality le 4 $
                        AND zcat.z gt 0.005)

for i = 0, n_elements(masks) -1 do $
    nred[i] = total(zcat.magr-zcat.magi gt redlim $
                     AND zcat.pgal gt 0.2 $
                     AND zcat.maskname eq masks[i])


 for i = 0, n_elements(masks) -1 do $
    nredbad[i] = total(zcat.magr-zcat.magi gt redlim $
                     AND zcat.pgal gt 0.2 $
                     AND zcat.maskname eq masks[i] $
                        AND zcat.zquality ge 1 AND zcat.zquality le 2)

for i = 0, n_elements(masks) -1 do $
    ngood[i] = total( zcat.pgal gt 0.2 $
                     AND zcat.maskname eq masks[i] $
                        AND zcat.zquality ge 3 AND zcat.zquality le 4 $
                        AND zcat.z gt 0.005)

 for i = 0, n_elements(masks) -1 do $
    nbad[i] = total( zcat.pgal gt 0.2 $
                     AND zcat.maskname eq masks[i] $
                        AND zcat.zquality ge 1 AND zcat.zquality le 2)



 for i = 0, n_elements(masks) -1 do $
    nobj[i] = total(zcat.pgal gt 0.2 $
                     AND zcat.maskname eq masks[i])


; print, 'blue fraction: ', median(nblue/nobj), ' +/- ', djsig(nblue/nobj)
; print, 'red fraction: ', median(nred/nobj), ' +/- ', djsig(nred/nobj)
; print, 'blue success: ', median(nbluegood/(nbluebad+nbluegood))
 bluesuccess = (nbluegood/(nbluebad+nbluegood))
 redsuccess = (nredgood/(nredbad+nredgood))
; print, 'red success:', median(nredgood/(nredbad+nredgood))

 success =  ngood/(nbad+ngood)
 maskname = masks

; write outputs
openw,10,'zcomp.dat'
for i = 0, n_elements(masks)-1 do $
 printf,10,maskname(i),success(i),bluesuccess(i), redsuccess(i)
close, 10


return
end
