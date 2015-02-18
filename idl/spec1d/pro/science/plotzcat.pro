;+
; NAME:
;   plotzcat
;
; PURPOSE:
;   make a plot in redshift space
;
; CALLING SEQUENCE:
;   plotzcat,zcat
; 
; INPUTS:
;   zcat -- zcat structure 
;
; COMMENTS:
;   TBD: needs to use comoving coordinates r(z) 
;
; REVISION HISTORY:
;   md 5nov023
;----------------------------------------------------------------------
pro plotzcat, zcat

  name = strmid(zcat[0].objname, 0, 2)

;setup colors
  rtiny   = [0, 1, 0, 0, 0, 1, 1, 1]
  gtiny = [0, 0, 1, 0, 1, 0, 1, 1]
  btiny  = [0, 0, 0, 1, 1, 1, 0, 1]
  tvlct, 255*rtiny, 255*gtiny, 255*btiny

tvlct, [255],[255],[255], !d.table_size-1

  zr = [0.7, 0.9]
  plot, zcat.z, zcat.ra*zcat.z, /nodata, xr=zr, yr=[-150., 100.], $
    xtit='Redshift z -0.7 -i*.2', ytit='Transverse Mpc', charsize=2., $
    charthick=1.5, title='Field ' + name
 
  coh = 3.e5/100. ;c/H_0
  meanra = mean(zcat.ra)
  meandec = mean(zcat.dec)
  transverse = (zcat.ra - meanra)*cosd(meandec)/!radeg*coh*zcat.z 
;rough transverse distance
  maxt = .33/!radeg*coh*1.5 ;maximum cross dimension

  for i=0, 3 do begin
    zstart = 0.7 + i*.2
    zend = zstart+0.2
    yr = [zstart, zend]/1.5*maxt ;transverse distance limits
    jj = where(zcat.z gt zstart and zcat.z le zend)
    zz = zcat[jj]
    tt = transverse[jj]
    jg = where(zz.dec gt meandec,  comp=jb)
    oplot, zz[jg].z -zstart+0.7, 75.-60.*i + tt[jg], $
      color=3, psym=1, syms=.5 ;red + if above mean dec
    oplot, zz[jb].z -zstart+0.7, 75.-60.*i + tt[jb], $
       psym=4, syms=.5 ; white if below mean dec
    plots, zr,   75.-60.*i + yr, line=1
    plots, zr  , 75.-60.*i - yr, line=1
    xyouts, .71, 90.-i*57., $
       'Zrange: '+string(zstart, format='(f4.2)')+':'+ string(zend, $
        format='(f4.2)'), charsize=2., charthick=1.5
  endfor
    













end
