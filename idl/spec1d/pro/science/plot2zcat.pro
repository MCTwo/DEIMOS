;+
; NAME:
;   plot2zcat
;
; PURPOSE:
;   make plots in redshift space that differentiate symbols for two
;   different class of objects
;
; CALLING SEQUENCE:
;   plotzcat,zcat1,zcat2,[title=title]
; 
; INPUTS:
;   zcat1 -- zcat structure - plotted in blue
;   zcat2 -- zcat structure - plotted in red
;
; KEYWORDS:
;   title -- if used, add to plottitle
;
; COMMENTS:
;   TBD: needs to use comoving coordinates r(z) 
;
; REVISION HISTORY:
;   md 5nov023
;----------------------------------------------------------------------
pro plot2zcat, zcat1, zcat2, title=title

  name = strmid(zcat1[0].objname, 0, 2)
  if keyword_set(title) then name = name + ' ' + title
  zcat = [zcat1,  zcat2] ;concatenate for boundary checking

;setup colors
  rtiny   = [0, 1, 0, 0, 0, 1, 1, 1]
  gtiny = [0, 0, 1, 0, 1, 0, 1, 1]
  btiny  = [0, 0, 0, 1, 1, 1, 0, 1]
;  tvlct, 255*rtiny, 255*gtiny, 255*btiny

;tvlct, [255],[255],[255], !d.table_size-1

  zr = [0.7, 0.9]
  plot, zcat.z, zcat.ra*zcat.z, /nodata, xr=zr, yr=[-150., 100.], $
    xtit='Redshift z -0.7 -i*.2', ytit='Transverse Mpc', charsize=2., $
    charthick=1.5, title='Field ' + name
 
  coh = 3.e5/100. ;c/H_0
  meanra = mean(zcat.ra)
  meandec = mean(zcat.dec)
  transverse1 = (zcat1.ra - meanra)*cosd(meandec)/!radeg*coh*zcat1.z 
  transverse2 = (zcat2.ra - meanra)*cosd(meandec)/!radeg*coh*zcat2.z 
;rough transverse distance
  maxt = .33/!radeg*coh*1.5 ;maximum cross dimensio
  for i=0, 3 do begin
    zstart = 0.7 + i*.2
    zend = zstart+0.2
    yr = [zstart, zend]/1.5*maxt ;transverse distance limits
    jb = where(zcat1.z gt zstart and zcat1.z le zend)
    zz = zcat1[jb]
    tt = transverse1[jb]
    oplot, zz.z -zstart+0.7, 75.-60.*i + tt, $
       psym=1, syms=.5 , color=3 ; blue for objects in zcat1


    jr = where(zcat2.z gt zstart and zcat2.z le zend)
    zz = zcat2[jr]
    tt = transverse2[jr]
    oplot, zz.z -zstart+0.7, 75.-60.*i + tt, $
       psym=4, syms=1. , color=1 ; red for objects in zcat2, bigger as well
    plots, zr,   75.-60.*i + yr, line=1
    plots, zr  , 75.-60.*i - yr, line=1
    xyouts, .71, 90.-i*57., $
       'Zrange: '+string(zstart, format='(f4.2)')+':'+ string(zend, $
        format='(f4.2)'), charsize=2., charthick=1.


  endfor
    













end
