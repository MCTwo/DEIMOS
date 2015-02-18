;+
; NAME:
;   zslice
;
; PURPOSE:
;   make a plot projected on the sky for a slice of redshift
;
; CALLING SEQUENCE:
;   zslice,zcat,zrange
; 
; INPUTS:
;   zcat -- input structure
;   zrange -- [zlow,zhigh] no cut made if zrange unspecified
;
; COMMENTS:
;
; REVISION HISTORY:
;   md 06nov02
;----------------------------------------------------------------------
pro zslice, zcat, zrange

  name = strmid(zcat[0].objname, 0, 2)

;setup colors
  rtiny   = [0, 1, 0, 0, 0, 1, 1, 1]
  gtiny = [0, 0, 1, 0, 1, 0, 1, 1]
  btiny  = [0, 0, 0, 1, 1, 1, 0, 1]
  tvlct, 255*rtiny, 255*gtiny, 255*btiny

  tvlct, [255],[255],[255], !d.table_size-1
  if n_params() lt 2 then zrange = [0., 100.] ;full range if not specified

  xr = [max(zcat.ra), min(zcat.ra)]
  plot, zcat.ra, zcat.dec, /nodata, xr=xr , $
    xtit='Right ascension', ytit='declination', charsize=2., $
    charthick=1.5, title='Field ' + name + '  Z-limits: [' + $
    string(zrange[0], format='(f4.2)') +','+  $
    string(zrange[1], format='(f4.2)')+']'
 

  j = where(zcat.z gt zrange[0] and zcat.z le zrange[1], nsel)

  if nsel gt 0 then $
     oplot, zcat[j].ra, zcat[j].dec, psym=4, syms=0.5

end
