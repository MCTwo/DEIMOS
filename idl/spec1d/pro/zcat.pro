;+
; NAME:
;   zcat
;
; PURPOSE:
;  Creates a new zcat out of the zspec .sav files.  Checks z's against
;  latest zresult files for discrepencies
;
; CALLING SEQUENCE:
;   zcat
;
; INPUTS:
;   none
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   zcat.yymmdd.fits       - Output structure zcat into zcat_archive dir
;             
; COMMENTS:
;  You will need to change the pointer zcat.latest.fits to point to
;  this new zcat file!
;
;
; BUGS:
;
; DATA FILES:
;   $DEEP2PRODUCTS/zspec_archive/zspec*.sav
;
; REVISION HISTORY:
;   ??  Written by DSM
;   09-Sep-2003  CC adds zresult files check
;------------------------------------------------------------------------------

function zcat_struct

   result = create_struct( $
    name = 'ZCAT', $
    'objno'      , 0L,  $
    'RA',         0.d, $
    'dec',       0.d, $
    'magB',      0., $
    'magR',      0., $
    'magI',      0., $
    'Pgal',      0., $
    'photoZ',    0., $
    'probgt07',  0., $
    'SFD_EBV',   0.0, $
; above from pcat, below from zresult
    'class'      ,  ' ', $
    'subclass'   ,  ' ', $
    'objname'    ,  ' ', $
    'maskname'   ,  ' ', $
    'slitname'   ,  ' ', $
    'date'       ,  ' ', $
    'mjd'        , 0.0, $
    'z'          , 0.0, $
    'zhelio'     , 0.0, $
    'z_err'      , 0.0, $
    'rchi2'      , 0.0, $
    'dof'        ,  0L, $
    'rchi2diff'  , 0.0, $
    'tfile'      ,  ' ', $
    'tcolumn'    , lonarr(10) - 1L, $
    'npoly'      ,  0L, $
    'theta'      , fltarr(10), $
    'vdisp'      , 0.0, $
    'vdisp_err'  , 0.0, $
    'zquality'   , 0 , $
    'comment'    , ' ', $
    'pca',        fltarr(5) -1000.$
   )

   return, result
end

;------------------------------------------------------

pro zcat,filename
;
;  construct zcat from the .sav files
;----------------

  ctbad=0
  filelist =  findfile('$DEEP2PRODUCTS/zspec_archive/zspec.*.sav', count=ct)
  print,filelist
  print,  ct,  ' files found.'
  if ct eq 0 then return

  topdir = getenv('D2_RESULTS')
  IF(topdir eq '') THEN print, 'Error: You need to set topdir.'

  for i=0, ct-1 do begin   
    restore,  filelist[i]
    zres = results.zresult[*]

    ;comments got stuck in results.comment way back in zspec
    zres.comment = results.comment

    ;create nlsky array to know if a z has been chosen with nlsky spec
    nlsky = (n_elements(nlsky) eq 0)? results.nonlocal_select : [nlsky,results.nonlocal_select]

;    ; Some entries are missing there dates... fix this
;    cnt = where(strcompress(zres.date, /rem) ne '')
;    date = zres[cnt[0]].date
;    cnt = where(strcompress(zres.date, /rem) eq '', c)
;    if c gt 0 then zres[cnt].date = date

    
;    zfile = findfile(topdir+'/zresult/zresult.' + zres[0].maskname $
;                     + '*.fits', count=ctt)
;  IF ctt lt 1 THEN begin 
;     print, 'Error: zcat file not found!'
;     return
;  endif else if ctt gt 1 then begin
;    print,  'More than one zcat found'
;    print,  zfile
;    read,  prompt='Enter index [0,1,..]:', indx
;  endif else if ctt eq 1 then $
;    indx= 0

;    hdr=headfits(zfile[indx],/silent)
;    mjd = sxpar(hdr,'MJD-OBS')
;    zres.date = mjd

    zz = replicate(zcat_struct(),n_elements(zres))
    struct_assign,zres,zz
    zres = zz
 
    if i eq 0 then zbest = zres else $
      zbest= [zres, zbest]
  endfor

  zzz = replicate(zcat_struct(), n_elements(zbest))
  struct_assign, zbest,  zzz ; copy structure

  zcat = zzz[sort(zzz.objname)] ;sort all by object name

  ;remove serendips, how did they get in here??
  objname = strcompress(zcat.objname, /rem)
  nonser = where(strmid(objname,0,1) ne 's')
  zcat = zcat[nonser] 

  ;set objs which used nlsky spectra to find z to q=0
  wh = where(nlsky eq 1,ct)
  zcat[wh].zquality = 0

  objects = zcat.objname
  pcats = strmid(objects, 0, 2) ;which pcat does this derive from?

  ;pcats[where(pcats gt '55')] = '-1' ;set bogus values to -1

  pcatfix = fix(pcats) ; some pcats ='s1' what does this mean?
  upcat = pcatfix[uniq(pcatfix)] ;get unique list

 
  print, 'drawing from following pcat files: ', upcat
  for i=0,  n_elements(upcat)-1 do begin
    if upcat[i] gt 0 then begin ;only deal with valid pcat files
      photo = getphot(upcat[i]) ;input photometry
      select = where(pcatfix EQ upcat[i], nselect) ;objects in this pcat file
      onumber = long(zcat[select].objname) ;get position in list
; N.B. this does NOT properly deal with faint objects ADDED to pcat
; list when merged objects were split
      print, nselect, ' objects found in field ', upcat[i]
      pointer = lonarr(nselect)
      for j=0,  nselect-1 do pointer[j] = (where(photo.objno EQ onumber[j]))[0]
      ;take only first entry if multiple matches

      psel = photo[pointer] ;selected objects from pcat
      ;copy over the relevant information
      zcat[select].RA = psel.RA
      zcat[select].dec = psel.dec
      zcat[select].objno = psel.objno
      zcat[select].magB = psel.magB
      zcat[select].magR = psel.magR
      zcat[select].magI = psel.magI
      zcat[select].Pgal = psel.Pgal
      zcat[select].photoZ = psel.photoZ
      zcat[select].probgt07 = psel.probgt07
      zcat[select].SFD_EBV = psel.SFD_EBV
    endif
  endfor


  for i=0L, n_elements(zcat)-1 do begin

    ;break up date into segments (yyyy mm dd)
    julian_day = zcat[i].mjd + double(2400000.5)
    vcorr = heliocentric(zcat[i].ra, zcat[i].dec, jd= julian_day, $
      longitude=360.-155.47, latitude= 19.826,  altitude=4205.)
    ;coordinates of Keck-II.
    zcat[i].zhelio = zcat[i].z +vcorr/299792. ;correct to heliocentric velocity

  endfor

  ;this is crude, fix later
  wh=where(zcat.zquality eq 255,ct)
  if ct gt 0 then zcat[wh].zquality = -1

  ;-----------
  ;match up z to latest zresults files, set q=2 if values are different
  print,'Updating zspec files to zresult data...'
 
  masklist = zcat[sort(zcat.maskname)].maskname
  masklist = masklist[uniq(masklist)]
 ; if n_elements(uniq(masklist)) ne n_elements(filelist) then begin
 ;     print, 'repeat zspec.sav files, fix and restart' 
 ;     return
 ; endif
      
  for m=0,n_elements(masklist)-1 do begin

      objs = where(zcat.maskname eq masklist[m])
      zres_file = findfile('$D2_RESULTS/zresult/zresult*'+strtrim(masklist[m],2)+'*')

      if n_elements(zres_file) gt 1 then begin
          ;pick most recent observation of mask
          if masklist[m] eq 4203 then zres_file=zres_file[0]
          zres_file = zres_file[n_elements(zres_file)-1]
      endif 

      zresults = mrdfits(zres_file[0],2,/silent) ;2nd HDU contains all 10 fits

      for r=0,n_elements(objs)-1 do begin
          k=objs[r]        
          ; we only care about q=3&4
          if zcat[k].zquality ne 3 and zcat[k].zquality ne 4 then continue
       
          if zcat[k].z lt .001 then begin ;this hangs up the code, its probably a star
              print,'star, skipping'
              continue
          endif

          wh = where(zresults.objname eq zcat[k].objname,ct)

          if ct eq 0 then begin
              print,'obj not in zres! setting q=2...'+strtrim(zcat[k].maskname,2)
              zcat[k].zquality = 2
              continue
          endif

          if ct ne 10 then $
              print, 'couldnt find all 10 fits...or more than 10 fits'+strtrim(ct,2)    

          zall = zresults[wh].z
          diff = abs(zall - zcat[k].z)/zcat[k].z
          mind = min(diff,cmin)
          goodfit = where(diff eq mind) ; now we know where the good fit came from
          IF(mind GT 0.01) THEN BEGIN 
              print, 'Miss-match >1%.  Object set to Q=2, diff='+$
                strtrim(mind,2)+' '+strtrim(zcat[k].maskname,2)
              zcat[k].zquality = 2         
              ctbad = ctbad+1
              continue
          ENDIF                                


          if n_elements(goodfit) gt 1 then begin
              gwh = where(zresults[goodfit].z_err eq min(zresults[goodfit].z_err))
              goodfit = goodfit[gwh]
              goodfit = goodfit[0]
          endif

          ;copy new zresult entries into the zcat
          tmp = zresults[wh]
          zcat[k].class = tmp[goodfit].class
          zcat[k].subclass = tmp[goodfit].subclass
          zcat[k].objname = tmp[goodfit].objname
          zcat[k].slitname = tmp[goodfit].slitname
;          zcat[k].maskname = tmp[goodfit].maskname
          zcat[k].date = tmp[goodfit].date
          zcat[k].mjd = tmp[goodfit].mjd
          zcat[k].z = tmp[goodfit].z
          zcat[k].z_err = tmp[goodfit].z_err
          zcat[k].rchi2 = tmp[goodfit].rchi2
          zcat[k].dof = tmp[goodfit].dof
          zcat[k].rchi2diff = tmp[goodfit].rchi2diff
          zcat[k].tfile = tmp[goodfit].tfile
          zcat[k].tcolumn = tmp[goodfit].tcolumn
          zcat[k].npoly = tmp[goodfit].npoly
          zcat[k].theta = tmp[goodfit].theta
          zcat[k].vdisp = tmp[goodfit].vdisp
          zcat[k].vdisp_err = tmp[goodfit].vdisp_err
      ;    zcat[k].zquality = tmp[goodfit].zquality
      ;    zcat[k].comment = tmp[goodfit].comment 

     endfor  
  endfor 

  print,'total # of objects set to q=2 == '+strtrim(ctbad,2)

  wh = where(nlsky eq 1,ct)
  print,strmid(ct,2)+' Objects set to Q=0 b/c nlsky was used for z'

  productsdir = getenv('DEEP2PRODUCTS')
  ;get current date in yymmdd
  date = systime()
  if strmid(date,4,3) eq 'Jan' then month = '01'
  if strmid(date,4,3) eq 'Feb' then month = '02'
  if strmid(date,4,3) eq 'Mar' then month = '03'
  if strmid(date,4,3) eq 'Apr' then month = '04'
  if strmid(date,4,3) eq 'May' then month = '05'
  if strmid(date,4,3) eq 'Jun' then month = '06'
  if strmid(date,4,3) eq 'Jul' then month = '07'
  if strmid(date,4,3) eq 'Aug' then month = '08'
  if strmid(date,4,3) eq 'Sep' then month = '09'
  if strmid(date,4,3) eq 'Oct' then month = '10'
  if strmid(date,4,3) eq 'Nov' then month = '11'
  if strmid(date,4,3) eq 'Dec' then month = '12'
  if strmid(date,8,1) eq ' ' then day = '0'+strmid(date,9,1) $
  else day = strmid(date,8,2)
  filename = 'zcat.'+strmid(date,22,2)+month+day+'.fits'
  filename = productsdir+'/zcat_archive/'+filename

  print,'Creating new zcat....'
  print,filename

  ffile = findfile(filename,count=ct)
  if ct gt 0 then filename = filename+'.tmp'
  mwrfits, zcat, filename, /create

  
  ;update zcat.latest.fits
  spawn,'/bin/rm $DEEP2PRODUCTS/zcat.latest.fits'
  spawn,'ln -s '+filename+$
    ' $DEEP2PRODUCTS/zcat.latest.fits'

  ;run zcomp and replace zcomp.dat
  print,'moving old zcomp.dat to zspec_archive/junkish/zcomp.dat.old'
  spawn,$
    '/bin/mv $DEEP2PRODUCTS/zspec_archive/zcomp.dat $DEEP2PRODUCTS/zspec_archive/junkish/zcomp.dat.old'
  zcomp
  spawn,'/bin/mv zcomp.dat $DEEP2PRODUCTS/zspec_archive/zcomp.dat'
  print,'zcomp.dat created and stored in $DEEP2PRODUCTS/zspec_archive/zcomp.dat'

END
