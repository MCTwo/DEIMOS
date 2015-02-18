;+
; NAME:
;  add_zquality
;
; PURPOSE:
;  Adds the zquality flags from pspec2.pro to the zcat file constructed 
;  using makezcat.pro.  The output file is identical to that of the
;  original zcat.fits, except that it has now filled in the entries for
;  the redshift quality, comments and canned comments tags.  It has
;  also overwritten the z tag with the pspec2 specified redshift.
;
; CALLING SEQUENCE:
;  add_zquality, 'zcat.fits'
;
; OUTPUT:
;  The output file is overwritten to the original filename
;
; STILL TO BE DONE:
;  1) This procedure does not currently add the canned comments results.
;     Need to convert the aeor string to a byte first.
;  2) Repeats are discriminated only on the basis of mask name.  The 
;     date of observation will also need to be incorporated for 
;     multiply observed masks.
;
;
;  REVISION HISTORY:  
;   DSM 11-25-2002
;
;-------------------------------------------------------------

FUNCTION pspec_struct
; 
; Holds all the data that is output by pspec2
;*********************

   result = create_struct( $
    name = 'PSPEC', $
    'objno'      , 0L, $
    'maskno'     , 0L, $
    'slitno'     , 0L, $    
    'z'          , 0.0, $
    'zquality',  0, $
    'ccomment',  ' ', $
    'mcomment',  ' ' $
   )

  return, result
END







FUNCTION read_redshifts, filename
;
; Read in the redshifts etc from the pspec2 output files
;******************************************


  ; Variables in file
  ;----------------
  buffer = ''
  objn = 0L
  slitn = 0
  zz =  0.0
  q =  0
  can = ''
  com = ''

  ; Extract mask name from filename
  ;------------------------------
  dirlist = STRSPLIT(filename, '/', /extract)
  ndir = N_ELEMENTS(dirlist)
  maskname = STRMID(dirlist[ndir-1],0,4)
  PRINT, 'Reading in from mask ', maskname

  OPENR, unit, filename, /GET_LUN
  READF, unit, buffer
  READF, unit, buffer

  ; Go through file to determine the number of lines
  ;-------------------------------------------------
  nlines=0
  WHILE NOT EOF(unit) DO BEGIN
   READF, unit, buffer
   IF buffer EQ '' THEN CONTINUE  ; Needed if blank line at end of file
   READS, buffer, objn, slitn, z, q, can, com, FORMAT='(I9,I6,F9.5,I4,A7,A0)'
   nlines=nlines+1
  ENDWHILE
  FREE_LUN, unit


  ; Allocate a structure
  ;--------------------
  zz = replicate(pspec_struct(),nlines)



  OPENR, unit, filename, /GET_LUN
  READF, unit, buffer
  READF, unit, buffer

  ; This time read in data
  ;------------------------
  FOR i=0, nlines-1 DO BEGIN
   READF, unit, buffer
   READS, buffer, objn, slitn, z, q, can, com, FORMAT='(I9,I6,F9.5,I4,A7,A0)'
   zz[i].objno      = objn
   zz[i].slitno     = slitn
   zz[i].maskno     = FIX(maskname)
   zz[i].z          = z
   zz[i].zquality   = q
   zz[i].ccomment   = STRCOMPRESS(can, /REMOVE_ALL)
   zz[i].mcomment   = STRTRIM(com,2)
  ENDFOR
  FREE_LUN, unit

  zz.mcomment = STRTRIM(zz.mcomment,2)
 RETURN, zz
END






PRO add_zquality, catname


  ; Read in the current zcat file
  ;-------------------------------
  zcat = MRDFITS(catname,1)



  ; Find the completed redshift files (ascii)
  ;------------------------------------------
  filelist = FINDFILE('/g/dsm/Redshifts/Complete/*redshifts.v1_0*', count=nfile)
  PRINT, nfile, ' redshift files found.'

  FOR i=0, nfile-1 DO BEGIN
    pspec = read_redshifts(filelist[i])
    pspec_all = (i eq 0) ? pspec : [pspec_all,pspec]
  ENDFOR



  ; Now append the pspec results to the zcat structure
  ;------------------------------------------------
  PRINT, 'Now matching to ', catname, '...'

  obj_names = pspec_all.objno
  mask_names = pspec_all.maskno
  counter = 0L
  rep_counter = 0L
  FOR i=0L, N_ELEMENTS(obj_names)-1 DO BEGIN
    cnt = WHERE(zcat.objno EQ obj_names[i],c)
    IF c LT 1 THEN BEGIN
      ; This shouldn't happen

      PRINT, 'ERROR:', c, ' instances of ', obj_names[i], ' found.'
      PRINT, 'Skipping this object.'
      CONTINUE
    ENDIF ELSE IF c GT 1 THEN BEGIN
      ; If the object is entered in zcat more than once, then
      ; we match it according to the name of the mask it was
      ; observed in.... will not work for repeated masks!!!

      cnt = WHERE(LONG(zcat.maskname) EQ mask_names[i] AND $
                     zcat.objno EQ obj_names[i],c)
      rep_counter = rep_counter+1
      IF c NE 1 THEN BEGIN
        PRINT, 'MAJOR ERROR: Please check the handling of repeated objects!'
        PRINT, '   Taking first....'
      ENDIF
    ENDIF
    counter = counter+1L
    zcat[cnt[0]].z        = pspec_all[i].z
    zcat[cnt[0]].zquality = pspec_all[i].zquality

    ; The strtrim still doesn't work, which is annoying...
    comm = pspec_all[i].mcomment
    zcat[cnt[0]].comment  = STRTRIM(comm,2)
;    zcat[cnt[0]].ccomment = pspec_all[i].ccomment  

    ; Heliocentric corrections
    date = zcat[cnt[0]].date
    date = strsplit(date, '-', /extract) 
    julian_day = julday(date[1], date[2], date[0])
    vcorr = heliocentric(zcat[cnt[0]].ra, zcat[cnt[0]].dec, jd= julian_day, $
      longitude=360.-155.47, latitude= 19.826,  altitude=4205.) 
    zcat[cnt[0]].zhelio = pspec_all[i].z +vcorr/299792
  ENDFOR

  PRINT, 'Total of ', counter, ' redshifts entered.'
  PRINT, ' of which ', rep_counter, ' were observed more than once.'

  ; Write new fits file
  ;---------------------
  mwrfits, zcat, catname, /create

END
