;+
; NAME:
;   QUICK_SCIENCE
; PURPOSE:
;   check for calibration files and start quicklook reductions
; CALLING SEQUENCE:
;   quick_science
; INPUTS:
;
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; RESTRICTIONS:
;   run in data directory
; EXAMPLES:
;
; COMMENTS:
;
; REVISION HISTORY:
;
;----------------------------------------------------------------------




pro quick_science, level,  newfile
   
err=''

  if n_elements(level) eq 0 then level = 2
  if n_elements(newfile) eq 0 then begin
     spawn, 'rsh polo wfi', newfile
     newfile = newfile[0]
     h = headfits('/net/polo'+newfile,err=err)
  endif else begin
     newfile = newfile[0]
     h = headfits('/net/polo'+newfile,err=err)
  endelse 

  if strlen(newfile) lt 5 or strlen(err) gt 0 then begin
      print,'Warning: Filename from wfi was blank or bad'
      mask='  '
      isscience=0
      wrongmode=0
      
  endif else begin
    
      mask = sxpar(h, 'SLMSKNAM')
      lamps = sxpar(h, 'LAMPS')
      exptime = sxpar(h, 'EXPTIME')
      obstype = sxpar(h, 'OBSTYPE')
      grating = sxpar(h, 'GRATENAM')
      object = sxpar(h, 'OBJECT')
      gratepos = sxpar(h, 'GRATEPOS')
      hatchpos = sxpar(h, 'HATCHPOS')
      if gratepos eq 3 then wave = sxpar(h, 'G3TLTWAV') $
         else wave = sxpar(h, 'G4TLTWAV') 
      mosmode=sxpar(h,'MOSMODE')


      isscience =  strpos(lamps, 'Off') ge 0 $ 
        AND exptime gt 250. $
        AND strpos(obstype, 'Object') ge 0 $
        AND wave gt 5000 $
        AND strpos(hatchpos, 'open') ge 0 $
        AND mosmode eq 'Spectral'

; to trap errors like 18 June 2004

      wrongmode =  strpos(lamps, 'Off') ge 0 $ 
        AND exptime gt 250. $
        AND strpos(obstype, 'Object') ge 0 $
        AND wave gt 5000 $
        AND strpos(hatchpos, 'open') ge 0 $
        AND mosmode NE 'Spectral'
  endelse


if wrongmode then begin

       spawn,'whoami',account
 
       spawn, '/home/kics/instr/bin/play_sound -v 99 -host pohue -account '+account+' /home/deepteam/sounds/doh.au &'

       spawn, '/home/kics/instr/bin/play_sound -v 99 -host hamoa -account '+account+' /home/deepteam/sounds/doh.au &'

    openw, 2, '/home/'+account+'/temp/error.txt'
              printf, 2
              printf, 2, 'Warning: Check if you are in spectral readout mode!'
    close,2

     spawnstring = 'cat /home/'+account+'/temp/error.txt | tkmessage -type warning &' 



     spawn, spawnstring
 endif

; and now back to the usual routine

  mask = strcompress(mask, /REMOVE)
; deal with SN or ERO masks


;  if strlen(mask) gt 6 then mask = strmid(mask, 0, 6)

     maskprocessed = n_elements(findfile(mask+'/arc*')) ge 4 $
       AND n_elements(findfile(mask+'/cal*')) ge 10


     stringsep = strsplit(strcompress(newfile, /REMOVE), '/', /extract)
     filestem = stringsep[n_elements(stringsep)-1]

     directory = '/net/polo'
     for i=0, n_elements(stringsep)-2 do directory = directory+'/'+stringsep[i]
     cd,directory
     cd,current=cwd

;     print, 'in directory '+directory
     print,  'processing file ',newfile
;       print, 'mask: ', mask
    

    if isscience AND maskprocessed then begin
       print, 'this is a science file!'
       
       cd, mask
; write quickSlit files
;            string = 'set clobber ; echo "deimos_2dreduce,
;            file='+filestem+', quick='+string(level < 9,
;            format='(i2)')+'" | idl >> science.log'
            
       deimos_2dreduce, file=filestem, quick=level
       if file_test('spSlit*.*') eq 1 then begin
          list=findfile('spSlit*.fits')
          spslit_combine,list, nlsky = 0
          slitfiles = findfile('slit*.fits', count=nfiles)
          
          do_extract, files=slitfiles, nsigma_optimal=1.75, nsigma_boxcar=1.1 

       endif
       quick_sciqa, quicklevel=level, file=filestem
       cd,cwd

    print, 'mask: ', mask, ' completed'


    endif


    quick_science, level

return
end


