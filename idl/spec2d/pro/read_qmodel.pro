pro read_qmodel,  qmodelpath, grating,slider,mu=mu,o3=o3,roll3=roll3

;+
; NAME:
;   read_qmodel
;
; PURPOSE:
;   Reads in an ASCII file detailing optical model parameters
;
; CALLING SEQUENCE:
;   read_qmodel, qmodelpath, grating,slider,mu=mu,o3=o3,roll3=roll3
;   
; INPUTS:
;   qmodelpath -- path to search
;   grating -- grating to check
;   slider  -- slider to check for
;
; OUTPUTS:
;   mu  -- mu value
;   o3  -- o3 value
;   roll3  -- roll3 value
;
; COMMENTS:
;
; REVISION HISTORY:
;   23may09 first version
;----------------------------------------------------------------------


  deimos_data = getenv('DEIMOS_DATA')+'/'
;  if deimos_data eq '/' then message, 'You need to set $DEIMOS_DATA!'
  if deimos_data eq '/' then deimos_data = '' ;not always set


  Foo = findfile(Strcompress(concat_dir(qmodelpath,$
   'qmodel.'+string(grating,format='(i4)')+'*slider*') $
    ,/remove), count=filect)

  if filect eq 0 then begin 
     print, 'Cannot find qmodel parameter file: ', qmodelfile
     retall
  endif 
  

  openr, unit, foo[0], /get_lun
  
  line='' 
  mu=-1
  o3=-1
  roll3=-1

  while NOT EOF(unit) do begin

    readf, unit, line

;    substring = strsplit(line, /extract) ;split into two or more 
    substring = strsplit(strcompress(line), $
                         '[^A-Za-z0-9._#/]+', /REGEX, /extract)

    case substring[0] of 
      'qmodel.mu': mu=float(substring[1])
      'qmodel.o3': o3=float(substring[1])
      'qmodel.roll3': roll3=float(substring[1])
       else: a = 0 ;junk statement
      
    endcase

  endwhile


  close, unit
  return
end 






