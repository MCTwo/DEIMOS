;+
; NAME:
;   get_airmass
;
; PURPOSE:
;   returns the average airmass that an object was observed at by
;       reading the plan file for that reduction and using the header values
;
; CALLING SEQUENCE:
;   airmass= get_airmass(planfile)
; 
; INPUTS:
;  planfile - string containing name of plan file. If not supplied,
;             then plan file in current directory is used
;                   
;
; OUTPUTS:
;   airmass - average airmass object was observed at 
;
;
; MODIFICATION HISTORY:
;   07oct02 A.Coil
;   07oct02 MD
;-

function get_airmass, planfile

   if n_params() eq 0 then planfile = '*plan' ;give a name!
   read_planfile, planfile, maskno, rawdatadir, outdatadir, flatnames, $
           arcnames, sciencenames

   nframes = n_elements(sciencenames)
   airmass_array = fltarr(nframes)

   for i = 0, nframes-1 do begin
      hdr = headfits_1(rawdatadir+sciencenames[i])
      airmass_array[i] = sxpar(hdr, 'AIRMASS')
   endfor

return,  mean(float(airmass_array))
end



