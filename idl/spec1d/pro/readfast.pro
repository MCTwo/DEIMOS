;+
; NAME:
;   READFAST
;
; PURPOSE:
;   to quickly read an ASCII table in regular columns. Data all
;   assumed float
;
; CALLING SEQUENCE:
;   readfast, filename, data, [header,] [ncols=ncols, skipline=skipline]
;
; INPUTS:
;   filename -- name of ascii file
;
; KEYWORD PARAMETERS:
;  ncols -- number columns to read
;  skipline -- number of lines at beginning of file to sky
;
; OUTPUTS:
;  data -- data array (nlines,ncols)
;  header -- lines skipped
;
; MODIFICATION HISTORY:
; jm00may24ucb
;-


; 
; routine to read in column-formated data very quickly.  at this point
; the routine is not very intelligent.  everything is assumed floating point.

pro readfast, filename, data, header, skipline=skipline, ncols=ncols

;	on_error, 2 ; return to user

;	t1 = systime(1)

        nrows=0L
        spawn, ['wc '+filename], stringout ; file size
        reads, stringout, nrows
        
        print, nrows, ' lines in file'
        openr, lun1, filename, /get_lun

        if keyword_set(skipline) then begin

            nrows = nrows-skipline 
            header = strarr(skipline)
            readf, lun1, header	; read the header

        endif

        data = fltarr(ncols,nrows)

        readf, lun1, data       ; read the data
        free_lun, lun1

;	print, systime(1)-t1

return
end

