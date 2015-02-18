;+
; NAME:
;          COUNTER
; PURPOSE:
;          Print a progress status to the screen on a single line.
;
; CALLING SEQUENCE:
;          COUNTER, NUMBER, OUTOF  [,INFOSTRING, /PERCENT, WAIT_TIME=variable]
;
; INPUTS:
;          NUMBER:  The current number. Usually the loop index variable
;          OUTOF:   The total number of iterations. The last loop index
;
; OPTIONAL INPUTS:
;
;          INFOSTRING: A string telling the user what is being
;                      counted e.g. 'Flat '
;
; KEYWORD PARAMETERS:
;         
;          PERCENT: Set to output update in percent completed 
;                   percent = rount(number/outof) * 100
;
;          WAIT_TIME:  Used for test and demo purposes only. See
;                      example below.
;
; OUTPUTS:
;          Status is printed to the screen and updated on a single line.
;
; SIDE EFFECTS:
;         This program takes much longer than a simple 
;         print statement. So use COUNTER judiciously. 
;         If your loop consists of only a couple 
;	  of relatively quick commands, updating the 
;	  status with this program could take up a 
;	  significant portion of the loop time!

; PROCEDURE:
;          Put counter statement inside your loop, preferably at the end.
;
; PROCEDURES CALLED:
;            
; EXAMPLE:
;          Try this to see how it works:
;
;          IDL> for i = 0,4 do counter,i,4,'test ',wait=.5
;
;
; MODIFICATION HISTORY:
;      Written by JohnJohn, Berkeley 06 January 2003
;-

pro counter,num,outof,infostring $
            ,wait_time = waittime $
            ,percent=percent      $
            ,clear=clear
on_error,2
clearline = fifteenb()          ;get "15b character to create a fresh line
if keyword_set(clear) then begin
    len = strlen(infostring)
    print,clearline,format='('+str(len)+'x, a)'
endif else begin
    if n_elements(infostring) eq 0 then infostring = 'Number ' 

    lenst = strtrim(strlen(infostring),2)
    leni = strtrim(strlen(strtrim(num,2)),2)
    leno = strtrim(strlen(strtrim(outof,2)),2)
    if keyword_set(percent) then begin
        per = strtrim(round(float(num)*100./outof),2)
        lenp = strtrim(strlen(strtrim(per,2)),2)
        form="($,a"+lenp+",' % Completed',a,a)"
        print,form = form, per, '         ', clearline
    endif else begin
        form="($,a"+lenst+",i"+leni+",' of ',i"+leno+",a,a)"
        print,form = form, infostring, num, outof, '         ',clearline
    endelse
endelse 
if n_elements(waittime) gt 0 then wait,waittime
end
