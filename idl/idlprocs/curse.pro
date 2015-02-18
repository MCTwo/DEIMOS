pro curse,image,diff=diff,delta=delta

;+
; NAME:
;          CURSE
;
; PURPOSE:
;          Print the position and optionally the value of a
;          mouse click on a plot. Wrapper for IDL's CURSOR 
;          procedure.
;
; CATEGORY:
;
; CALLING SEQUENCE:
;
;          CURSE [,IMAGE]
;
; OPTIONAL INPUTS:
;
;          IMAGE - 2-D array currently displayed. Z value
;                  is printed along with position.
;
; OUTPUTS:
;
;          Information on your screen
;
; PROCEDURE:
;          Plot or display data. Execute curse and click
;          on the image. When done, click negative X value.
;
; RESTRICTIONS:
;          This is best used for CCD or spectral data where
;          there are no negative ordinate or X-axis values.
;          
;          Don't worry about the name. It's just short for
;          CURSOR but way cooler sounding. Author not responsible
;          for actual curses.
;
; MODIFICATION HISTORY:
;          Written sometime in Dec 2002 by JohnJohn
; 02.07.2003  JJ  Modified so any click outside of the
;                 plot range causes CURSE to exit
;-

dum = min(!y.crange, ylo)
if not keyword_set(diff) then begin
    repeat begin
        cursor,x,y,/data,/up
        cool = x gt !x.crange[0] and y gt !y.crange[ylo]
        if cool then begin
            print,'x = ',strcompress(x,/rem)
            print,'y = ',strcompress(y,/rem)
            if n_elements(image) gt 0 then $
              print,'z = ',strcompress(image[x,y],/rem)
            print
        endif
    endrep until not cool
endif else begin
    x = fltarr(2)
    y = x
    delta = x
    for i = 0,1 do begin
        cursor,x1,y1,/data,/up
        x[i] = x1
        y[i] = y1
    endfor
    delta[0] = x[1] - x[0]
    delta[1] = y[1] - y[0]
    print
    print,'Delta x = ',strcompress(delta[0],/rem)
    print,'Delta y = ',strcompress(delta[1],/rem)
    r = sqrt(delta[0]^2 + delta[1]^2)
    print,'Delta r = ',strcompress(r,/rem)
    print
endelse

end
