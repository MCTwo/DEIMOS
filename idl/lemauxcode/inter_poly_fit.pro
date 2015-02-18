Function inter_poly_fit,x,y,degree,chisq=chisq,yfit=yfit,yerror=yerror,$
                        sig=sig,measure_errors=measure_errors,include=include

;+---------------------------------------------------------------------
;
; INTER_POLY_FIT	7/2002 GDB
;
; An interactive version of the instrinsic IDL routine POLY_FIT.
; Calling sequence, inputs, outputs and keywords are the same in so
; far as they have been included.
;
; The data points and the fit are plotted to the current device.  The
; user may select points to leave out from the fit and change the
; order of the polynomial.  If no degree is specified at the start, a
; defualt of 1 is used.
;
; Note: If degree is a named variable then it will be returned as the
;       final choice of the degree of fit made by the user.
;
; INPUTS/KEYWORDS/OTPUTS - See poly_fit
;
;       include   - Array with the same length as x and y denoting
;                   points to include in fit.  Points where 
;                   include = 1 will be included, while points where
;                   include = 0 will be excluded.  On output, include
;                   will contain values as set by the use.
;
;       Keywords 'chisq' and 'yerror' will return values that apply
;       only to points included in the fit.
;
; HISTORY
;	Written 7/16/2002 GDB
;       Revised 8/23/2004 GDB
;       INCLUDE keyword added 2/9/2005 GDB
;----------------------------------------------------------------------

if (n_params() lt 2) then begin
   print,'CALLING SEQUENCE: result = inter_poly_fit(x,y,degree,chisq=chisq,'
   print,'                             yfit=yfit,yerror=yerror,sig=sig,'
   print,'                             measure_errors=measure_errors'
   print,'                             include=include'
   return,-1
endif

;;; If no degree specified, start with default = 1.

if (n_params() eq 2) then degree=1

;;; Plot fit in order of increasing x
xorder = sort(x)

;;; Start by fitting all points, unlesss the include keyword
;;; is set.

if not(keyword_set(include)) then include = 1 + 0*x

done=0

while not(done) do begin

   ;;; 1. Compute the fit to the current subset of points.

   use  = where(include eq 1)
   xsub = x(use)
   ysub = y(use)

   if keyword_set(measure_errors) then begin
      me_use = measure_errors(use)
      pfit = poly_fit(xsub,ysub,degree,/double,chisq=chi2,yerror=yerr,$
                      sig=sigma,measure_errors=me_use)
   endif else begin
      pfit = poly_fit(xsub,ysub,degree,/double,chisq=chi2,yerror=yerr,$
                      sig=sigma)
   endelse

   ;;;; 2. Compute the fitted y values for all x values.
   
   yfit    = y
   if (degree eq 0) then yfit(*)=pfit else yfit(*)=pfit(0)
   
   if (degree gt 0) then begin
      for i=1,degree do yfit = yfit + (x^i)*pfit(i)
   endif
   
   ;;; Compute the chi-squared value for all points
   chi2tot = total((y-yfit)^2)
   
   ;;; 3. Plot the input data points and the new fit.
   ;;;    Put a box around the unused points.

   xlo    = min(x)
   xhi    = max(x)
   xrange = xhi-xlo
   xedge  = 0.05*xrange

   ylo    = min([min(y(use)),min(yfit)])
   yhi    = max([max(y(use)),max(yfit)])
   yrange = yhi-ylo
   yedge  = 0.2*yrange

   degprint     = strtrim(string(degree),2)
   chi2print    = strtrim(string(chi2),2) 
   chi2totprint = strtrim(string(chi2tot),2) 

   plot,x,y,psym=2,xrange=[xlo-xedge,xhi+xedge],yrange=[ylo-yedge,yhi+yedge],$
      xstyle=1,ystyle=1,xtitle='CHANGE DEGREE',ytitle='DONE',charsize=1.3,$
      title='Degree = '+degprint+'  ChiSq = '+chi2totprint+'  ('+chi2print+')'

   if keyword_set(measure_errors) then doerrors,x,y,measure_errors

   oplot,x(xorder),yfit(xorder)

   notused = where(include eq 0)
   if (total(notused) ne -1) then oplot,x(notused),y(notused),psym=6,symsize=1.5

   ;;; 4. Prompt user to make changes, if desired.

   print,'Click on a point to add or remove it from fit.'
   print,'Otherwise, click "CHANGE DEGREE" or "DONE"'

   cursor,cx,cy,/down,/data

   if ( (cy ge ylo-yedge) and (cy le yhi+yedge) and $
        (cx ge xlo-xedge) and (cx le xhi+xedge) ) then begin
      ;;; Determine to which point cursor is nearest
      ls = where( abs(x-cx) eq min(abs(x-cx)) )
      match = ls(0)
      ;;; Toggle the matching point
      if (include(match) eq 1) then include(match)=0  else  include(match)=1
      ;;; Do not allow only one point remaining!
      use = where(include eq 1,n_use)
      if (n_use lt 2) then begin
         print,'**Must keep at least two points!**'
         include(match)=1
      endif
   endif      

   if (cy lt ylo-yedge) then begin
      read,degree,prompt='Enter degree for fit: '
   endif

   if (cx le xlo-xedge) then done=1  

endwhile

chisq  = chi2tot
yerror = yerr
sig    = sigma

return,pfit
end
