pro fix_nspec_quadnoise, filenames

;******************************************************************
;
; NAME: fix_nspec_quadnoise
;
; PURPOSE:  when the 8 line pattern noise shows up in NIRSPEC data,
;		this routine will attempt to remove it
;
; CALLING SEQUENCE: fix_nspec_quadnoise, 'filenames'
;
; INPUTS: filenames is a single column text file
;
; OPTIONAL INPUT PARAMETERS:
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;
; OUTPUTS: 1 fits file per input...with a "c" prefix added to the name
;
; EXAMPLE: fix_nspec_quadnoise, 'filenames'
;
; NOTES:
;
; RESTRICTIONS:
;
; MODIFICATION HISTORY:
;	2007 Mar 06 Jim Lyke W.M. Keck Observatory	Original
;	2008 Apr 22 Alice Shapley UCLA changed algorithm for case of row 6/8 being the bad one
;
;******************************************************************


; pattern repeats every 8 lines:
rep = 8

readcol, filenames, filelist, format="A"
for k=0, n_elements(filelist)-1 do begin
  a=readfits(filelist[k], h)
  xsz = (size(a))[1]
  ysz = (size(a))[2]

  ; here are the definitions to break an image into quadrants
  ; generally, only one quadrant is bad...
  quad1 = a(xsz/2.:xsz-1, ysz/2.:ysz-1)
  quad2 = a(xsz/2.:xsz-1, 0:ysz/2.-1)
  quad3 = a(0:xsz/2.-1, 0:ysz/2.-1)
  quad4 = a(0:xsz/2.-1, ysz/2.:ysz-1)

  ; all quadrants should be the same size
  quadx = (size(quad1))[1]
  quady = (size(quad1))[2]

  num_bad_lines = quady / rep

  ; correct quadrant I.  If you change quadrant, you also have to
  ; change the statement that copies the fixed quadrant back into
  ; the image
  b = quad1

  for i=0, (quady)/rep-1 do begin
     b1=median(b[400:499,i*rep+0])
     b2=median(b[400:499,i*rep+1])
     b3=median(b[400:499,i*rep+2])
     b4=median(b[400:499,i*rep+3])
     b5=median(b[400:499,i*rep+4])
     bbad=median(b[400:499,i*rep+5])
     b7=median(b[400:499,i*rep+6])
     b8=median(b[400:499,i*rep+7])
     ;print, b1, b2, b3, b4, b5, bbad, b7, b8	
     bgoods=[b1, b2, b3, b4, b5, b7, b8]
     bgood=median(bgoods)
     bdiff=bbad - bgood
     b[0:quadx-1,i*rep+5]=b[0:quadx-1,i*rep+5]-bdiff
  endfor

  ; copy back into original image
  ; for quad1 (use definitions at the top...)
  a(xsz/2.:xsz-1, ysz/2.:ysz-1) = b

  dotpos=strpos(filelist[k], '.', /reverse_search)
  rootfn=strmid(filelist[k],  0, dotpos)  
  ;writefits, rootfn+'c.fits', a, h
  writefits, 'c'+rootfn+'.fits', a, h

endfor


end
