Function rflam,filename,hdr,lambda,float=float,invcm=invcm,$
               vacuum=vacuum,_extra=ex

;+---------------------------------------------------------------------
;
; RFLAM		6/20002
;
; Uses READFITS to read in a fits file.  Additionally, uses the header
; keywords CTYPE1, CRVAL1, and CDELT1 to construct a wavelength array.
; Wavelength spacing may either be a constant (CTYPE1 = 'LAMBDA'),
; logarithmic (CTYPE1 = 'LINEAR'), or in 'MULTISPE' format (CTYPE1 =
; 'MULTISPE').,
;
; Note: There is a conflict between (at least) the way MAKEE programs
; and HST denote linear wavelength scaling, the former by setting
; CTYPE1 = 'LAMBDA' and the latter by setting CTYPE1 = 'LINEAR'.  In
; contrast, MAKEE programs use 'LINEAR' to denote log-linear scaling.
; In order to resolve this conflict, any spectrum which has 'CRVAL1' >
; 10 is assumed to have true linear wavelength scaling.  Otherwise the
; first pixel would have a wavelength > 10^10.
;
; INPUTS
;	filename	- Name of fits file
;
; KEYWORDS
;	float		- If set, lambda is returned as a 
;		 	  single precision floating
;			  point array.  Default is double
;			  precision.  Note: Since header
;			  keywords are often double
;			  precision, the calculations for
;			  lambda are done in double precision
;		   	  arithmetic regarless of the output
;			  precision.
;       invcm           - Convert wavelength units
;                         from inverse cm to angstroms and reverse
;                         array order to run from short to long
;                         wavelengths.
;	_extra		- Accepts READFITS keywords
;       vacuum          - Convert wavelengths from air to vaccum
;	
; OUTPUTS
;	<result>	- Floating point array of values read
;			  from fits file.
;	hdr		- Fits header (string)
;	lambda		- Wavelength array.  Default is double
;			  precision floating point, unless 'float'
;                         keyword is set..
;
; HISTORY
;	Written 6/19/2002 GDB
;       Added failsafe linear wavelength recognition 
;          11/19/2002 GDB
;	DOUBLE keyword added 2/19/2003 GDB
;       INVCM keyword added 6/21/2004 GDB
;       DOUBLE replaced with FLOAT keyword 5/5/2005 GDB
;       Revised to give lambda = -1 when not enough header 
;          information available  5/17/2005 GDB
;       Added MULTIPE format capability and VACUUM keyword 9/21/2006 GDB
;----------------------------------------------------------------------- 

if (n_params() eq 0) then begin
   print,'CALLING SEQUNCE: result = rflam(filename,hdr,lambda,/float,'
   print,'                                /vacuum,/invcm,_extra)'
   return,-1
endif

;;; Read in the data
data = readfits(filename,hdr,_extra=ex)

;;; Default is no wavelengths
lambda = -1

;;; Get wavelength solution type
type   = sxpar(hdr,'CTYPE1',count=type_count)
if (type_count eq 0) then begin
   print,'RFLAM: Warning - CTYPE1 not specified.'
   type = 'UNKNOWN'
endif

if (type eq 'MULTISPE') then begin

   ;;; Compute wavelengths for 'MULTISPE' type

   ; Read in all the WAT2_??? keywords into one long string.
   i = 1
   done = 0
   all_wat2 = ''
   while not(done) do begin
      if (i le 9) then istring = '00'+strtrim(i,2)
      if (i ge 10 and i le 99) then istring = '0'+strtrim(i,2)
      if (i ge 100 and i le 999) then istring = strtrim(i,2)
      if (i ge 1000) then begin
         print,'RFLAM: Too many WAT2_??? keywords.  Wavelengths not computed.'
         return,data
      endif
      temp_wat2 = sxpar(hdr,'WAT2_'+istring,count=wat2_count)
      if (wat2_count ne 0) then begin
         all_wat2 = all_wat2 + temp_wat2
         i = i+1
      endif else begin
         done = 1
      endelse
   endwhile
   ; Compute wavelengths for each row.  The entry for each row
   ; should start with spec? = "Some# 1 0 lam0 dlam ...", where
   ; Some# is the order number(?).  The entries should proceed
   ; from red to blue.
   sz = size(data)
   n_orders     = sz(2)
   pixperorder  = sz(1)
   pixarr       = dindgen(pixperorder)
   lambda       = dblarr(pixperorder,n_orders)
   lastposition = 0
   for i=1,n_orders do begin
      entry_start = strpos(all_wat2,'spec'+strtrim(i,2),lastposition)
      lam0start   = 4 + strpos(all_wat2,'1 0 ',entry_start)
      dlamstart   = 1 + strpos(all_wat2,' ',lam0start)
      dlamend     = -1 + strpos(all_wat2,' ',dlamstart)
      lam0end     = dlamstart - 2
      lam0length  = lam0end - lam0start + 1
      dlamlength  = dlamend - dlamstart + 1
      lam0string  = strmid(all_wat2,lam0start,lam0length)
      dlamstring  = strmid(all_wat2,dlamstart,dlamlength)
      lam0        = double(lam0string)
      dlam        = double(dlamstring)
      lambda(*,i-1) = lam0 + dlam*pixarr
      lastposition = dlamend
   endfor

endif else begin

   ;;; Compute wavelength for all other types

   refpix = sxpar(hdr,'CRPIX1',count=refpix_count)
   if (refpix_count eq 0) then begin
      print,'RFLAM: Warning - CRPIX1 not specified.  Using 1.'
      refpix = 1
   endif
   lam0   = sxpar(hdr,'CRVAL1',count=lam0_count)
   delta  = sxpar(hdr,'CDELT1',count=delta_count)
   ; Alternate wavelength spacing keyword
   if (delta eq 0) then delta = sxpar(hdr,'CD1_1',count=delta_count)
   
   ;;; Make sure enough information is available to compute wavelengths
   if (lam0_count eq 0 or delta_count eq 0) then begin
   
      lambda = -1
   
   endif else begin
   
      ;;; Make sure that the reference pixel is the first pixel
   
      if (refpix ne 1) then print,'Warning: reference pixel is not first pixel.'
   
      ;;; Construct wavelength array
   
      type  = strtrim(type,2)
      lam0  = double(lam0)
      delta = double(delta)
   
      n     = n_elements(data)
      wave  = dblarr(n)
   
      for i=long(0),long(n-1) do wave(i) = lam0 + (i - (refpix-1))*delta
   
      knowtype = 0.  ;;; Default is 'CTYPE1' value not recognized
   
      ;;; Is the wavelength scaling most likely linear?  If so,
      ;;; then overide type.
      
      if (lam0 gt 10) then type = 'LAMBDA'
        
      ;;; Default is that pixel size is constant in wavelength
      
      if (type eq 'LAMBDA' or type eq 'PIXEL') then begin 
         lambda    = wave
         knowtype  = 1
      endif
      
      ;;; Pixel size may alternatively be constant in velocity
      
      if (type eq 'LINEAR') then begin
         lambda   = 10^wave
         knowtype = 1
      endif
      
      if not(knowtype) then begin
         print,'Error: Unknown value of CTYPE1 in header.'
         return,-1
      endif
      
      ;;; Convert from inverse centimeters to angstroms (optional)
      
      if keyword_set(invcm) then begin
         lambda = 1d8 / lambda
         ; Need to reverse order?
         if (lambda(0) gt lambda(1)) then begin
            lambda = reverse(lambda)
            data   = reverse(data)
         endif
      endif
      
      ;;; Convert wavelength array to single precision (optional).
      if keyword_set(float) then lambda=float(lambda)
   
   endelse

endelse

;;; Convert wavelengths to vacuum (optional)
if keyword_set(vacuum) then airtovac,lambda

return,data
end
