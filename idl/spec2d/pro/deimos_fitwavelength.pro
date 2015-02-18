; NAME:
;   deimos_fitwavelength
;
; PURPOSE:
;   finds wavelengths in one 2d rectified slitlet
;
; CALLING SEQUENCE:
;    deimos_fitwavelength, rect_arc,rect_arcivar,lamps,slitcoords, $
;       model_lambda, wave,  [ flat=flat], $
;       [polyflag=polyflag, polyx=polyx, s_coeff=s_coeff, dlam= ]
;
; INPUTS:
;    rect_arc -- 2-d spectrum of arc lamp
;    rect_arcsat -- 2-d saturation mask (1=bad) for rect_arc
;    rect_arcivar -- invvar of same
;    slitcoords -- structure detailing bluslit information for this slitlet
;    model_lambda -- optical model wavelength solution for this slitlet
;    lamps  -- structure detailing spectral linelist
;    chipno -- which chip number? (1-8)
;    grating -- grating [lines/mm]
;    lambda_c -- central wavelength [Ang]
;    
; KEYWORDS:
;    flat   -- if set, normalize by flat
;    anamorph -- anamorphic factor, defaults to 1.6
;    polyflag -- set if using special polynomial expansion of 2d
;                Lambda
;    ymid     -- row to use as 'central'
;    (other keywords passed to deimos_arcfit)
;
; OUTPUTS:
;   wave   -- lambda for each pixel in 2-d image
;   polyx  -- Legendre Polynomial coefficients if polyflag set
;   s_coeff-- slit tilt coefficients if polyflag set
;   dlam   -- mean lambda shift per row
;   sigma  -- RMS uncertainty in wavelength solution (in Angstroms)
; COMMENTS:
;   called from deimos_spslit, or from quicklook 
;    
;
; REVISION HISTORY:
;
;       Wed Feb 20 17:52:02 2002, Douglas Finkbeiner (dfink)
;             Split of from 2dtest, from 2001-Aug-22
;       15apr02 split off from deimos_spslit MD
;
;
;----------------------------------------------------------------------


pro deimos_fitwavelength, rect_arc, rect_arcsat, rect_arcivar, lamps, $
          chipno, grating, $
          lambda_c, slitcoords, model_lambda, wave,  flat=flat, $
          anamorph=anamorph, plot=plot, polyflag=polyflag,  $
	ymid=ymid, $
          polyx=polyx, s_coeff=s_coeff, dlam = dlam, dline=dline, $
          dlsig=dlsig, wset=wset, sigma=sigma, dirtyarc=dirtyarc

  if NOT keyword_set(anamorph) then begin 
     print, '  WARNING:  You must set anamorphic factor!'
     anamorph = 1.6
  endif 

;!except=2

  if n_elements(dirtyarc) eq 0 then dirtyarc=rect_arc

  scale_mm_per_asec = 0.73
  scale_pix_per_asec = 8.52
  
  wid_asec = slitcoords.slitwidth/scale_mm_per_asec
  print, 'Slit width [asec]: ', wid_asec
  wid_pix = wid_asec*scale_pix_per_asec/anamorph

  
  sizey=(size(rect_arc))[2]
  sizex=(size(rect_arc))[1]
  ymid = sizey/2

  ignorevig = (rect_arcsat eq 4b OR rect_arcsat eq 6b)

  ;help, rect_arcivar, ignorevig
; don't use the middle row if e.g. it's on a bad column, or the
; vignetted regions at the ends of the mask
  nbadpix=total((rect_arcivar eq 0) AND (ignorevig eq 0),1)
  srt=nbadpix[sort(nbadpix)]

  ;print, 'The number of bad pixels are', nbadpix		;Added BL 5/13

  whok=where(nbadpix le ((srt[0.3*sizey]+30) < sizex/2),okct)
  if okct gt 0 then begin
	minfrommid=min(abs(whok-ymid),newymid) 
	print, 'Minimum offset is ', min(whok-ymid)		;Added BL 5/13
	print, 'The y slitsize is ', sizey
	endif else begin
    
	minfrommid = 0

	endelse
  if minfrommid gt 0 and sizey ne 108 then print,'Using row ',whok[newymid],' instead of central row ',ymid
  if okct gt 0 then ymid=whok[newymid]
  if minfrommid gt 0 and sizey eq 108 then begin
  print, where(nbadpix lt 10)
  print,'Gotcha bitch, using row 90 instead of central row', ymid ;whok[newymid]+ymid*2,' instead of central row ',ymid
  atv, rect_arc
  spec2 = dirtyarc[*,ymid]
  spec2 = reform(spec2, n_elements(spec2))
  atvspec = dblarr(108, 4096)
  for i=0, n_elements(spec2)-1 do begin 
  	atvspec[*, i] = spec2[i]
  endfor
  atv, atvspec
  
  if okct gt 0 then ymid=90.
  endif

  print, 'The middle y is', ymid
  spec = dirtyarc[*,ymid]

  spec = reform(spec, n_elements(spec))

  if NOT keyword_set(grating) then message, 'Must set grating value!'
  deimos_wavematch, model_lambda, spec, lamps, chipno, grating, lambda_c, $
    arcline_x, slitwid=round(wid_pix),wset=wset
       ;, arcsat=rect_arcsat

  if n_elements(arcline_x) LT 3 then begin 
     message, 'ABORT: No arclines found - setting lambda model to zeros.', /info
     print,'ABORT: No arclines found - setting lambda model to zeros.'
     wave = rect_arc*0
     return
  endif 
     
; -------- Wavelength solution
     sigma =  deimos_arcfit(rect_arc, rect_arcsat, arcline_x, lamps, wset, $
          lamdif=lamdif, ncoeff=4, slitwid=wid_pix,  arcivar=rect_arcivar,  $
          plot=plot, dlam=dlam, dline=dline, dlsig=dlsig, $
          polyflag=polyflag, polyx=polyx, s_coeff=s_coeff, $
 	  maxerr=0.15*(1200./float(grating)), anamorph=anamorph, $
                          ymid=ymid)
; 	  maxerr=0.2, anamorph=anamorph,ymid=ymid )
     
;!except=1

; to try as an experiment when I find a bad slit:
;     if n_elements(polyx) le 3 then polyx=wset.coeff

badfit=0
if n_elements(polyx) lt 6 then badfit = 1
if badfit eq 0 and n_elements(polyx) gt 1 then badfit = (polyx[1] lt 0)


; if deimos_arcfit did poorly, use tweaked optical model guess
if badfit then begin
        polyx=wset.coeff
	badtilt=(model_lambda.tiltx)[0] eq 0.
        if badtilt eq 0 then $ 
          s_coeff=model_lambda.tiltx * (max(abs(model_lambda.tiltx)) lt 0.01)
        dlam=fltarr(sizey)
        dline=0.
        dlsig=0.
	; flag values: 1: used tweaked optical model; 
        ; 10: no model tilt available 
	sigma=1.+9.*badtilt
endif



return
end










