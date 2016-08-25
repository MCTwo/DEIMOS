;+
; NAME:
;   deimos_trace_crude
;
; PURPOSE:
;   Make a crude trace of the slit "edges" on a flat
;
; CALLING SEQUENCE:
;   deimos_trace_crude, im, tset1, tset2, imgbin, $
;        invvar=invvar, nmed=nmed, ybin=ybin, ncoeff=ncoeff
;
; INPUTS:
;   im         - 2D spectral image
;
; OPTIONAL KEYWORDS:
;   invvar     - inverse variance of image
;   nmed       - y-direction median of nmed pixels before detection
;   ybin       - bin in y by a factor of ybin to speed up.  [4]
;   ncoeff     - number of coefficients to use in traceset fit. 
;   gap        - typical gap between slits, in pix                
;
; OUTPUTS:
;   tset1      - trace set structure for "left" traces
;          tset1 = $
;            { func    :    'legendre'  , $
;              xmin    :    0  , $
;              xmax    :    4095   , $
;              coeff   :    array[ncoeff, ntrace] $
;            }
;   tset2      - trace set structure for "right" traces
;   imgbin     - image binned in y by a factor of ybin
;
;
; COMMENTS:
;   An inverse variance array should be passed, but if it is not the
;   procedure pretends it was all 1.  
;
;   This routine makes heavy use of the SDSS trace_crude procedure and 
;     several traceset handling routines by D. Schlegel and S. Burles.
;
; BUGS:
;   treats invvar like a bad pixel mask (1=good 0=bad) for now. 
;
;   Because of the way trace_crude works, we don't want to pass it
;   ivar=0, but rather ivar = small for the bad columns.  If we
;   interpolate over them with djs_mask_interp the trace stays on
;   track, but the errors are appropriately larger. 
;
; EXAMPLES:
;
; REVISION HISTORY:
;   01-Dec-2000  Written by D. Finkbeiner, Berkeley
;   12-Jul-2002  demand derivative of ratio to running median also be
;      above a certain threshold - DPF 
;   17-Sep-2002  swap tset1, tset2 to compensate for earlier sign
;                 error - MD
;   14-Oct-2002  Test if invvar_in set instead of invvar
;                 and use edgeivar (with nonzero floor) - DPF
;   
;-
;------------------------------------------------------------------------------
pro deimos_trace_crude, im, tset1, tset2, imgbin, $
        invvar=invvar_in, nmed=nmed, ybin=ybin, ncoeff=ncoeff, model=model, $
                        gap=gap

  nx = (size(im, /dimens))[0]
  ny = (size(im, /dimens))[1]

  if not keyword_set(ncoeff) then ncoeff = 4
  if not keyword_set(ybin) then ybin = 4   ; y bin factor
  if n_elements(gap) eq 0 then gap=4

  imgbin = rebin(im, nx, ny/ybin)

  if not keyword_set(invvar_in) then invvar = $
    byte(imgbin-imgbin)+1B $      ;set to one for now
    else invvar = rebin((invvar_in NE 0), nx, ny/ybin)
  if not keyword_set(nmed) then nmed = 9 ;median nmed pixels in disp. dir. 
     
; Median filter the entire image along columns by NMED rows
  imgmed = imgbin
  if (nmed GT 1) then $   
    for ix=1, nx-1 do imgmed[ix,*] = median(transpose(imgmed[ix,*]), nmed)

; Here we are assuming that the slits are not more than about
;   35/2 pix apart...

  ystart = ny/ybin/2


    useful = where(smooth(reform(imgmed[*, ystart]), 101) gt 1000., usect)
    if usect lt 100 then useful = lindgen(nx)
    useflag=intarr(nx)
    useflag[useful]=1
   

; make sure we still can find derivatives if the mask doesn't extend all
; the way to edge
  baselevel=djs_median(imgmed[*,ystart], width=8*fix(gap)+1,boundary='reflect')
  print,'Smoothing to fix base level: ',8*fix(gap)+1,' pix'

  if usect lt n_elements(baselevel) then begin 
      lessuseful=1-dilate(useflag eq 0,intarr(101)+1)
      lessuseful=lessuseful*useflag
      if total(lessuseful) gt 100 and total(lessuseful eq 0) gt 0 then $
    baselevel[where(lessuseful eq 0)]=min(baselevel[where(lessuseful eq 1)])
  endif

  rat = imgmed[*,ystart]/baselevel
  drat = shift(rat, -1)-shift(rat, 1)
  drat_cut = 0.15 ; 15 percent change per 2pix required. 

  edge = shift(imgmed, -1)-shift(imgmed, 1)
  edgeivar = float(shift(invvar, -1) AND shift(invvar, 1))
  edgeivar = edgeivar > 1.0E-3


 
  yrow = edge[*, ystart]

   thresh = 5*djsig(yrow[useful])
   ;print, yrow[useful]
   if thresh gt 10000 then thresh = 7800.		; Added BCL 7/29/16, for issues with flats that have >55000 counts near edges
  print, 'WARNING: Threshold counts for separation between slits set very high, likely spurious high values near edge of chip, trying a more reasonable value...'
  print, 'Slit edge trace thresh: ', thresh



; Get the "crude" trace
  xset = trace_crude((-edge)>0, edgeivar, ystart=ystart, radius=4, $
                   thresh=thresh, yset=yset, xerr=xerr, maxshifte=0.05*ybin)
  ind = xset[ystart, *]

  dratcompare=drat[round(ind)]<drat[round(ind)-1]
  dratcompare=dratcompare<drat[round(ind)+1]
  w = where(dratcompare LT -drat_cut, ngood)


  if ngood eq 0 then message, 'we are in deep trouble...'
  xset = xset[*, w]
  xerr = xerr[*, w]
  yset = yset[*, w]
  xy2traceset, yset*ybin, xset, tset2, ncoeff=ncoeff, invvar=1./xerr^2, $
    xmax=ny-1,/SILENT
;reversed, 9/17/02 MD



; -------- 
  xset = trace_crude(edge>0, edgeivar, ystart=ystart, radius=4, $
                   thresh=thresh, yset=yset, xerr=xerr, maxshifte=0.05*ybin)
  ind = xset[ystart, *]


 dratcompare=drat[round(ind)]>drat[round(ind)-1]
  dratcompare=dratcompare>drat[round(ind)+1]
  w = where(dratcompare GT drat_cut, ngood)


  if ngood eq 0 then message, 'we are in deep trouble...'
  xset = xset[*, w]
  xerr = xerr[*, w]
  yset = yset[*, w]
  xy2traceset, yset*ybin, xset, tset1, ncoeff=ncoeff, invvar=1./xerr^2, $
    xmax=ny-1,/SILENT
;reversed, 9/17/02 MD

; Convert x and y positions to a traceset (tset) structure
;  in this case we want to fit x=f(y)


  return
end
