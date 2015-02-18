;+
; NAME:
;  genrandom
;
; PURPOSE:
;  generates a random set of points in a given pcat field, with mask
;  for observed region for holes in pcat coverage. 3d version if
;  phi_r is passed
;
; CALLING SEQUENCE:
;   random = genrandom(nran,pcatfield,[windowf=windowf,phi_r=phi_r,$ 
;            ralim=ralim, declim=declim])
; 
; INPUTS:
;   nran - number of random points desired
;   pcatfield -- 2 digit pcat field number, eg. 11 or 42
;
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;   windowf -- if specified, the window function of observation, from 
;              makewindowf.pro
;   declim  -- 2vector needed with windowf
;   ralim   -- 2vector needed with windowf
;   phi_r   -- radial selection function.  If not provided,
;                         all points will have z=0
;              (phi_r is a structure with tags 'z', 'r' (comoving
;              distance)  and 'phir')
;
; OUTPUTS:
;   random  -- structure containing relevant information for random catalog
;
; COMMENTS:
;
; REVISION HISTORY:
;   created Oct00 md
;   modified to 3d with detailed window function, 11dec02
;
;----------------------------------------------------------------------

function genrandom, nran,photfile, windowf=windowf,phi_r=phi_r,$ 
               ralim=ralim, declim=declim


  if n_elements(photfile) eq 0 then photfile=42
  photo = getphot(photfile, outlabel, header)
  filestring = string(photfile, format='(i2)')

  s1=string(photfile/10, format='(i1)')
  s2=string(photfile mod 10, format='(i1)')

  directory=getenv('IDLDEEPPHOTOMETRY_DIR') +'field' +s1 + '/'+ s2+ '/quilt/'
;mask center-  get from mask database
  files = findfile(directory  + 'R.000.fits.Z', count=ct)
  if (ct EQ 0) then begin
         message, 'No data files found in path'
         return, 0 
  endif

  hdr = headfits(files[0])
  extast, hdr,  astr ; get WCS information for .000 file (lower left)
  coh =  299792.458/100.


  tmp=minmax(photo.xs(0))
  maxra=tmp(1)+.02
  minra=tmp(0)-.02
  tmp=minmax(photo.xs(1))
  maxdec=tmp(1)+.02
  mindec=tmp(0)-.02

  seed=-70991
  object=create_struct('ra',0.,'dec',0.,'xs',[0.,0.],'flag',0b, $
          'z', 0.0, 'rz', 0.0, 'rx', 0.0, 'ry', 0.0, 'sf', 0.0, 'weight', 0.0) 
  rand=replicate(object,nran) ;copy structure definition

;  delta=5.75e-05 
  racen = astr.crval[0]  ;center of field
  deccen = astr.crval[1]

  rand.xs[0]=randomu(seed,nran)*(maxra-minra) + minra ;ra distribution
  rand.xs[1]=randomu(seed,nran)*(maxdec-mindec) + mindec 
  rand.dec =  +rand.xs[1] +deccen
  rand.ra  =  -rand.xs[0]/cos(rand.dec/!radeg) +racen  ;get ra, dec

;
; mask out objects in position of bad spots in pcat file
; return list of acceptable objects 
  flag=maskquick(rand.xs,filestring, 'BRI') 

  rand.flag=flag
  rand = rand[where(rand.flag le 1)] ; keep objects defined in 2 or more bands
  nnow = n_elements(rand)
  rand.weight = 1. ;weight factor for correlations (not needed for randoms)

;
; use supplied window function to restrict randoms to observed areas
; 
  if keyword_set(windowf) then begin 
    if NOT keyword_set(ralim) then message, 'You must specify ralim'
    if NOT keyword_set(declim) then message, 'You must specify declim'
    nra = (size(windowf, /dimen))[0]
    ndec = (size(windowf, /dimen))[1]
    ra_inc = (ralim[1]-ralim[0])/nra
    dec_inc = (declim[1]-declim[0])/ndec

; first trim the windowf to region spanned by pcat data by defining a
; polygon of the valid region
    pcatra = minmax(photo.ra)
    pcatdec = minmax(photo.dec)
    xrange = [pcatra[0], pcatra[1], pcatra[1], pcatra[0] ] ;limits of pcat RA
    yrange = [pcatdec[0], pcatdec[0], pcatdec[1], pcatdec[1] ] ;dec limits
    xwindow = (fix((xrange - ralim[0])/ra_inc)  < nra-1) >  0 ;get mask x,y
    ywindow = (fix((yrange - declim[0])/dec_inc) <  ndec-1) >  0
    keep = polyfillv(xwindow, ywindow, nra, ndec) ;get valid range
    temp = windowf*0. 
    temp[keep] =  temp[keep] +1. ;set region to retain
    windowf = windowf*0. + windowf*temp  ;set all bad regions to 0
    delvarx, temp ;release memory

; now apply windowf to randomly drawn points
    xwindow = (fix((rand.ra - ralim[0])/ra_inc)  < nra-1) >  0 ;get mask x,y
    ywindow = (fix((rand.dec - declim[0])/dec_inc) <  ndec-1) >  0
    weight = windowf(xwindow, ywindow) ;prob of including each point
    keep = weight gt randomu(seed, nnow) ;get next set of random numbers
    rand = rand[where(keep eq 1)] ;keep only objects in region
  endif

  if keyword_set(phi_r) then begin  ;set redshift distribution if provided a
; selection function
    nnow = n_elements(rand)
    nz = n_elements(phi_r.z)  
    zmin = phi_r.z[0]
    zmax = phi_r.z[nz-1]
    rmin = phi_r.r[0]
    rmax = phi_r.r[nz-1]
;    dr = phi_r.r[1] - phi_r.r[0] ;increment in r
    cvolume = phi_r.r*0. ;cumulative volume, starting at 0 for 1st point

; selectable volume in the survey
    for i=1, nz-1 do cvolume[i] = cvolume[i-1] + .5* $
       ( (phi_r.r[i-1])^2*phi_r.phir[i-1] + (phi_r.r[i])^2*phi_r.phir[i])
    cvolume = cvolume/max(cvolume) ;normalize to cumulative fraction
    sss = spl_init(cvolume,phi_r.r) 
;setup spline interpolation of r versus volume
    rand.rz = spl_interp(cvolume, phi_r.r, sss, randomu(seed, nnow)) 
; points should be uniformly distributed in dV
    rand.z  = r2z(rand.rz) ;convert to z, model dependent but setup prior to
                         ;running this routine 
; next do spline in reverse, so as to get sf for each point
    sssr = spl_init(phi_r.r, phi_r.phir)
    rand.sf =  spl_interp(phi_r.r,  phi_r.phir, sssr, rand.rz) 
    rand.rz = rand.rz*coh ;convert distance to Mpc


  endif

return, rand
end
















