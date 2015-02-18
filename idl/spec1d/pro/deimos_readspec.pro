;+
; NAME:
;   deimos_readspec
;
; PURPOSE:
;   generic routine for reading deimos 1-D outputs
;
; CALLING SEQUENCE:
;   deimos_readspec, maskname, slitname, silent=silent, $
;     boxcar=boxcar,topdir=topdir, zans=zans, $
;     lambda=lambda, spec=spec, ivar=ivar, zall=zall
; 
; INPUTS:
;   maskname   - mask name (string)
;   slitname   - slit name (string)
;
; OPTIONAL INPUTS:
;   topdir     - if set, use this instead of $D2_RESULTS
;  	
; KEYWORDS:
;   silent     - read silently
;   boxcar     - use the boxcar extraction instead of the optimal (default)
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   spec       - spectrunm  [counts]
;   lambda     - wavelength [angstrom]
;   ivar       - inverse variance [same units as spec]
;   zans       - structure output by reduce1d
;   zall       - all redshifts and slitnames in file
;
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
;   still lots of hardwired evil until we have run the 1D code on 
;    enough stuff to test this code.
;
; REVISION HISTORY:
;   2002-Oct-01  Written by D. Finkbeiner
;   2002-Oct-09  added synspec stuff
;   2002-Oct-18  DSM - hardwired to v0_9
;   2002-Oct-23  AC - added the boxcar keyword and changed default to
;                     read the optimal extraction; also looks in current
;                     directory first for the file before going to d2_results
;----------------------------------------------------------------------
pro deimos_readspec, maskname, slitname, silent=silent, $
                     boxcar = boxcar, topdir=topdir, zans=zans, $
                     lambda=lambda, spec=spec, ivar=ivar, $
                     synspec=synspec, zbest=zbest, fulldir=fulldir

  if keyword_set(boxcar) then boxcar = 1 else boxcar = 0
  if NOT keyword_set(topdir) then $
    topdir = getenv('D2_RESULTS')
  if topdir eq '' then print, 'You should set $D2_RESULTS'


  deep_dir = getenv('DEEP_DIR')
  if deep_dir eq '' then print, 'You should set $DEEP_DIR'

; -------- read zresult file
  fname = 'zresult.'+maskname+ '*.fits'
  zdir = concat_dir(topdir, 'zresult')
  zbest = mrdfits(concat_dir(zdir, fname), 1, /silent) 
;1st HDU is best, 2nd is all
  w = where(long(zbest.slitname) EQ long(slitname), nslit)
  if nslit EQ 0 then begin
     print, 'No slit number ', slitname
     spec = 0
     return
  endif 
  if nslit GT 1 then print, 'More than one zresult entry with that slit number!  Taking first...'
  zans = zbest[w[0]]

; -------- get slitobj file
  slitformat = '(I3.3)'
  slitstr = slitname
  if size(slitname, /tname) EQ 'STRING' then begin 
     if strlen(slitname) LT 3 then $
       slitstr=string(long(slitname), format=slitformat)
  endif else begin 
     slitstr=string(slitname, format=slitformat)
  endelse

; version = whatever -- or can get from topdir
;  datadir = concat_dir(concat_dir(topdir, version), maskname)

; DSM--- hardwired version v0_9
;********************************
;  datadir = concat_dir(topdir, maskname)
;  datadir = concat_dir(concat_dir(topdir, maskname), 'v0_9')
;  filelist=findfile(getenv('D2_RESULTS')+'/'+maskname+'/*/*.plan')
  filelist=findfile(topdir+'/'+maskname+'/*/obj_info*.fits*')
  filelist=filelist[0]
  if keyword_set(fulldir) then begin
      filelist = findfile(concat_dir(fulldir[0], 'obj_info*.fits*'))
      filelist = filelist[0]
  endif
  dirlist = strsplit(filelist, '/', /extract)
;  position=strpos(filelist,maskname,/reverse_search)
;  date=strmid(filelist,position-10,9)
  ndir = n_elements(dirlist)
  date = dirlist[ndir-2]
  datadir=concat_dir((topdir+'/'+maskname),date)
  if keyword_set(fulldir) then datadir = fulldir[0]

  filespec = 'spec1d.'+maskname+'.'+slitstr+'*.fits'
  
; look first in the current directory for the file, if not found
;      then use the D2_RESULTS directory

  specname = findfile(djs_filepath(filespec, root_dir='./'), count=nfile)
  if nfile eq 0 then specname = findfile(djs_filepath(filespec, root_dir=datadir), count=nfile)
  if nfile EQ 0 then begin 
     print, 'file not found: '+filespec
     spec = 0
     return
  endif 

; -------- read in data with fill_gap
  data = fill_gap(specname[0],/tweak);,  optimal = (boxcar eq 0))

; -------- correct A band

  hdr = HEADFITS(specname[0], ext=1, /SILENT) 
  airmass = sxpar(hdr, 'AIRMASS')            
  remove_telluric, data, airmass
      
  spec = data.spec
  lambda = data.lambda
  ivar = data.ivar

; ------ convert wavelengths from air to vacuum
  lambda = data.lambda
  AIRTOVAC, lambda
  data.lambda = float(lambda)

; -------- synspec (synthesized from eigenspectra and PCA fit)  
  if arg_present(synspec) then begin 
     eigen_name = concat_dir(concat_dir(deep_dir, 'spec1d/templates'), $
                             zans.tfile)
     eigen_name = strcompress(eigen_name, /rem)
     eig = mrdfits(eigen_name, 0, h, /silent)
     w = where(zans.tcolumn GE 0, nw)
     syn = 0

  ; Temporary fix to deal with stars
  IF nw LT 1 THEN BEGIN
    nw = 1
    w = [0]
    zans.tcolumn[w[0]] = 10
  ENDIF


     for i=0, nw-1 do begin
       jj = zans.tcolumn[w[i]]
       syn = syn+zans.theta[i]*eig[*, jj]
     endfor 
     IF zans.npoly gt 0 THEN begin ;add in the polynomial terms
        npoints = n_elements(eig[*, 0])
        parray = poly_array(npoints, zans.npoly)
        for i=0, zans.npoly-1 do syn = syn + zans.theta[nw+i]*parray[*, i]
     ENDIF

     coeff0 = sxpar(h, 'COEFF0')
     coeff1 = sxpar(h, 'COEFF1')
     
     elambda = 10.d^(coeff0+dindgen((size(eig, /dimens))[0])*coeff1)
; ??? interpol is pretty stupid - use combine1fiber in idlspec2d (SDSS)
     synspec = interpol(syn, elambda, lambda/(zans.z+1))
  endif 

  return
end


