;+
;
; NAME
;      fill_gap.pro
;
; PURPOSE
;      Stitches together the red and blue portions of an object
;      spectrum. 
;
;
; SYNTAX
;      fill_gap, file1d, [nbuf=nbuf, npoly=npoly, nbins=nbins, $
;                binsz=binsz, /horne, /boxsprof, /optimal, /boxcar, $
;                /doplot, /debug, /nlsky, header=header,/tweak]
;
; INPUTS
;      file1d = a parameter giving a spec1d file name. 
;      nbuf = this optional argument allows the user to specify the
;             buffer size at the edge of the red and blue spectra for
;             which the inverse variance will be set to zero. The
;             default value is 0. That is, the inverse variance for
;             the first and last nbuf pixels of the blue and red spectra
;             is set to 1E-20.
;      npoly = this optional parameter allows the user to specify the
;              order of polynomial to fit in order to determine the
;              continuum levels for the red and blue
;              spectra. Specifically, npoly gives the number of
;              coefficients in the polynomial. That is, a zeroth-order
;              polynomial corresponds to npoly=1 and so on for higher
;              orders. The default is npoly=1.
;      nbins = this parameter in conjunction with the parameter binsz
;              allows the user to specify the pixel range at the end
;              of the blue spectrum and beginning of red spectrum
;              which will be used in determining the continuum
;              levels. In each of nbins sets of binsz pixels, a median
;              value is calculated and these are fit according to a
;              polynomial function. The default is nbins=5.
;      binsz = this optional parameter works with nbins to determine
;              the amount of pixels used to determine the continuum
;              levels of the blue and red spectra. The default is
;              binsz=100.
;
; KEYWORDS
;      boxcar = if this keyword is set, then the 1-d spectrum
;               extracted using the boxcar extraction will be pulled
;               from the spec1d file. Otherwise, fill_gap will by 
;               default work on the horne extracted spectrum.
;      optimal = if this keyword is set, then the 1-d spectrum
;                extracted using the optimal extraction will be pulled
;                from the spec1d file. 
;      boxsprof = if this keyword is set, then the boxsprof extraction
;                 is returned. 
;      horne = if this keyword is set, then the horne extraction is
;              retrieved. Note that fill_gap will by default work on
;              the horne extracted spectrum regardless of whether or
;              not this keyword is set. 
;      ----> to check what types of extractions are contained in a
;            given spec1d file, use the GODDARD routine fits_help.pro:
;          IDL> fits_help, 'spec1dfile'
;
;      nlsky = if this keyword is set, then the 1-d spectrum extracted
;              using the non-local sky subtracted data will be pulled
;              from the spec1d file.
;      doplot = if this keyword is set, then a plot is created and sent
;               to the ps/ directory in the current working
;               directory. However, if the keyword debug is set (along
;               with the plot keyword), then the plot is sent to the
;               current plotting device. The plot displays the ends of
;               the blue and red spectra, the median points, the
;               polynomial fit, and the continuum level approximations. 
;      debug = if this keyword is set along with the doplot keyword,
;              then the plot is sent to the current plotting
;              device. If debug is set without setting doplot, then
;              nothing happens!
;      telluric = if this keyword is set, then the a-band atmospheric
;                 absorption band is corrected using the routine
;                 remove_telluric.pro.
;      tweak - apply skytweak to this file's wavelengths
;      header = returns the FITS header of the spec1d file
;               corresponding to the extraction tpye requested. in
;               case of an error, the header is returned as a int(0).
;
; OUTPUTS
;      ss1d = a structure containing the combined red and blue
;             spectrum. The tags are as follows:
;           ss1d.spec = the flux values.
;           ss1d.ivar = the corresponding iverse variance values.
;           ss1d.lambda = the wavelength values in linear lambda form.
;             note that if the fill_gap routine is unable to return a
;             1-d spectrum, then a value of zero is returned.
;
; PROCEDURES CALLED 
;      djs_median
;      makearr
;      ivarsmooth
;      mcc_polyfit
;      getcolor
;      fits_open/fits_close
;      remove_telluric 
;      fix_response
; EXAMPLES
;
;
; COMMENTS
;      Call made to FITS_INFO, so system variables
;         DEFSYSV,'!TEXTOUT',1
;         DEFSYSV,'!TEXTUNIT',0
;      are set by this routine. ----> No longer true!
;
; HISTORY
;      Created August 7, 2002 by mcc.
;      Revised August 15, 2002 by mcc - modified routine to do a
;         polynomial fit to determine the continuum levels of the red
;         and blue portions of the spectrum. Also, added plotting
;         capabilities to enable easier diagnosis of problems with
;         matching continuum levels.  
;      Revised August 29, 2002 by mcc - modified to accomodate new
;         format for spec1d files. Files will now contain both the
;         optimal and boxcar extractions. I added an /optimal keyword
;         to fill_gap to designate which one is pulled from the file
;         and has its gap filled and returned.
;      Revised January 29, 2003 by mcc - modified routine be able to
;         handle the non-local sky 1-d extractions. 
;
;-

function make_ss1d, spec=spec, ivar=ivar, lambda=lambda, $
                    ssBd=ssBd, ssRd=ssRd, gapsz=gapsz
  npix = 4096
  zero1d = {spec:fltarr(npix), ivar:fltarr(npix), lambda:fltarr(npix), $
            crmask:intarr(npix), bitmask:intarr(npix), $
            ormask:intarr(npix), infomask:intarr(npix), $
            nbadpix:intarr(npix), $
            objpos:float(0), fwhm:float(0), nsigma:float(0), $
            r1:long(0), r2:long(0), $
            skyspec:fltarr(npix), ivarfudge:float(0)}
  ss1 = replicate(zero1d, 1)
  ss2 = replicate(zero1d, 1)
  struct_assign, ssBd, ss1
  struct_assign, ssRd, ss2
  crmask = [ss1.crmask,fltarr(gapsz),ss2.crmask]
  bitmask = [ss1.bitmask,fltarr(gapsz),ss2.bitmask]
  ormask = [ss1.ormask,fltarr(gapsz),ss2.ormask]
  infomask = [ss1.infomask,fltarr(gapsz),ss2.infomask]
  nbadpix = [ss1.nbadpix,fltarr(gapsz),ss2.nbadpix]
  objpos = mean([ss1.objpos, ss2.objpos])
  fwhm = mean([ss1.fwhm, ss2.fwhm])
  nsigma = mean([ss1.nsigma, ss2.nsigma])
  r1 = mean([ss1.r1, ss2.r1])
  r2 = mean([ss1.r2, ss2.r2])
  skyspec = [ss1.skyspec, fltarr(gapsz), ss2.skyspec]
  ivarfudge = mean([ss1.ivarfudge, ss2.ivarfudge])
  ss1d = {spec:spec, ivar:ivar, lambda:lambda, $
          crmask:crmask, bitmask:bitmask, ormask:ormask, $
          infomask:infomask, nbadpix:nbadpix, $
          objpos:objpos, fwhm:fwhm, nsigma:nsigma, $
          r1:r1, r2:r2, skyspec:skyspec, ivarfudge:ivarfudge}
  return, ss1d
end



function fill_gap, file1d, nbuf=nbuf, npoly=npoly, nbins=nbins, $
                   binsz=binsz, optimal=optimal, boxcar=boxcar, $
                   boxsprof=boxsprof, horne=horne, $
                   nlsky=nlsky, spectra=spectra, flats=flats, $
                   doplot=doplot, debug=debug, gapsize=gapsize, $
                   silent=silent, telluric=telluric, header=header, $
                   tweak=tweak

; check that the file1d parameter (spec1d file name) was passed.
  if n_params() lt 1 then $
    message, 'Incorrect number of paramaters!'

; check if user wants it silent.
  if n_elements(silent) gt 0 then silent = silent[0] ge 1 $
  else silent = 0

; check if we should make a telluric correction.
  if n_elements(telluric) gt 0 then telluric = telluric[0] ge 1 $
  else telluric = 0

; check which keywords were set.
  if n_elements(nlsky) then nlsky = nlsky[0] ge 1 else nlsky = 0
  if n_elements(boxcar) then boxcar = boxcar[0] ge 1 else boxcar = 0
  if n_elements(optimal) then optimal = optimal[0] ge 1 else optimal = 0
  if n_elements(boxsprof) then boxsprof = boxsprof[0] ge 1 else boxsprof = 0
  if n_elements(horne) then horne = horne[0] ge 1 else horne = 0
  if (optimal or boxsprof or boxcar) eq 0 then horne = 1
  if n_elements(tweak) eq 0 then tweak=0


; determine data directory
  tmp=str_sep(file1d,'/') 
  nparts=n_elements(tmp)
  path=''

  for i=0,nparts-2 do path=path+tmp[i]+'/'

  if strlen(path) eq 0 then cd,'.',current=path


; open the spec1d FITS file.
  fits_open, file1d, fcb
  extnames = strcompress(fcb.extname, /rem)
  number_of_extensions = n_elements(extnames)
; find the desired FITS extension(s) based on the keywords set by the
; user.
  if nlsky then begin
      if boxcar then dex = where(strpos(extnames, 'Boxcar-NL') ne -1, extnum)
      if boxsprof then dex = where(strpos(extnames, 'Bxspf-NL') ne -1, extnum)
      if optimal then dex = where(strpos(extnames, 'Optimal-NL') ne -1, extnum)
      if horne then dex = where(strpos(extnames, 'Horne-NL') ne -1, extnum)
  endif else begin
      if boxcar then dex = where(strpos(extnames, 'Boxcar') ne -1 and $
                                 strpos(extnames, 'NL') eq -1, extnum)
      if boxsprof then dex = where(strpos(extnames, 'Bxspf') ne -1 and $
                                 strpos(extnames, 'NL') eq -1, extnum)
      if optimal then dex = where(strpos(extnames, 'Optimal') ne -1 and $
                                 strpos(extnames, 'NL') eq -1, extnum)
      if horne then dex = where(strpos(extnames, 'Horne') ne -1 and $
                                 strpos(extnames, 'NL') eq -1, extnum)
  endelse
; close the FITS file.
  fits_close, fcb

; if the desired extraction type (extension) was not found, then
; print an error message and return a value of 0.
  if extnum eq 0 then begin
      if nlsky and not(silent) then begin
          if boxcar then print, '(fill_gap.pro) No non-local sky, ' + $
            'boxcar extraction found!' 
          if boxsprof then print, '(fill_gap.pro) No non-local sky, ' + $
            'boxsprof extraction found!'
          if optimal then print, '(fill_gap.pro) No non-local sky, ' + $
            'optimal extraction found!'
          if horne then print, '(fill_gap.pro) No non-local sky, ' + $
            'horne extraction found!'
      endif else begin
          if not(silent) then begin
              if boxcar then $
                print, '(fill_gap.pro) No boxcar extraction found!' 
              if boxsprof then $
                print, '(fill_gap.pro) No boxsprof extraction found!'
              if optimal then $
                print, '(fill_gap.pro) No optimal extraction found!'
              if horne then $
                print, '(fill_gap.pro) No horne extraction found!'
          endif
      endelse
; read the header of the spec1d file.
      if number_of_extensions gt 0 then $
        header = headfits(file1d, ext=1, /silent) $
      else header = 0
      return, 0
  endif

; print an error message if > 2 extensions were found!
  if extnum gt 2 then begin
      if not(silent) then $
        print, '(fill_gap.pro) More than 2 extensions found!!!'
  endif
  
; read-in the blue and red portions of the 1-d spectra from the spec1d
; file. also, kill off any NaN's before they are a problem!
  ss1 = mrdfits(file1d, dex[0], hdrB, /silent)
  header = hdrB

  if tweak then ss1.lambda=applytweaks(ss1.lambda,hdrB,path)

  notfinite = where(finite(ss1.spec)*finite(ss1.ivar) eq 0, nct)
  if (nct gt 0) then begin
      ss1.spec[notfinite] = 0.
      ss1.ivar[notfinite] = 0.
  endif
  if extnum gt 1 then begin
     ss2 = mrdfits(file1d, dex[1], hdrR, /silent)

       if tweak then ss2.lambda=applytweaks(ss2.lambda,hdrR,path)
     notfinite = where(finite(ss2.spec)*finite(ss2.ivar) eq 0, nct)
     if (nct gt 0) then begin 
        ss2.spec[notfinite] = 0.
        ss2.ivar[notfinite] = 0.
     endif
  endif

; now take the first nbuf and last nbuf pixels in each spectrum and set
; their inverse variance to ~zero (1E-20). this is done to de-weight
; any spikes found at ends of slits. note that the default value for
; nbuf is 0 pixels.
  if keyword_set(nbuf) then nbuf = nbuf[0] else nbuf = 0
  n1 = n_elements(ss1.ivar)
  if nbuf gt 0 then begin
      ender = fltarr(nbuf) + 1E-20
      ss1.ivar = [ender, ss1.ivar[nbuf:n1-nbuf-1], ender] 
      if extnum gt 1 then begin
          n2 = n_elements(ss2.ivar)
          ss2.ivar = [ender, ss2.ivar[nbuf:n2-nbuf-1], ender] 
      endif
  endif
; try to find bogus spike in spectrum occuring away from the ends. if
; found, set ivar to ~zero.
;  medflx = djs_median(ss1.spec, width=5, boundary='reflect')
;  badpts = where(abs(ss1.spec - medflx) gt 10/sqrt(median(ss1.ivar)), bcnt)
;  if bcnt gt 0 then ss1.ivar[badpts] = 1E-20
;  if extnum gt 1 then begin
;    medflx = djs_median(ss2.spec, width=5, boundary='reflect')
;    badpts = where(abs(ss2.spec - medflx) gt 10/sqrt(median(ss2.ivar)), bcnt)
;    if bcnt gt 0 then ss2.ivar[badpts] = 1E-20
;  endif
  
; if a sky-line falls on the end of the slit, then mask that region by
; setting ivar to zero there. recall the the sky spectra has been
; normalized to counts per pixel per hour.
  skyline = 400.0
  sb = 10
  n1 = n_elements(ss1.spec)
  if max(ss1.skyspec[0:sb] gt skyline) then begin
      if not(silent) then $
        print, '(fill_gap.pro) Masking sky line residual at end of slit...'
      ss1.ivar[0:sb] = 0.0
  endif
  if max(ss1.skyspec[n1-1-sb:n1-1]) gt skyline then begin
      if not(silent) then $
        print, '(fill_gap.pro) Masking sky line residual at end of slit...'
      ss1.ivar[n1-1-sb:n1-1] = 0.0
  endif
  if extnum gt 1 then begin
      n2 = n_elements(ss2.spec)
      if max(ss2.skyspec[0:sb]) gt skyline then begin
          if not(silent) then $
            print, '(fill_gap.pro) Masking sky line residual at end of slit...'
          ss2.ivar[0:sb] = 0.0
      endif
      if max(ss2.skyspec[n2-1-sb:n2-1]) gt skyline then begin
          if not(silent) then $
            print, '(fill_gap.pro) Masking sky line residual at end of slit...'
          ss2.ivar[n2-1-sb:n2-1] = 0.0
      endif
  endif

; remove the a-band atmospheric absorption.
  if telluric then begin
      airmass = sxpar(hdrB, 'AIRMASS')
      remove_telluric, ss1, airmass,silent=silent
      fix_response,ss1
      if extnum gt 1 then begin 
          remove_telluric, ss2, airmass,silent=silent
	  fix_response,ss2
      endif
  endif

; make sure that both sides of the spectrum are good. one might be
; filled with all zeros due to an error in extract1d.
  if extnum gt 1 then begin
      if total(ss2.lambda) eq 0 and $
        total(ss1.lambda) eq 0 then begin
          print, 'No valid 1-d data for this slit!'
          return, 0.
      endif
      if total(ss2.lambda) eq 0 then extnum = 1
      if total(ss1.lambda) eq 0 then begin
          extnum = 1
          ss1 = ss2
      endif
  endif

; check that the wavelength values are positive!
  ltz1 = where(ss1.lambda lt 0.0, ltz1cnt)
  if extnum gt 1 then ltz2 = where(ss2.lambda lt 0.0, ltz2cnt) $
  else ltz2cnt = 1
  if ltz1cnt gt 0 and ltz2cnt gt 0 then begin
      if not(silent) then print, 'Negative wavelength values!'
      return, 0 
  endif
  if ltz1cnt gt 0 and ltz2cnt eq 0 then begin
      if not(silent) then print, 'Negative wavelength values!'
      ss1 = ss2
  endif
  if ltz1cnt eq 0 and ltz2cnt gt 0 then begin
      if extnum gt 1 then print, 'Negative wavelength values!'
      extnum = 1
  endif

; if there are two spectra in spec1d file, then stitch them 
; together. otherwise just return the single spectrum.
  if extnum eq 1 then begin
      gapsize = 0
; check that the spec1d structure isn't filled with all zeros. recall
; that do_extract will return a structure of all zeros if it runs into
; a snag.
      if n_elements(uniq(ss1.lambda))-1 le 0 then begin
          if not(silent) then print, 'No unique lambda values!!'
          return, 0
      endif
      return, ss1
  endif else begin
; the stitching of the spectra can be done 3 ways. the first way
; (/flats) is by using hard-coded values to rescale the red and blue
; spectra. these values are based upon measuring the relative counts
; in the flat frames. this is done assuming that the dominant effect
; in scaling differences is due to efficiency of the detector which is
; best measured with the strong counts in the flats. and the second
; way (/spectra) uses the spectra to measure the continuum level at
; the gap and scale the spectra to match there. the third way is the
; default setting in which no re-scaling is done prior to stitching
; the spectra together.

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%% /FLATS %%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if keyword_set(flats) then begin
          chipB = STRING(SXPAR(hdrB, 'CHIPNO'), format='(I1.0)')
          chipR = STRING(SXPAR(hdrR, 'CHIPNO'), format='(I1.0)')
          case chipB of
              '1':begin
                  if chipR ne '5' then $
                    print, '(fill_gap.pro) ERROR: Chips do not match!'
                  chip_ratio = 1.05
              end
              '2':begin
                  if chipR ne '6' then $
                    print, '(fill_gap.pro) ERROR: Chips do not match!'
                  chip_ratio = 1.01
              end
              '3':begin
                  if chipR ne '7' then $
                    print, '(fill_gap.pro) ERROR: Chips do not match!'
                  chip_ratio = 1.05
              end
              '4':begin
                  if chipR ne '8' then $
                    print, '(fill_gap.pro) ERROR: Chips do not match!'
                  chip_ratio = 0.97
              end
          endcase
; extract the slit number and mask number from the header of the
; spec1d file.
          slitn = string(sxpar(hdrB, 'SLITNO'), format='(I3.3)')
          mask = sxpar(hdrB, 'SLMSKNAM')
; trim the mask name. TBD: do this properly by checking if this is a
; DEEP2-1HS mask.
          if strlen(mask) gt 4 then mask = strmid(mask, 0, 4)
; find and read-in the corresponding calibSlit files.
          blufile = findfile('calibSlit.' + mask + '.' + slitn + 'B*', count=nblu)
          redfile = findfile('calibSlit.' + mask + '.' + slitn + 'R*', count=nred)
          if nred gt 0 and nblu gt 0 then begin
              blu = mrdfits(blufile[0], 1, /silent)
              red = mrdfits(redfile[0], 1, /silent) 
; measure the ratio in the flats and compare to the hard-coded
; values.
              nrows = n_elements(blu.rawflat[0,*])
              nlam = n_elements(blu.rawflat[*,0])
              flat_blu = median(blu.rawflat[nlam-30:nlam-6,5:nrows-6])
              flat_red = median(red.rawflat[5:30,5:nrows-6])
              flat_ratio = flat_red / flat_blu
              if abs(flat_ratio-chip_ratio)/chip_ratio gt 0.2 then $
                print, '(fill_gap.pro) ERROR: Chips ratio off by more than 20%!'
          endif
          ss1.spec = ss1.spec * chip_ratio
; lastly, fill the gap between the red and blue spectra. 
          nlam = n_elements(ss1.lambda)
          dlam = (ss1.lambda[nlam-1] - ss1.lambda[0]) / nlam
          gapsize = round( (ss2.lambda[0] - ss1.lambda[nlam-1]) / dlam )
          if gapsize lt 0 then begin
              print, '(fill_gap.pro) ERROR: wavelength overlap; no gap to fill!'
              return, 0.
          endif
          lambda1 = ss1.lambda[nlam-1] + dlam
          lambda2 = ss2.lambda[0] - dlam
          gap_lam = makearr(gapsize, lambda1, lambda2)
          gap_ivar = fltarr(gapsize); + 1E-20
          
	  blu_lev = median(ss1.spec[nlam-1:nlam-50])	;changed BL 6/2/2008, was derived from a single array value at the end of the blue chip spectrum
          red_lev = median(ss2.spec[0:nlam+50])
          gap_spec = findgen(gapsize)/(gapsize-1) * $
            (red_lev-blu_lev) + blu_lev
	  spec = [ss1.spec, gap_spec, ss2.spec]
          lambda = [ss1.lambda, gap_lam, ss2.lambda]
          ivar = [ss1.ivar / chip_ratio^2, gap_ivar, ss2.ivar]
;          ss1d = {spec:spec, lambda:lambda, ivar:ivar, $
;                  flat_ratio:flat_ratio, chip_ratio:chip_ratio}
          ss1d = make_ss1d(spec=spec, ivar=ivar, lambda=lambda, $
                           ssBd=ss1, ssRd=ss2, gapsz=gapsize) 
          return, ss1d
      endif

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%% /SPECTRA %%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if keyword_set(spectra) then begin
; FIRST CHECK if KEYWORD npoly WAS SET. 
          if keyword_set(npoly) then npoly = npoly[0] $
          else npoly = 1
;;; NOW CHECK THAT THE SPECTRA ARE IN CORRECTLY ORDERED (ss1=blue,
;;; ss2=red).
          if ss1.lambda[0] gt ss2.lambda[0] then begin
;;; if NOT, SWAP 'EM.
              ss1_temp = ss1
              ss1 = ss2
              ss2 = ss1_temp
          endif
;;; DETERMINE THE MEDIAN LEVELS AT THE EDGES OF THE TWO SPECTRA.
;;; FIRST SMOOTH THE FLUX ARRAYS TO REMOVE EFFECTS FROM EMISSION LINES
;;; OR BOGUS SPIKES. 
          if keyword_set(nbins) then nbins = nbins[0] $
          else nbins = 5
          if keyword_set(binsz) then binsz = binsz[0] $
          else binsz = 100
          med1 = FLTARR(nbins)
          med2 = FLTARR(nbins)
          ss1smooth = ivarsmooth(ss1.spec, ss1.ivar, 7)
          ss2smooth = ivarsmooth(ss2.spec, ss2.ivar, 7)
;    ss1smooth = djs_median(ss1.spec, width=7)
;    ss2smooth = djs_median(ss2.spec, width=7)
          for i=0,nbins-1 DO begin
              med1[i] = MEDIAN( ss1smooth[n1-((i+1)*binsz):n1-1-(i*binsz)] )
              med2[i] = MEDIAN( ss2smooth[0+(i*binsz):0+((i+1)*binsz)-1] )
          endfor
          med1 = REVERSE(med1)
          
;;; NOW FIT A POLYNOMIAL OF ORDER npoly TO THESE POINTS.
          xdata = FINDGEN(nbins) * binsz + binsz / 2
          pow = FINDGEN(npoly)
          mcc_polyfit, xdata, med1, pow, a=a1
          mcc_polyfit, xdata, med2, pow, a=a2
;;; AND DETERMINE THE NUMBER OF PIXEL POSITIONS NEEDED TO FILL GAP. 
          nlam = n_elements(ss1.lambda)
          dlam = (ss1.lambda[nlam-1] - ss1.lambda[0]) / nlam
          npix = (ss2.lambda[0] - ss1.lambda[nlam-1]) / dlam / 2
;;; FROM POLYFIT, DETERMINE EXTRAPOLATED LEVELS AT THE MIDDLE OF THE
;;; GAP. 
          xlev1 = (nbins * binsz) + (npix / binsz) ; - (binsz/2)
          xlev2 = (-npix / binsz) ; + (binsz / 2)
          lev1 = 0.
          lev2 = 0.
          for i=0,npoly-1 DO begin
              lev1 = lev1 + a1[i] * xlev1^(pow[i])
              lev2 = lev2 + a2[i] * xlev2^(pow[i])
          endfor
;;; CHECK THAT THE LEVELS ON THE RED AND BLUE SIDE ARE BOTH POSITIVE.
          if lev1 LT 0 then print, '(fill_gap.pro) ERROR: Blue level is negative!'
          if lev2 LT 0 then print, '(fill_gap.pro) ERROR: Red level is negative!'
          newlev = ABS(lev2 - lev1)/2. + MIN([lev1,lev2])
;;; And then DETERMINE THE FACTORS BY WHICH WE MUST MULTIPLY IN ORDER
;;; TO MATCH THE SPECTRA AT THE newlev.
          x1 = newlev / lev1
          x2 = newlev / lev2
;---------------------------
          slitn = STRING(SXPAR(hdrB, 'SLITNO'), forMAT='(I3.3)')
          mask = SXPAR(hdrB, 'SLMSKNAM')
          ;if STRLEN(mask) gt 4 then mask = STRMID(mask, 0, 4)	;edited out for 1604 data
          mask = strcompress(mask, /remove_all) 	;added for 1604 data
	  blufile = 'calibSlit.' + mask + '.' + slitn + 'B*'
          redfile = 'calibSlit.' + mask + '.' + slitn + 'R*'
	  ;CHANGE THIS IF NOT USING ORELSE DATA!!
          blufile = '/Volumes/Data2/orelse/lemaux/deimos/ORELSEmasks/' + mask + '/*/' + blufile 
	  redfile = '/Volumes/Data2/orelse/lemaux/deimos/ORELSEmasks/' + mask + '/*/' + redfile 
	  blu = MRDFITS(blufile, 1, /SILENT)
          red = MRDFITS(redfile, 1, /SILENT) 
          nrows = n_elements(blu.rawflat[0,*])
          nlam = n_elements(blu.rawflat[*,0])
          flat_blu = MEDIAN(blu.rawflat[nlam-30:nlam-6,5:nrows-6])
          flat_red = MEDIAN(red.rawflat[5:30,5:nrows-6])
          flat_ratio = flat_red / flat_blu
;---------------------------

;;; if KEYWORD plot IS SET, then MAKE PLOT.
          if keyword_set(doplot) then begin
              n = FLOOR(nbins*binsz+npix)
              xfit1 = FINDGEN(n) 
              yfit1 = FLTARR(n)
              xfit2 = FINDGEN(n) + (nbins*binsz) + npix
              yfit2 = FLTARR(n)
              for i=0,npoly-1 DO begin
                  yfit1 = yfit1 + a1[i] * xfit1^(pow[i])
                  yfit2 = yfit2 + a2[i] * xfit1^(pow[i])
              endfor
              xfit1 = FINDGEN(n)
              xfit2 = FINDGEN(n) + (nbins*binsz) + npix
              xpts1 = [ (FINDGEN(nbins)*binsz)+(binsz/2), (nbins*binsz)+npix]
              xpts2 = [ (nbins*binsz)+npix, $
                        FINDGEN(nbins)*binsz+(nbins*binsz)+(2*npix)+(binsz/2) ]
              ypts1 = [med1, lev1]
              ypts2 = [lev2, med2]
              ydat1 = ss1smooth[n1-(nbins*binsz):n1-1]
              ydat2 = ss2smooth[0:nbins*binsz-1]
              ydat = [ydat1, ydat2]
              if NOT(keyword_set(debug)) then begin
                  filenm = 'ps/' + STRMID(file1d, 0, 16) + 'match.ps'
                  SET_PLOT, 'PS'
                  DEVICE, FILE=filenm, /LANDSCAPE, /COLOR
              endif
              PLOT, FINDGEN(nbins*binsz), ydat1, THICK=1, $
                XTITLE="PIXELS", /XSTY, XR=[0, (2*nbins*binsz)+(2*CEIL(npix))], $
                /YSTY, YR=[MIN(ydat),MAX(ydat)], XTHICK=2, YTHICK=2, $
                YTITLE="Counts", TITLE=file1d, CHARSIZE=1.2
              OPLOT, FINDGEN(nbins*binsz) + (nbins*binsz) + (2*npix), $
                ydat2, THICK=1
              OPLOT, xfit1, yfit1, THICK=2, COLOR=getcolor('green')
              OPLOT, xfit2, yfit2, THICK=2, COLOR=getcolor('red')
;;;-----------
              if keyword_set(flat_ratio) then $
                OPLOT, xfit2, yfit2*flat_ratio, THICK=3, LINESTY=2, $
                COLOR=getcolor('blue')
;;;-----------
              OPLOT, xpts1, ypts1, PSYM=2, SYMSIZE=2, COLOR=getcolor('green')
              OPLOT, xpts2, ypts2, PSYM=2, SYMSIZE=2, COLOR=getcolor('red')
              LEGend, ['Blue Level: ' + STRCOMPRESS(STRING(lev1), /REMOVE_ALL), $
                       'Red Level: ' + STRCOMPRESS(STRING(lev2), /REMOVE_ALL), $
                       'New Level: ' + STRCOMPRESS(STRING(newlev), /REMOVE_ALL)]
              if NOT(keyword_set(debug)) then begin
                  DEVICE, /CLOSE
                  SET_PLOT, 'X'
              endif
          endif
;;; AND then FILL-IN THE GAP BETWEEN THE BLUE AND RED SPECTRA.
;;; DETERMINE THE LENgtH OF THE BLUE LAMBDA ARRAY.
          nlam = n_elements(ss1.lambda)
;;; DETERMINE THE WAVELENgtH BINSIZE. 
          dlam = (ss1.lambda[nlam-1] - ss1.lambda[0]) / nlam
;;; DETERMINE THE NUMBER OF ARRAY POSITIONS NEEDED TO FILL THE GAP IN
;;; WAVELENgtH SPACE BETWEEN THE SPECTRA. 
          gapsize = ROUND( (ss2.lambda[0] - ss1.lambda[nlam-1]) / dlam )
;;; CATCH POSSIBLE ERROR.
          if gapsize LT 0 then begin
              print, '(fill_gap.pro) ERROR: wavelength overlap; no gap to fill!'
              return, 0.
          endif
          if gapsize EQ 0 then begin
              spec = [ss1.spec*x1, ss2.spec*x2]
              lambda = [ss1.lambda, ss2.lambda]
              ivar = [ss1.ivar / x1^2, ss2.ivar / x2^2]
;;; CREATE THE OUTPUT STRUCTURE.
;              ss1d = {spec:spec, lambda:lambda, ivar:ivar, $
;                      blevel:lev1, rlevel:lev2, midpt:newlev}
              ss1d = make_ss1d(spec=spec, ivar=ivar, lambda=lambda, $
                               ssBd=ss1, ssRd=ss2, gapsz=gapsize) 
;;; AND return THE STRUCTURE.
              return, ss1d
          endif
;;; CREATE A LAMBDA ARRAY TO RUN FROM lambda1 TO lambda2.
          lambda1 = ss1.lambda[nlam-1] + dlam
          lambda2 = ss2.lambda[0] - dlam
          gap_lam = makearr(gapsize, lambda1, lambda2)
;;; CREATE A CORRESPONDING ARRAY TO FILL THE ivar. 
          gap_ivar = fltarr(gapsize) + 1E-20
;;; CREATE AN ARRAY TO FILL THE GAP IN THE SPECTRA. SCALE THE FLUX EQUAL
;;; TO THE PRODUCT OF THE LOCAL MEDIANS.
          gap_spec = fltarr(gapsize) + newlev
;;; CREATE NEW spec, ivar, and lambda ARRAYS.
          spec = [ss1.spec*x1, gap_spec, ss2.spec*x2]
          lambda = [ss1.lambda, gap_lam, ss2.lambda]
          ivar = [ss1.ivar / x1^2, gap_ivar, ss2.ivar / x2^2]
;;; CREATE THE OUTPUT STRUCTURE.
;          ss1d = {spec:spec, lambda:lambda, ivar:ivar, $
;                  blevel:lev1, rlevel:lev2, midpt:newlev}
          ss1d = make_ss1d(spec=spec, ivar=ivar, lambda=lambda, $
                           ssBd=ss1, ssRd=ss2, gapsz=gapsize) 
;;; AND return THE STRUCTURE.
          return, ss1d
      endif

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%% DEFAULT %%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      nlam = n_elements(ss1.lambda)
      dlam = (ss1.lambda[nlam-1] - ss1.lambda[0]) / nlam
      gapsize = ROUND( (ss2.lambda[0] - ss1.lambda[nlam-1]) / dlam )
      if gapsize le 0 then begin
          if not(silent) then $
            print, '(fill_gap.pro) ERROR: wavelength overlap; no gap to fill!'
          return, 0.
      endif
      if gapsize gt 1000. then begin
          if not(silent) then $
            print, '(fill_gap.pro) ERROR: gap too large to fill!'
          return, 0.
      endif
      if keyword_set(nbins) then nbins = nbins[0] $
          else nbins = 5
          if keyword_set(binsz) then binsz = binsz[0] $
          else binsz = 100
      if keyword_set(npoly) then npoly = npoly[0] $
          else npoly = 1     
      lambda1 = ss1.lambda[nlam-1] + dlam
      lambda2 = ss2.lambda[0] - dlam
      gap_lam = makearr(gapsize, lambda1, lambda2)
      gap_ivar = fltarr(gapsize) + 1.0E-30
      
      blu_lev_array = fltarr(nbins)  ;changed BL 6/2/2008 for same reason as change in blu_lev/red_lev above (see /flats section), but this is much more elaborate
      blu_lam_array  = fltarr(nbins)
      red_lev_array = fltarr(nbins)       
      red_lam_array = fltarr(nbins)
      	  for i=0, nbins-1 do begin
          blu_lev_array[i] = median(ss1.spec[nlam-(nbins-i)*binsz-1:nlam-(nbins-1-i)*binsz-1])
	  blu_lam_array[i] = mean(ss1.lambda[nlam-(nbins-i)*binsz-1:nlam-(nbins-1-i)*binsz-1])
          red_lev_array[i] = median(ss2.spec[i*binsz:(i+1)*binsz])
          red_lam_array[i] = mean(ss2.lambda[i*binsz:(i+1)*binsz])
	  endfor 
      
      
      lev_array = [blu_lev_array, red_lev_array]
      lam_array = [blu_lam_array, red_lam_array]
      lev_poly = polyfit(lam_array, lev_array, npoly)
      ;blu_lev = median(ss1.spec[nlam-51:nlam-1]) ;edited out BL 6/2/2008 for same reason as above (see /flats section)
      ;red_lev = median(ss2.spec[0:50])
      flam = 0.   
      
         for i=0, npoly do begin
         flam = flam + lev_poly[i]*gap_lam^(i)
      	 endfor
      
      gap_spec = fltarr(gapsize) + flam ;findgen(gapsize)/(gapsize-1)*(red_lev-blu_lev) + blu_lev
      ;gap_ivar = 1/gap_spec^2				;changed 4/7/09 BL, making the inverse variance associated with gap counts Poissonian, is bullshit, but if I'm going to use the interpolation over the gap for measurements, this is a reasonable estimate of the error. Previous method was to set all ivar counts to 1e-30 (variance = 1e30)
      spec = [ss1.spec, gap_spec, ss2.spec]
      lambda = [ss1.lambda, gap_lam, ss2.lambda]
      ivar = [ss1.ivar, gap_ivar, ss2.ivar]
      ss1d = make_ss1d(spec=spec, ivar=ivar, lambda=lambda, $
                       ssBd=ss1, ssRd=ss2, gapsz=gapsize) 
      return, ss1d
      
 
  endelse
  

end



