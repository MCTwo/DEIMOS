;+
; NAME:
;   pspec
;
; PURPOSE:
;   Routine for plotting spectra from DEIMOS 1D spectro outputs.
;
; CALLING SEQUENCE:
;   pspec, maskname, slitname, znum=znum, nsmooth=nsmooth, $
;    zline=zline, nosyn=nosyn, noerr=noerr, $
;    psfile=psfile, xrange=xrange, yrange=yrange, noerase=noerase, $
;    restframe=restframe, zwarning=zwarning, topdir=topdir, $
;    _EXTRA=KeywordsForSplot
;
; INPUTS:;
;   maskname   - mask name (string or integer)
;
; OPTIONAL INPUTS:
;   slitname   - slit name (string or integer)
;   znum       - If set, then return not the best-fit redshift, but the
;                ZUM-th best-fit; e.g., set ZNUM=2 for second-best fit.
;   nsmooth    - If set, then boxcar smooth both the object and synthetic
;                spectra with a width equal to NSMOOTH.
;   zline      - If set, then overplot the emission line fits.
;   nosyn      - If set, then do not overplot the synthetic fit spectrum.
;   noerr      - If set, then do not overplot the error vector.
;   psfile     - NOT IMPLEMENTED ---
;                If set, then send plot to a PostScript file instead of
;                to the SPLOT interactive widget.  The PostScript file name
;                can be set explicitly, e.g. with PSFILE='test.ps'.  Or if
;                you simply set this as a flag, e.g. with /PSFILE, then the
;                default file name is "spec-pppp-mmmmm-fff.ps",
;                where pppp=plate number, mmmmm=MJD, fff=fiber ID.
;                If FIBERID is specified, then put all plots in a single file
;                named "spec-pppp-mmmmm.ps".
;   restframe  - If set, then plot the wavelengths in the rest frame,
;                e.g. divide the wavelengths by (1+z).
;   zwarning   - If set, then only select those non-sky spectra where the
;                ZWARNING flag has been set; can be used with or without
;                specifying slits with slitname
;   topdir     - Top-level directory for data; default to the environment
;                variable $D2_RESULTS
;   _EXTRA     - Kewords for SPLOT, such as XRANGE, YRANGE, THICK.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The data are read with deimos_readspec.  See the documentation for that
;   routine to see how to set environment variables that describe where
;   the data files are.
;
; EXAMPLES:
;   Plot the spectrum of mask 4206, slit 23 using the SPLOT plotting tool:
;     IDL> pspec, 4206, 23
;
;   The spectrum is shown in white, the errors in red (except masked points
;   are set to zero), and the best-fit eigenspectrum in blue. The mouse
;   buttons will zoom in (left), recenter (center), or zoom out (right).
;   The frame can be saved as a PostScript file by selecting File->WriteEPS
;   from the left-hand corner. 
;
;   Make the same plot, but boxcar-smooth the spectrum and limit the
;   wavelength range to [4000,5000] Angstroms:
;     IDL> pspec, 4206, 23, nsmooth=10, xrange=[5000,6000]
;
;   Loop through all the spectra for mask 4206, interactively:
;     IDL> pspec, 4206
;
; PROCEDURES CALLED:
;   dfpsclose
;   dfpsplot
;   djs_icolor()
;   djs_oplot
;   djs_plot
;   djs_xyouts
;   readspec
;   soplot
;   splot
;   sxyouts
;   textoidl()
;
; INTERNAL SUPPORT ROUTINES:
;   plotspec1
;
; REVISION HISTORY:
;   01-Sep-2000  Written by D. Schlegel, Princeton (as plotspec.pro)
;   09-Oct-2002  Hacked up by D. Finkbeiner for use with DEIMOS
;-
;------------------------------------------------------------------------------
pro plotspec1, maskname, slitname, znum=znum, nsmooth=nsmooth, $
 zline=q_zline, nosyn=nosyn, noerr=noerr, ormaskname=ormaskname, andmaskname=andmaskname, $
 psfile=psfile, xrange=passxr, yrange=passyr, noerase=noerase, $
 restframe=restframe, topdir=topdir, EXTRA=KeywordsForSplot

   cspeed = 2.99792458e5
   textcolor = 'green'
   linecolor = 'magenta'

;   readspec, maskname, slitname, mjd=mjd, znum=znum, flux=objflux, $
;    wave=wave, plug=plug, zans=zans, topdir=topdir, /silent
   deimos_readspec, maskname, slitname, lambda=lambda, spec=spec, ivar=ivar, zans=zans

; return if no data
   if NOT keyword_set(spec) then return
   objflux = spec
   objerr = sqrt(1./(ivar + (ivar eq 0))) < 999
   if (keyword_set(restframe)) then lambda = lambda / (1. + zans.z)
   if (NOT keyword_set(spec)) then begin
      print, maskname, slitname, $
       format='("Spectrum not found for maskname=", A, "slit =", A)'
      return
   endif
;   if (NOT keyword_set(noerr)) then $
;    readspec, maskname, slitname, mjd=mjd, flerr=objerr, topdir=topdir, /silent
   if (NOT keyword_set(nosyn)) then $
     deimos_readspec, maskname, slitname, spec=spec, zans=zans, synspec=synspec

;   if (keyword_set(zans) AND keyword_set(q_zline)) then $
;    readspec, maskname, slitname, mjd=mjd, zline=zline, lineflux=lineflux, $
;     topdir=topdir, /silent

   if (keyword_set(nsmooth)) then begin
      if (nsmooth GT 1) then begin
         objflux = ivarsmooth(objflux, ivar, nsmooth)
         if (keyword_set(synspec)) then $
          synspec = smooth(synspec, nsmooth)
         if (keyword_set(lineflux)) then $
          lineflux = smooth(lineflux, nsmooth)
      endif
   endif


   csize = 1.75
   if (keyword_set(passyr)) then begin
      yrange = passyr
      ymin = yrange[0]
      ymax = yrange[1]
   endif else begin
      if (keyword_set(synspec)) then begin
         yrange = minmax([synspec, objflux]) 
      endif  else $
       yrange = minmax(objflux)

      if (yrange[0] EQ yrange[1]) then yrange = minmax(objflux)
      ymin = (1.3 * yrange[0] - 0.3 * yrange[1]) < 0
      ymax = -0.3 * yrange[0] + 1.3 * yrange[1]
      if (ymax EQ ymin) then ymax = ymin + 1
      yrange = [ymin, ymax]
   endelse
   if (keyword_set(passxr)) then xrange = passxr $
    else xrange = minmax(lambda)
   if (keyword_set(ormaskname) OR keyword_set(andmaskname)) then $
    xrange[1] = 1.15 * xrange[1] - 0.15 * xrange[0]

   title = 'Mask ' + maskname + '   Slit ' + strtrim(string(slitname),2) 
   if (keyword_set(restframe)) then xtitle = 'Rest-Frame Wavelength [Ang]' $
    else xtitle = 'Observed Wavelength [Ang]'
;   ytitle = TeXtoIDL('Flux [10^{-17} erg/s/cm^2/Ang]')
   ytitle = 'counts/hour'
   if (keyword_set(psfile)) then begin
      djs_plot, lambda, objflux, xrange=xrange, yrange=yrange, $
       xtitle=xtitle, ytitle=ytitle, $
       title=title, charsize=csize, _EXTRA=KeywordsForSplot, /xstyle, /ystyle
      if (NOT keyword_set(noerr)) then $
       djs_oplot, lambda, objerr, color='red', _EXTRA=KeywordsForSplot
      if (keyword_set(synspec)) then $
       djs_oplot, lambda, synspec, color='blue', _EXTRA=KeywordsForSplot, lw=2
   endif else begin
      if (NOT keyword_set(noerase)) then $
       splot, lambda, objflux, xrange=xrange, yrange=yrange, $
        xtitle=xtitle, ytitle=ytitle, $
        title=title, charsize=csize, _EXTRA=KeywordsForSplot $
      else $
       soplot, lambda, objflux, _EXTRA=KeywordsForSplot
      if (NOT keyword_set(noerr)) then $
       soplot, lambda, objerr, color='red', _EXTRA=KeywordsForSplot
      if (keyword_set(synspec)) then $
       soplot, lambda, synspec, color='blue', _EXTRA=KeywordsForSplot, lw=2
   endelse

   xpos = 0.9 * !x.window[0] + 0.1 * !x.window[1]
   dypos = 0.05 * (!y.window[0] - !y.window[1])
   ypos = !y.window[1] + 1.5 * dypos
   
   if (keyword_set(zans)) then begin
      zans = zans[0]
      cz = zans.z * cspeed
      if (abs(cz) LT 3000) then $
        zstring = '   cz=' + string(cz,format='(f6.0)') + ' km/s' $
       else $
        zstring = '   z=' + string(zans.z,format='(f8.5)')
;      if (zans.zwarning NE 0) then $
;       zstring = zstring + '  ZWARNING=' + strtrim(string(zans.zwarning),2) + ''
;       zstring = zstring + ' (' $
;        + sdss_flagname('ZWARNING', zans.zwarning, /concat) + ')'
      if (keyword_set(znum)) then $
       zstring = zstring + ' (fit #' + strtrim(string(znum),2) + ')'

      if (keyword_set(psfile)) then $
       xyouts, xpos, ypos, zans.class + ' ' + zans.subclass + zstring, $
        charsize=csize, color=djs_icolor(textcolor), /normal $
      else $
       sxyouts, xpos, ypos, zans.class + ' ' + zans.subclass + zstring, $
        charsize=csize, color=textcolor, /normal

      ypos = ypos + dypos
      
      if (keyword_set(psfile)) then begin 
         xyouts, xpos, ypos, $
           TeXtoIDL('\Chi^2_r =' + strtrim(string(zans.rchi2, format='(f6.2)'),2)), $
           charsize=csize, color=djs_icolor(textcolor), /normal 
      endif else begin 
         sxyouts, xpos, ypos, font=-1, $
           TeXtoIDL(' \chi^2_\nu =' + strtrim(string(zans.rchi2, format='(f8.2)'),2)), $
           charsize=csize, color=textcolor, /normal
         sxyouts, xpos, ypos+dypos, font=-1, $
           TeXtoIDL('\Delta\chi^2 =' + strtrim(string(zans.rchi2diff*zans.dof, $
                                                        format='(f8.2)'),2)), $
           charsize=csize, color=textcolor, /normal
      endelse 
   endif

   print, zstring

return

   if (keyword_set(lineflux)) then begin
      if (keyword_set(psfile)) then $
       djs_oplot, lambda, lineflux, color=linecolor, _EXTRA=KeywordsForSplot, lw=2 $
      else $
       soplot, lambda, lineflux, color=linecolor, _EXTRA=KeywordsForSplot, lw=2

      linelambda = zline.linelambda $
       * (1 + zline.linez * (keyword_set(restframe) EQ 0))
      ; Convert line sigma from km/sec to Angstroms
      linesigma = linelambda * zline.linesigma / cspeed
      linepeak = zline.linecontlevel + zline.linearea / (sqrt(2*!pi) * linesigma)
      for iline=0, n_elements(zline)-1 do begin
         if (zline[iline].linearea_err GT 0) then begin
            if (keyword_set(psfile)) then $
             xyouts, linelambda[iline], linepeak[iline], $
              '  '+zline[iline].linename, orient=90, $
              charsize=0.75*csize, color=djs_icolor(linecolor) $
            else $
             sxyouts, linelambda[iline], linepeak[iline], $
              '  '+zline[iline].linename, orient=90, $
              charsize=0.75*csize, color=linecolor
         endif
      endfor
   endif

   if (keyword_set(ormaskname)) then begin
      plotspec_maskname, lambda, ormaskname, psfile=psfile, $
       psym=1, symsize=0.6, color=orcolor, nolabel=keyword_set(andmaskname)
   endif

   if (keyword_set(andmaskname)) then begin
      plotspec_maskname, lambda, andmaskname, psfile=psfile, $
       psym=6, symsize=0.6, color=andcolor
   endif

   if (keyword_set(primtarget)) then begin
      ypos = ypos + dypos
      if (keyword_set(psfile)) then $
       xyouts, xpos, ypos, 'PRIMTARGET = ' + primtarget, $
        charsize=csize, color=djs_icolor(textcolor), /normal $
      else $
       sxyouts, xpos, ypos, 'PRIMTARGET = ' + primtarget, $
        charsize=csize, color=textcolor, /normal
   endif

   if (keyword_set(sectarget)) then begin
      ypos = ypos + dypos
      if (keyword_set(psfile)) then $
       xyouts, xpos, ypos, 'SECTARGET = ' + sectarget, $
        charsize=csize, color=djs_icolor(textcolor), /normal $
      else $
       sxyouts, xpos, ypos, 'SECTARGET = ' + sectarget, $
        charsize=csize, color=textcolor, /normal
   endif

   return
end



;------------------------------------------------------------------------------



pro pspec, maskname, slitname, znum=znum, nsmooth=nsmooth, $
 zline=zline, nosyn=nosyn, noerr=noerr, $
 psfile=psfile, xrange=xrange, yrange=yrange, noerase=noerase, $
 restframe=restframe, zwarning=zwarning, topdir=topdir, $
 _EXTRA=KeywordsForSplot
  
  if (n_params() LT 1) then begin
     doc_library, 'pspec'
     return
  endif
  
  quiet = !quiet
  !quiet = 1
  
  nmaskname = n_elements(maskname) 
  if nmaskname eq 0 then message, 'You must pass maskname'
    
  nslit = n_elements(slitname) 
  
  if nslit le 1 then begin 
     nslit = 150   ; ???
     masklist = replicate(maskname, nslit)
     slitlist = indgen(nslit)
  endif 
  if size(masklist, /tname) NE 'STRING' then $
    masklist=string(masklist, format='(I4.4)')
      
     ;----------
     ; Loop over each plot
  
  islit = 0
  while (islit LT nslit) do begin
     
      ;----------
      ; Open the PostScript file if appropriate

;     if (keyword_set(psfile)) then begin
;        if (NOT keyword_set(q_onefile)) then begin
;           if (size(psfile,/tname) EQ 'STRING') then $
;             psfilename = psfile $
;           else $
;             psfilename = string(masknamelist[islit], $
;              slitname[islit], format='("spec-",i4.4,"-",i5.5,"-",i3.3,".ps")')
;        endif
;        
;        if (NOT keyword_set(q_onefile) OR islit EQ 0) then begin
;           dfpsplot, psfilename, /color, /square
;        endif
;     endif
     
; convert to string
     if size(slitlist, /tname) NE 'STRING' then $
       slitlist=string(slitlist, format='(I3.3)')
     plotspec1, masklist[islit], slitlist[islit], $
       znum=znum, nsmooth=nsmooth, zline=zline, nosyn=nosyn, noerr=noerr, $
       psfile=psfile, $
       xrange=xrange, yrange=yrange, noerase=noerase, $
       restframe=restframe, topdir=topdir, _EXTRA=KeywordsForSplot
     
     if (islit LT nslit-1) then begin
        if (keyword_set(nsmooth)) then $ 
          sstring = ' (currently=' + strtrim(string(nsmooth),2) + ')' $
        else $
          sstring = ''
        
        print, 'Press b=back one slit'
        print, '      m=select new maskname'
        print, '      n=select new slitname'
        print, '      q=quit (and enter interactive mode for this plot)'
        print, '      s=change boxcar smoothing' + sstring
        print, '      x=change X plotting range'
        print, '      y=change Y plotting range'
        print, '      z=change which PCA-fit to overplot'
        print, '      any other key=forward'
        
        cc = strupcase(get_kbrd(1))
        print,cc
        case cc of
           'B': islit = (islit - 1) > 0
           'M': begin
              maskname = ''  ; force type string
              read, maskname, prompt='Enter new maskname: '
              nslit = 150
              help, maskname
              masklist = replicate(maskname, nslit)
              slitlist = indgen(nslit)
           end
           'N': begin
              read, newslit, prompt='Enter new slitname: '
              islit = ((long(newslit)) > 0) < (nslit-1)
           end
           'Q': islit = nslit
           'S': begin
              read, nsmooth, prompt='Enter boxcar smoothing width (0=none): '
              nsmooth = (long(nsmooth) > 0) < 99
           end
           'X': begin
              read, xmin, xmax, prompt='Enter new X range values (0 0=full range): '
              if (xmin EQ 0 AND xmax EQ 0) then xrange = 0 $
              else xrange = [xmin, xmax]
           end
           'Y': begin
              read, ymin, ymax, prompt='Enter new Y range values (0 0=full range): '
              if (ymin EQ 0 AND ymax EQ 0) then yrange = 0 $
              else yrange = [ymin, ymax]
           end
           'Z': begin
              read, znum, prompt='Enter 1=best redshift, 2=2nd best, ...: '
              znum = long(znum) > 0
           end
           else: islit = islit + 1
        endcase
     endif else begin
        islit = nslit
     endelse

  endwhile

  !quiet = quiet
  return
end
