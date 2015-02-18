
; file = spec1d FITS file.
; zresult = zresult structure for an object.
; nlsky = use non-local sky extraction.

pro vdisp_wrapper, file, zresult=zresult, nfind=nfind, $
                   airmass=airmass, nlsky=nlsky

; check if the nlsky keyword was set.
  if keyword_set(nlsky) then nlsky=nlsky[0] else nlsky = 0

; grab the boxcar extraction from the spec1d file.
  if nlsky then ss1d = fill_gap(file, /boxsprof, /nlsky,/tweak) $
  else ss1d = fill_gap(file, /boxsprof,/tweak)

; check that the fill_gap routine succeeded in filling the gap.
  if size(ss1d, /tname) ne 'STRUCT' then begin
      print, 'Skipping slit ' + strcompress(zresult[0].slitname, /rem) + $
        '...no vdisp measured!'
      zresult.vdisp = 999.0
      zresult.vdisp_err = 999.0
  endif else begin
; correct for telluric absorption bands (do this before shifting to
; the vacuum wavelengths!).
      remove_telluric, ss1d, airmass
      
; convert lambda values from air to vacuum.
      lambda = ss1d.lambda 
      airtovac, lambda
; replace the lambda values and convert from double to float.
      ss1d.lambda = float(lambda)

; set some parameters for the zrefind.
      pspace = 1
      width = 3 * pspace
      nfind0 = 1
      npoly = 0
      columns = [0,1,2]
; call the vdispfit routine interatively for each z value...refind the
; z value since the boxcar spectrum will likely vary from the optimal.
      nz = n_elements(zresult)
      for ii=0,nfind-1 do begin
          zmin = zresult[ii].z - 0.03
          zmax = zresult[ii].z + 0.03
          ztune = zfind(ss1d, zmin=zmin, zmax=zmax, /linear_lambda, $
                        eigenfile=eigenfile_gal, eigendir=eigendir, $
                        npoly=npoly, pspace=pspace, width=width, $
                        columns=columns, nfind=nfind0)
          if finite(ztune.z[0]) then begin
              print, 'Difference between optimal-boxcar redshifts: ' + $
                strcompress(zresult[ii].z - ztune.z, /rem)
              zobj = ztune.z
          endif else begin
              print, 'Unable to determine z for boxcar extraction!'
              zobj = zresult[ii].z
          endelse
          vdispfit, ss1d, npoly=2, zobj=zobj, sigma=sigma, $
            sigerr=sigerr
          zresult[ii].vdisp = sigma
          zresult[ii].vdisp_err = sigerr
      endfor
  endelse


end

