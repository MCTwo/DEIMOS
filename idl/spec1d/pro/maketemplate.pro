;+
; NAME: maketemplate
;
; PURPOSE:
;   Make an artificial spectral template for any spectral line list,
;   in order to test the redshift fitting process.
;
; CALLING SEQUENCE: 
; maketemplate, linelistfile,templatefile [,minlamb=minlamb, maxlamb=maxlamb, $
;               logresol=logresol]
; 
; INPUTS:
;   linelistfile -- the data file containing the line list. 
;                   Format should be like
;                   ~renbin/deep/DEEPlinelist.dat
;                   The data file should have four columns:
;                 name -- strings without space in it
;                 lambda -- wavelength,floating point
;                 width -- FWHM of the line
;                 fvalue -- peak flux of the line
;
; KEYWORDS:
;   minlambda -- minimal wavelength,in Angstrom 
;   maxlambda -- maximal wavelength,in Angstrom
;   logresol -- required resolution in log(lambda) space for the
;               output template,units in log10(Angstrom)
;   vdispersion -- internal velocity dispersion to assume for template
;                  (0 default)
;   template -- template created
;   header -- FITS header of output file
;  
; OUTPUTS:
;   templatefile -- filename of the artificial spectral template
; 
; OPTIONAL OUTPUTS:
;
; RESTRICTIONS:
;
; EXAMPLES:
;   maketemplate,'DEEPlinelist.dat','faketemplate.fits', $
;                minlambda=2000,maxlambda=7000,logresol=2.e-5
; COMMENTS:
;   it is good to keep the minlambda and maxlamba to have substantial
;   buffer beyond last lines, so that dofarr does not vary overmuch
;   during the reduce1d process.
;
; REVISION HISTORY:
;       Wed Jul 24 14:46:07 2002, Renbin Yan
;           Sep 27,2002 MD -- added internal dispersion keyword
;----------------------------------------------------------------------

function gauss_inte, wave, center, fwhm, height
  sig = fwhm/(2.*sqrt(alog(4.)))
  return, height*exp(-1./2.*(wave-center)^2/sig^2)
end

pro maketemplate, linelistfile, templatefile, minlambda=minlambda,$
     maxlambda=maxlambda, logresol=logresol,  vdispersion=vdispersion, $
     template=result,  header=hdr

  readcol, linelistfile, name, lambda, width, fvalue, format='A F F F'
  airtovac, lambda ;convert wavelengths to vacuum

  vdisp = 0.
  if (NOT keyword_set(minlambda)) then minlambda = 3000.
  if (NOT keyword_set(maxlambda)) then maxlambda = 7500.
  if (NOT keyword_set(logresol)) then logresol = 2.e-5
  if (keyword_set(vdispersion)) then vdisp = vdispersion

;factor to add in quadrature to include internal FWHM
  voverc = vdisp/299000.*(2.*sqrt(alog(4.))) 
  

  ind = where(lambda gt minlambda and lambda lt maxlambda)
  name = name[ind]
  lambda = lambda[ind]
  width = width[ind]
  fvalue = fvalue[ind]



  num = ceil( (alog10(maxlambda)-alog10(minlambda))/logresol )
  result = fltarr(num)

  nline = n_elements(name)


  for kk=0, nline-1 do begin
     center = lambda[kk]
     fwhm = sqrt(width[kk]^2 + (center*voverc)^2) ;net FWHM
     height = fvalue[kk]
     bottom = alog10(center-2.5*width[kk]) > alog10(minlambda)
     top = alog10(center+2.5*width[kk]) <  alog10(maxlambda)
     bottom = floor((bottom-alog10(minlambda))/logresol)
     top = ceil((top-alog10(minlambda))/logresol)
     logwavearrary = (findgen(top-bottom+1)+bottom)*logresol+alog10(minlambda)
     result[bottom:top] = result[bottom:top] + gauss_inte(10^logwavearrary, $
          center, fwhm, height)
  endfor

;Write output template file
  mkhdr, hdr, result, /extend
  sxaddpar, hdr, 'OBJECT', 'GALAXY  ', before='EXTEND'
  sxaddpar, hdr, 'COEFF0', alog10(minlambda), after='OBJECT'
  sxaddpar, hdr, 'COEFF1', logresol, after='COEFF0'
  writefits, templatefile,result, hdr

end




