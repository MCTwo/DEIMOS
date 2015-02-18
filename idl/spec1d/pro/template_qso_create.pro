;+
; NAME:
;    template_qso_create
;
; PURPOSE:
;    combines SDSS template with emission line template at higher resolution
; 
; CALLING SEQUENCE:
;    template_qso_create,filename
;
; INPUTS:
;    filename  -- name of output file
;
; COMMENTS:
;    this procedure uses the high SN SDSS template of absorption, with
;    very wide spectral coverage, and adds to it a high resolution
;    emission line template.  The two templates together should be a
;    good template for DEEP2 analysis.
;
; REVISION HISTORY:
;   MD 6jun03 -modified from template_create
;
;----------------------------------------------------------------------
pro template_qso_create, fileout

  eigendir = concat_dir( GETENV('IDLSPEC1D_DIR'), 'templates' )
  eigenfile = 'spEigenQSO.fits'
  sdss = readfits(eigendir +'/' + eigenfile, hdr)

  c0 = sxpar(hdr, 'COEFF0')
  c1 = sxpar(hdr, 'COEFF1')

;  sdss_select = sdss[*, 0]+ sdss[*, 2] ;sum of 2 eigenspectra
;  sdss_select = sdss_select/mean(sdss_select) ;normalize
  sdss_select = sdss[*, 0]
  nspec = n_elements(sdss_select)

;  sdss_uv = eigendir+ '/lrg.0701a.05.faruv.spave1' ;file of aver. spectra
;  readfast, sdss_uv, uv_in, skip=3, ncol=2 ;read this data file
;  uvlam = reform(uv_in[0, *])
;  uvspec = smooth(reform(uv_in[1, *]), 3) ;deconstruct and smooth a bit
  

  logl = findgen(nspec)*c1 +c0

; now splice the two SDSS spectra together. At 3600AA they
; close to equal
;  minold = max(where(10.^logl lt 3600.))
;  maxnew = max(where(uvlam le 3600.))
;  minnew = min(where(uvlam ge 2700.)) ; end of believable spectrum
;  start = minold- (maxnew-minnew) ;start of composite array with data

;  sdss_composite = sdss_select*0. ;set new array
;  sdss_composite[minold:nspec-1] = 2.*sdss_select[minold:nspec-1] ;visible
;  sdss_composite[start:minold-1] = uvspec[minnew:maxnew-1]
;  sdss_composite[0:start-1] = sdss_composite[start] ;extend to 0 lambda
  offset = 3000
  sdss_composite = sdss_select[offset:*]

  dloglambda = 2.e-5
;  inputfile = concat_dir(GETENV('IDLSPEC1D_DIR'), 'etc')+ '/DEEPlinelist.dat' 
  c0 = c0-2.5*dloglambda  +offset*c1 ;adjust c0 for expansion, offset

  lambda0 = 10.^c0 ;initial wavelength, including expansion of scale
; and presumption that SDSS spectra are pointing at central wavelengths  
  lambda1 = 10.^logl[nspec-1]  ;final wavelength
; create emission line template
;  maketemplate,  inputfile, 'junk.fits', minl=lambda0, maxl=lambda1, $
;         logres=dloglambda, vdisp=60., template=emitemplate, header=headf

;  biglogl = c0 + findgen(n_elements(emitemplate))*dloglambda
  scale = c1/dloglambda

  nspec = long(n_elements(sdss_composite)*fix(scale))

  biglogl = c0 +findgen(nspec)*dloglambda
 
  logl = logl[offset:*]
  sdss_c = sdss_composite

  sspline = spl_init(logl, sdss_c, yp0=0., ypn_1=0.) 
;set spline to have constant extrapolation outside of bounds
  sdss_refined = spl_interp(logl, sdss_c, sspline, biglogl)
  sdss_refined = sdss_refined/mean(sdss_refined) ;normalize

  result = float(sdss_refined)

  fxaddpar, hdr, 'COEFF0', c0
  fxaddpar, hdr, 'COEFF1', dloglambda
  writefits, fileout, result, hdr

end


