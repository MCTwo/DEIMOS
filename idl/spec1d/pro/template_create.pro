;+
; NAME:
;    template_create
;
; PURPOSE:
;    combines SDSS template with emission line template at higher resolution
; 
; CALLING SEQUENCE:
;    template_create,filename
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
;   MD 29sep02
;   MD 16dec02, corrected start point by interpolation error, added
;   new SDSS spectrum in UV (see Eisenstein etal astroph/0212087)
;----------------------------------------------------------------------
pro template_create, fileout

  eigendir = concat_dir( GETENV('IDLSPEC1D_DIR'), 'templates' )
  eigenfile = 'spEigenGal.fits'
  sdss = readfits(eigendir +'/' + eigenfile, hdr)

  c0 = sxpar(hdr, 'COEFF0')
  c1 = sxpar(hdr, 'COEFF1')

  sdss_select = sdss[*, 0]+ sdss[*, 2] ;sum of 2 eigenspectra
  sdss_select = sdss_select/mean(sdss_select) ;normalize
  nspec = n_elements(sdss_select)

  sdss_uv = eigendir+ '/lrg.0701a.05.faruv.spave1' ;file of aver. spectra
  readfast, sdss_uv, uv_in, skip=3, ncol=2 ;read this data file
  uvlam = reform(uv_in[0, *])
  uvspec = smooth(reform(uv_in[1, *]), 3) ;deconstruct and smooth a bit
  

  logl = findgen(nspec)*c1 +c0

; now splice the two SDSS spectra together. At 3600AA they
; close to equal
  minold = max(where(10.^logl lt 3600.))
  maxnew = max(where(uvlam le 3600.))
  minnew = min(where(uvlam ge 2700.)) ; end of believable spectrum
  start = minold- (maxnew-minnew) ;start of composite array with data

  sdss_composite = sdss_select*0. ;set new array
  sdss_composite[minold:nspec-1] = 2.*sdss_select[minold:nspec-1] ;visible
  sdss_composite[start:minold-1] = uvspec[minnew:maxnew-1]
  sdss_composite[0:start-1] = sdss_composite[start] ;extend to 0 lambda


  dloglambda = 2.e-5
  inputfile = concat_dir(GETENV('IDLSPEC1D_DIR'), 'etc')+ '/DEEPlinelist.dat' 

  c0 = c0-2.5*dloglambda ;adjust c0 for expansion

  lambda0 = 10.^c0 ;initial wavelength, including expansion of scale
; and presumption that SDSS spectra are pointing at central wavelengths  
  lambda1 = 10.^logl[nspec-1]  ;final wavelength
; create emission line template
  maketemplate,  inputfile, 'junk.fits', minl=lambda0, maxl=lambda1, $
         logres=dloglambda, vdisp=60., template=emitemplate, header=headf

  biglogl = c0 + findgen(n_elements(emitemplate))*dloglambda

  sspline = spl_init(logl, sdss_composite, yp0=0., ypn_1=0.) 
;set spline to have constant extrapolation outside of bounds
  sdss_refined = spl_interp(logl, sdss_composite, sspline, biglogl)
  sdss_refined = sdss_refined * 10^biglogl ;Change units from ergs/A to photon-counts/A. 
  sdss_refined = sdss_refined/mean(sdss_refined) ;normalize

  ka = mrdfits(eigendir+'/spDeepKAv60.fits',0)
  astar = ka[*,1]/mean(ka[*,1])

  result = [[sdss_refined], [emitemplate], [astar]]
  sxaddpar,headf,'UNITS','photon-counts',after='COEFF1'

  writefits, fileout, result, headf

end


