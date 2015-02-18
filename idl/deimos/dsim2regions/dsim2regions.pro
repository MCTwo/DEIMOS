;+
; NAME:
;	DSIM2REGIONS
;
; PURPOSE:
;	This procedure will generate a DS9 regions file which can be
;	overlain atop a WCS-enhanced image to show the DEIMOS FOV,
;	DEIMOS TV guider FOV, locations of selected targets and
;	locations of slits.
;
; CATEGORY:
;	Observing planning.
;
; CALLING SEQUENCE:
;       dsim2regions, infile, outfile
;
; INPUTS:
;	infile:	FITS file produced by DSIMULATOR
;
; KEYWORD PARAMETERS:
;       SLITS: if set, show slits as red polygons and alignment boxes and blue
;       boxes
;
;       TARGETS: if set, show targets as green diamonds, whether or not slits
;       were assigned.  Also show guider alignment stars as yellow diamonds.
;       
;       TV: if set, show TV fov as yellow polygon
;
;       MASK: if set, show accessible mask area as red polygon
;
;       LABEL: if set, then label slits with slit number and targets
;       with target name
;
; OUTPUTS:
;       This procedure generates a ds9 regions file which can be read
;       into ds9 and plotted atop an image of the corresponding field.
;       http://hea-www.harvard.edu/RD/ds9/ref/region.html
;
; REQUIRED ROUTINES:
;       - Requires routines from the IDL astronomy library at:
;       http://idlastro.gsfc.nasa.gov/homepage.html
;
; REQUIRED FILES:
;       - deimos_mask.dat: file containing DEIMOS mask FOV
;       - deimos_tv.dat: file containing DEIMOS TV FOV
;
; EXAMPLE:
;       1) Generate a ds9 regions file 'foo.regions' for mask
;       'foo.fits', marking slits, targets, DEIMOS FOV and TV FOV, but
;       do not labels slits:
;               dsim2regions, 'foo.fits', 'foo.regions'
;
;       2) Same as (1), but label slits and targets:
;               dsim2regions, 'foo.fits', 'foo.regions', /LABEL
;
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2006-Sep-16	GDW	Original version
; 	2008-Apr-06	GDW	Minor mods plus doco
;-

;-----------------------------------------------------------------------
pro unflatten, dx, dy, ra_center, dec_center, ra, dec
;-----------------------------------------------------------------------
;+
; NAME:
;	UNFLATTEN
;
; PURPOSE:
;       This routine converts a set of arcsecond offets from the field
;       center into RA and Dec.
;
; CATEGORY:
;       Astrometry
;
; CALLING SEQUENCE:
;       UNFLATTEN, dx, dy, ra_center, dec_center, ra, dec
;
; INPUTS:
;       dx:             RA offset from field center (arcsec)
;       dy:             Dec offset from field center (arcsec)
;	ra_center:      RA of the field center [deg]
;       dec_center:     Dec of the field center [deg]
;	
; OUTPUTS:
;	ra:             right ascension [deg]
;       dec:            declination [deg]
;
; EXAMPLE:
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2006-Sep-02	GDW	Original version
;-
;-----------------------------------------------------------------------

secperdeg = 3600.

sindeccent = sin(dec_center*!dtor)
cosdeccent = cos(dec_center*!dtor)

;; convert offset from arcsec to deg...
xflat = dx/secperdeg
yflat = dy/secperdeg

;; convert offset to radians...
x = xflat/!radeg
y = yflat/!radeg                   

d = atan(sqrt(x^2 + y^2))
theta = atan(-x,y)
dec = asin(cos(d)*sindeccent + sin(d)*cosdeccent*cos(theta))
delta_ra = asin(sin(theta)*sin(d)/cos(dec))*!radeg
ra = ra_center + delta_ra
dec = dec*!radeg

end


;-----------------------------------------------------------------------
pro plot_outline, infile, unit, outline, type, color
;-----------------------------------------------------------------------
; Purpose:
;       Helper routine to plot data from file.
;-----------------------------------------------------------------------

    ;; get mask outline (arcsec offsets from field center at PA=0)...
    readcol, outline, format='(f,f)', x, y
    x = -x - 8.
    y = y + 22.

    ;; read input...
    table = mrdfits(infile,3,header)

    ;; get mask PA...
    ra0  = table.ra_pnt
    dec0 = table.dec_pnt
    pa   = table.pa_pnt - 90.

    ;; rotate outline...
    u=x*cos(pa*!dtor)-y*sin(pa*!dtor)
    v=x*sin(pa*!dtor)+y*cos(pa*!dtor)

    ;; convert to radec...
    unflatten, u, v, ra0, dec0, ra, dec

    ;; plot mask as polygon...
    printf, unit, '# BEGIN ', type
    n = n_elements(ra)
    printf, unit, format='(a,$)', 'polygon'
    for i=0,n-1 do begin
        printf, unit, format='(f12.6, "d ", f12.6, "d ", $)', ra[i], dec[i]
    endfor
    printf, unit, ' # color = ', color

    printf, unit, '# END ', type
end

;-----------------------------------------------------------------------
pro dsim2regions, infile, outfile, SLITS=slits, TARGETS=targets, $
                  TV=tv, MASK=mask, LABEL=label
;-----------------------------------------------------------------------

;; check keywords...
do_slits   = keyword_set(slits)
do_targets = keyword_set(targets)
do_tv      = keyword_set(tv)
do_mask    = keyword_set(mask)
do_label   = keyword_set(label)

;; default is to show everything...
if not do_slits and not do_targets and not do_tv and not do_mask then begin
    do_slits = 1
    do_targets = 1
    do_tv = 1
    do_mask = 1
endif

;; definitions...
quote = '"'
maskfile = 'deimos_mask.dat'
tvfile = 'deimos_tv.dat'

;; assign colors...
color_align  = "magenta"
color_slit   = "green"
color_tv     = "red"
color_mask   = "white"
color_target = "yellow"
color_guide  = "blue"

;; open output file...
openw, unit, outfile, /get_lun
printf, unit, '# ', outfile
printf, unit, '# ds9 regions file generated from ', infile
printf, unit, '# by dsim2regions.pro'
printf, unit, ''

;; process slits...
if do_slits then begin

    printf, unit, '# BEGIN SLITS'

    ;; read input...
    table = mrdfits(infile,4,header)

    ;; count slits...
    n = n_elements(table.slitname)

    ;; loop over slits...
    printf, unit, 'global color=red'
    printf, unit, 'j2000'

    ;; generate regions...
    for i=0,n-1 do begin
        ra     = (table.slitra)[i]
        dec    = (table.slitdec)[i]
        length = (table.slitlen)[i]
        width  = (table.slitwid)[i]
        pa     = (table.slitlpa)[i]+90.
        id     = (table.slitname)[i]
        type   = (table.slittyp)[i]

        ;; assign special color to alignment boxes...
        if strcmp(type,'a',/fold) eq 1 then begin
            color = color_align
        endif else begin
            color = color_slit
        endelse

        ;; output slit...
        printf, unit, format='("box", $)'
        printf, unit, format='(" ", f12.6, "d", $)', ra
        printf, unit, format='(" ", f12.6, "d", $)', dec
        printf, unit, format='(" ", f5.1, a, $)', length, quote
        printf, unit, format='(" ", f5.1, a, $)', width, quote
        printf, unit, format='(" ", i4, "d", $)', nint(pa)
        printf, unit, format='(" # color = ", a, $)', color
        if do_label then $
          printf, unit, format='(", text = {", a, "}", $)', id
        printf, unit, ''
    endfor

    printf, unit, '# END SLITS'

endif

;; process targets...
if do_targets then begin

    printf, unit, '# BEGIN TARGETS'

    ;; read input...
    table = mrdfits(infile,1,header)

    ;; count slits...
    n = n_elements(table.object)

    ;; generate regions...
    for i=0,n-1 do begin
        
        id      = (table.object)[i]
        ra      = (table.ra_obj)[i]
        dec     = (table.dec_obj)[i]
        equinox = (table.equinox)[i]
        type    = (table.objclass)[i]

        ;; assign special color to alignment boxes...
        if strcmp(type,'a',1,/fold) eq 1 then begin
            color = color_target
        endif else if strcmp(type,'g',1,/fold) eq 1 then begin
            color = color_guide
        endif else begin
            color = color_target
        endelse

        ;; output slit...
        printf, unit, format='("diamond point", $)'
        printf, unit, format='(" ", f12.6, "d", $)', ra
        printf, unit, format='(" ", f12.6, "d", $)', dec
        printf, unit, format='(" # color = ", a, $)', color
        if do_label then $
          printf, unit, format='(", text = {", a, "}", $)', id
        printf, unit, ''
    endfor

    printf, unit, '# END TARGETS'

endif

;; process mask...
if do_mask then plot_outline, infile, unit, maskfile, 'MASK', color_mask

;; process tv...
if do_tv then plot_outline, infile, unit, tvfile, 'TV', color_tv

;; close file...
free_lun, unit

end
