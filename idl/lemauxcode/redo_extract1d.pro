
pro rename_files, maskname

maskname = strcompress(maskname, /rem)

slit = findfile('slit.*.fits.gz', count=nslit)
spslit = findfile('spSlit.*.fits.gz', count=nspslit)
spec1d = findfile('spec1d.*.fits', count=nspec1d)

for ii=0,nslit-1 do begin
    parts = strsplit(slit[ii], '.', /extract)
    spstr = 'mv ' + slit[ii] + ' ' + $
      'slit.' + maskname + '.' + parts[2] + '.fits.gz'
    spawn, spstr
endfor

for ii=0,nspslit-1 do begin
    parts = strsplit(spslit[ii], '.', /extract)
    spstr = 'mv ' + spslit[ii] + ' ' + $
      'spSlit.' + maskname + '.' + parts[2] + '.fits.gz'
    spawn, spstr
endfor

for ii=0,nspec1d-1 do begin
    parts = strsplit(spec1d[ii], '.', /extract)
    spstr = 'mv ' + spec1d[ii] + ' ' + $
      'spec1d.' + maskname + '.' + parts[2] + '.' + parts[3] + '.fits'
    spawn, spstr
endfor


end

pro redo_extract1d, oldspec1dfile, slitfiles, objpos, fwhm

bfile = slitfiles[0]
rfile = slitfiles[1]

posB = objpos[0]
posR = objpos[1]

fwhmB = fwhm
fwhmR = fwhm
avgfwhmB = fwhm
avgfwhmR = fwhm

nsig_box = 1.1
nsig_opt = 1.75

blu_opt = extract1d(bfile, posB, fwhmB, /horne, $
			   nsigma=nsig_opt)
blu_box = extract1d(bfile, posB, avgfwhmB, $
			   /boxsprof, nsigma=nsig_box) 
red_opt = extract1d(rfile, posR, fwhmR, /horne, $
			   nsigma=nsig_opt)
red_box = extract1d(rfile, posR, avgfwhmR, $
			   /boxsprof, nsigma=nsig_box)

specfile = 'new' + oldspec1dfile

hdr = headfits(oldspec1dfile, ext=1, /silent)
hdrB = copy_header(hdr, 'Bxspf-R')
sxaddpar, hdrB, 'objpos', posB, $
'Object Position used in extraction', after='objno'
sxaddpar, hdrB, 'ext_fwhm', avgfwhmB, $
'FWHM employed in object extraction', after='cor_fwhm'

mwrfits, blu_box, specfile, hdrB, /silent, /create

hdrR = copy_header(hdr, 'Bxspf-R')
sxaddpar, hdrR, 'objpos', posR, $
'Object Position used in extraction', after='objno'
sxaddpar, hdrR, 'ext_fwhm', avgfwhmR, $
'FWHM employed in object extraction', after='cor_fwhm'

mwrfits, red_box, specfile, hdrR, /silent


sxaddpar, hdrB, 'EXTNAME', 'Horne-B', 'Extension Name'
sxaddpar, hdrR, 'EXTNAME', 'Horne-R', 'Extension Name' 

mwrfits, blu_opt, specfile, hdrB, /silent 
mwrfits, red_opt, specfile, hdrR, /silent

end
