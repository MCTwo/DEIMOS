pro maketargetable
	
	a = mrdfits('target.SC1NM1.fits',1)
	b = mrdfits('target.SC1NM2.fits',1)
	c = mrdfits('target.SC2NM1.fits',1)
	d = mrdfits('target.SC2NM2.fits',1)
	e = mrdfits('target.CE1.fits',1)
	f = mrdfits('target.FG1.fits',1)
	g = mrdfits('target.FG2.fits',1)
	h = mrdfits('target.GHF1.fits',1)
	i = mrdfits('target.GHF2.fits',1)
	j = mrdfits('target.16XR1.fits',1)
	k = mrdfits('target.16XR2.fits',1)
	l = mrdfits('target.16XR3.fits',1)

	; get rid of stars	
	aa = a[where(a.slitlen ne 4.)]
	bb = b[where(b.slitlen ne 4.)]
	cc = c[where(c.slitlen ne 4.)]
	dd = d[where(d.slitlen ne 4.)]
	ee = e[where(e.slitlen ne 4.)]
	ff = f[where(f.slitlen ne 4.)]
	gg = g[where(g.slitlen ne 4.)]
	hh = h[where(h.slitlen ne 4.)]
	ii = i[where(i.slitlen ne 4.)]
	jj = j[where(j.slitlen ne 4.)]
	kk = k[where(k.slitlen ne 4.)]
	ll = l[where(l.slitlen ne 4.)]

	; start printing
	openw, lun, 'sc1604.slitinfo.2008mar20.dat', width=200, /get_lun
	
	for i=0, n_elements(aa)-1 do begin
	printf, lun, 'SC1NM1    ', aa[i].objname, '   ', aa[i].slit, '   ', aa[i].slitlen, '   ', aa[i].PA, '   ', aa[i].RA, '   ', aa[i].Dec, '   ', aa[i].imag
	endfor

	for i=0, n_elements(bb)-1 do begin
        printf, lun, 'SC1NM2    ', bb[i].objname, '   ', bb[i].slit, '   ', bb[i].slitlen, '   ', bb[i].PA, '   ', bb[i].RA, '   ', bb[i].Dec, '   ', bb[i].imag
	endfor

	for i=0, n_elements(cc)-1 do begin
        printf, lun, 'SC2NM1    ', cc[i].objname, '   ', cc[i].slit, '   ', cc[i].slitlen, '   ', cc[i].PA, '   ', cc[i].RA, '   ', cc[i].Dec, '   ', cc[i].imag
	endfor

	for i=0, n_elements(dd)-1 do begin
        printf, lun, 'SC2NM2    ', dd[i].objname, '   ', dd[i].slit, '   ', dd[i].slitlen, '   ', dd[i].PA, '   ', dd[i].RA, '   ', dd[i].Dec, '   ', dd[i].imag	
	endfor

	for i=0, n_elements(ee)-1 do begin
        printf, lun, 'CE1       ', ee[i].objname, '   ', ee[i].slit, '   ', ee[i].slitlen, '   ', ee[i].PA, '   ', ee[i].RA, '   ', ee[i].Dec, '   ', ee[i].imag
	endfor

	for i=0, n_elements(ff)-1 do begin
        printf, lun, 'FG1       ', ff[i].objname, '   ', ff[i].slit, '   ', ff[i].slitlen, '   ', ff[i].PA, '   ', ff[i].RA, '   ', ff[i].Dec, '   ', ff[i].imag
	endfor

	for i=0, n_elements(gg)-1 do begin
        printf, lun, 'FG2       ', gg[i].objname, '   ', gg[i].slit, '   ', gg[i].slitlen, '   ', gg[i].PA, '   ', gg[i].RA, '   ', gg[i].Dec, '   ', gg[i].imag
	endfor

	for i=0, n_elements(hh)-1 do begin
        printf, lun, 'GHF1      ', hh[i].objname, '   ', hh[i].slit, '   ', hh[i].slitlen, '   ', hh[i].PA, '   ', hh[i].RA, '   ', hh[i].Dec, '   ', hh[i].imag
	endfor

	for i=0, n_elements(ii)-1 do begin
        printf, lun, 'GHF2      ', ii[i].objname, '   ', ii[i].slit, '   ', ii[i].slitlen, '   ', ii[i].PA, '   ', ii[i].RA, '   ', ii[i].Dec, '   ', ii[i].imag
	endfor

	for i=0, n_elements(jj)-1 do begin
	printf, lun, '16XR1     ', jj[i].objname, '   ', jj[i].slit, '   ', jj[i].slitlen, '   ', jj[i].PA, '   ', jj[i].RA, '   ', jj[i].Dec, '   ', jj[i].imag        
	endfor

	for i=0, n_elements(kk)-1 do begin
        printf, lun, '16XR2     ', kk[i].objname, '   ', kk[i].slit, '   ', kk[i].slitlen, '   ', kk[i].PA, '   ', kk[i].RA, '   ', kk[i].Dec, '   ', kk[i].imag
	endfor

        for i=0, n_elements(ll)-1 do begin
        printf, lun, '16XR3     ', ll[i].objname, '   ', ll[i].slit, '   ', ll[i].slitlen, '   ', ll[i].PA, '   ', ll[i].RA, '   ', ll[i].Dec, '   ', ll[i].imag	
	endfor

free_lun, lun

end
