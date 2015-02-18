

; Definitions of Bandpass/sideband from Fisher et al
FUNCTION measure_all_ews

	zz = read_kzcat()
	roi = where(zz.zquality ge 3)
	nz = n_elements(roi)

	v = [0.D, 0.D]
	ews = {objno: 0L, oii: v, halpha: v, hbeta: v, hgamma: v, hdelta: v, Dn4000: v}
	ews = replicate(ews, nz)

	tick
	for i = 0, nz - 1 do begin
		zt = zz[roi[i]]
		ews[i].objno = long(zt.objname)
		npk_readspec, strcompress(zt.objname,/rem), strcompress(zt.maskname,/rem),$
			strcompress(zt.slitname,/rem), /silent,$
			lambda = lambda, spec = ss, ivar = ivar
		spec = {lambda: lambda/(1.+zt.z), spec: ss}

		ews[i].oii = measure_simple_ew(spec,$
			[[3696.3,3716.3],[3738.3,3758.3]],$
			[3716.3,3738.3])
		ews[i].hdelta = measure_simple_ew(spec, $
			[[4017.,4057.],[4153.,4193.]], $
			[4083.5,4122.250])
		ews[i].hgamma = measure_simple_ew(spec, $
			[[4242.,4282.],[4404.,4444.]], $
			[4319.75,4363.5])
		ews[i].hbeta = measure_simple_ew(spec, $
			[[4799.,4839.],[4886.,4926.]], $
			[4847.875,4876.625])

		; Defined in Kauffmann astro-ph/0204055
		ews[i].dn4000 = measure_simple_ratio(spec, $
			[3850,3950],[4000,4100])

		if i mod 100 eq 0 then begin
			print, 'At i ', i, ' out of ', nz
			tock
			tick
		endif
	endfor
	return,ews
end

