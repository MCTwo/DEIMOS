pro crystalthymee, cplot=cplot, makedat=makedat

	readcol, 'BL12a4.15_0.1.dat', format='D,D,D', lstar0, phistar0, arcnum0, /silent
	readcol, 'BL12a4.25_0.1.dat', format='D,D,D', lstar1, phistar0, arcnum0, /silent
	readcol, 'BL12a4.35_0.1.dat', format='D,D,D', lstar2, phistar0, arcnum0, /silent
	readcol, 'BL12a4.45_0.1.dat', format='D,D,D', lstar3, phistar0, arcnum0, /silent
	readcol, 'BL12a4.55_0.1.dat', format='D,D,D', lstar4, phistar0, arcnum0, /silent
	readcol, 'BL12a4.65_0.1.dat', format='D,D,D', lstar5, phistar0, arcnum0, /silent 
	readcol, 'BL12a4.75_0.1.dat', format='D,D,D', lstar6, phistar1, arcnum1, /silent
	readcol, 'BL12a4.85_0.1.dat', format='D,D,D', lstar7, phistar2, arcnum2, /silent
	readcol, 'BL12a4.95_0.1.dat', format='D,D,D', lstar8, phistar3, arcnum3, /silent
	readcol, 'BL12a5.05_0.1.dat', format='D,D,D', lstar9, phistar4, arcnum4, /silent
	readcol, 'BL12a5.15_0.1.dat', format='D,D,D', lstar10, phistar5, arcnum5, /silent
	readcol, 'BL12a5.25_0.1.dat', format='D,D,D', lstar11, phistar6, arcnum6, /silent
	readcol, 'BL12a5.35_0.1.dat', format='D,D,D', lstar12, phistar7, arcnum7, /silent
	readcol, 'BL12a5.45_0.1.dat', format='D,D,D', lstar13, phistar8, arcnum8, /silent
	readcol, 'BL12a5.55_0.1.dat', format='D,D,D', lstar14, phistar9, arcnum9, /silent
	readcol, 'BL12a5.65_0.1.dat', format='D,D,D', lstar15, phistar10, arcnum10, /silent
	readcol, 'BL12a5.75_0.1.dat', format='D,D,D', lstar16, phistar11, arcnum11, /silent
	readcol, 'BL12a5.85_0.1.dat', format='D,D,D', lstar17, phistar12, arcnum12, /silent
	readcol, 'BL12a5.95_0.1.dat', format='D,D,D', lstar18, phistar13, arcnum13, /silent
	readcol, 'BL12a6.05_0.1.dat', format='D,D,D', lstar19, phistar14, arcnum14, /silent
	readcol, 'BL12a6.15_0.1.dat', format='D,D,D', lstar20, phistar15, arcnum15, /silent
	readcol, 'BL12a6.25_0.1.dat', format='D,D,D', lstar21, phistar16, arcnum16, /silent
	readcol, 'BL12a6.35_0.1.dat', format='D,D,D', lstar22, phistar17, arcnum17, /silent
	readcol, 'BL12a6.45_0.1.dat', format='D,D,D', lstar23, phistar18, arcnum18, /silent

	lstar = lstar0*3.26d*10.d^(42)
	phistar = phistar0*5.5d*10.d^(-3)
	num = fltarr(n_elements(arcnum0))
	
	for j=0, n_elements(arcnum0)-1 do begin
		num[j] = 12*15.4*(arcnum0[j] + arcnum1[j] +  arcnum2[j] + arcnum3[j] + arcnum4[j] + arcnum5[j] + arcnum6[j] + arcnum7[j] + arcnum8[j] + arcnum9[j] + arcnum10[j] + arcnum11[j] + arcnum12[j] + arcnum13[j] + arcnum14[j] + arcnum15[j] + arcnum16[j] + arcnum17[j] + arcnum18[j]+arcnum14[j] + arcnum15[j] + arcnum16[j] + arcnum17[j] + arcnum18[j])	;numbers are per arcmin of slit, 15.4 arcmins of slit per mask
	endfor

	if n_elements(makedat) then begin

	openw, lun, 'totalnumbers_a16_4.65z6.45.dat', /get_lun
	
	for i=0, n_elements(arcnum0)-1 do begin
		printf, lun, lstar[i], '   ', phistar[i], '   ', num[i]
	endfor
	free_lun, lun	

	endif

	
	if n_elements(cplot) then begin
	entry_device = !d.name
        loadcolors

	device, decomposed=0
        set_plot, 'PS'
        device, /color, xoffset=0.2, yoffset=0.2, xsize=8, ysize=8, /inches, filename='cl1604LAE_a12_MS.ps'

        logphistar = alog10(phistar)
        loglstar = alog10(lstar)

        levels = [1,4,5,10,20,30,36,40,50,100,200,500]
        c_labels = [1,0,1,1,1,1,0,1,1,1,1]
        c_colors= [7,13,13,13,13,13,7,7,7,7,7]
	xtitle = textoidl('log(L_{*}) erg s^{-1}')
        ytitle = textoidl('log(\Phi_{*}) Mpc^{-3}')
        contour, num, loglstar, logphistar, /irregular, /cell_fill, levels=levels, c_colors=c_colors, xrange =[41.2, 43.2], charsize=2, xthick=6, ythick=6, charthick=5, thick=4, xstyle=1, yrange =[-4.45, -2.05], ystyle=1, xtitle=xtitle, ytitle=ytitle
        contour, num, loglstar, logphistar, /irregular, levels=levels, c_charsize=1.5, c_thick=5, c_charthick=4, c_labels=c_labels, /overplot
        xyouts, 41.4, -4.2, textoidl('\alpha = -1.2'), charsize=2.2, charthick=4, color=0
        device, /close_file
        set_plot, entry_device
        endif

	
end






