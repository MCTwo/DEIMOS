pro deepthymee, cplot=cplot, makedat=makedat

	readcol, 'Nzlog4.25_0.1_E.dat', format='D,D,D', deeplstar0, deepphistar0, deeparcnum0, /silent 
	readcol, 'Nzlog4.35_0.1_E.dat', format='D,D,D', deeplstar1, deepphistar1, deeparcnum1, /silent
	readcol, 'Nzlog4.45_0.1_E.dat', format='D,D,D', deeplstar2, deepphistar2, deeparcnum2, /silent
	readcol, 'Nzlog4.55_0.1_E.dat', format='D,D,D', deeplstar3, deepphistar3, deeparcnum3, /silent
	readcol, 'Nzlog4.65_0.1_E.dat', format='D,D,D', deeplstar4, deepphistar4, deeparcnum4, /silent
	readcol, 'Nzlog4.75_0.1_E.dat', format='D,D,D', deeplstar5, deepphistar5, deeparcnum5, /silent
	readcol, 'Nzlog4.85_0.1_E.dat', format='D,D,D', deeplstar6, deepphistar6, deeparcnum6, /silent

        readcol, '../BL12a4.25_0.1.dat', format='D,D,D', lstar0, phistar0, arcnum0, /silent
        readcol, '../BL12a4.35_0.1.dat', format='D,D,D', lstar1, phistar1, arcnum1, /silent
        readcol, '../BL12a4.45_0.1.dat', format='D,D,D', lstar2, phistar2, arcnum2, /silent
        readcol, '../BL12a4.55_0.1.dat', format='D,D,D', lstar3, phistar3, arcnum3, /silent
	readcol, '../BL12a4.65_0.1.dat', format='D,D,D', lstar4, phistar4, arcnum4, /silent
        readcol, '../BL12a4.75_0.1.dat', format='D,D,D', lstar5, phistar5, arcnum5, /silent
        readcol, '../BL12a4.85_0.1.dat', format='D,D,D', lstar6, phistar6, arcnum6, /silent
        readcol, '../BL12a4.95_0.1.dat', format='D,D,D', lstar7, phistar7, arcnum7, /silent
        readcol, '../BL12a5.05_0.1.dat', format='D,D,D', lstar8, phistar8, arcnum8, /silent
        readcol, '../BL12a5.15_0.1.dat', format='D,D,D', lstar9, phistar9, arcnum9, /silent
        readcol, '../BL12a5.25_0.1.dat', format='D,D,D', lstar10, phistar10, arcnum10, /silent
        readcol, '../BL12a5.35_0.1.dat', format='D,D,D', lstar11, phistar11, arcnum11, /silent
        readcol, '../BL12a5.45_0.1.dat', format='D,D,D', lstar12, phistar12, arcnum12, /silent
        readcol, '../BL12a5.55_0.1.dat', format='D,D,D', lstar13, phistar13, arcnum13, /silent
        readcol, '../BL12a5.65_0.1.dat', format='D,D,D', lstar14, phistar14, arcnum14, /silent
        readcol, '../BL12a5.75_0.1.dat', format='D,D,D', lstar15, phistar15, arcnum15, /silent
        readcol, '../BL12a5.85_0.1.dat', format='D,D,D', lstar16, phistar16, arcnum16, /silent
        readcol, '../BL12a5.95_0.1.dat', format='D,D,D', lstar17, phistar17, arcnum17, /silent
        readcol, '../BL12a6.05_0.1.dat', format='D,D,D', lstar18, phistar18, arcnum18, /silent
        readcol, '../BL12a6.15_0.1.dat', format='D,D,D', lstar19, phistar19, arcnum19, /silent
        readcol, '../BL12a6.25_0.1.dat', format='D,D,D', lstar20, phistar20, arcnum20, /silent
        readcol, '../BL12a6.35_0.1.dat', format='D,D,D', lstar21, phistar21, arcnum21, /silent
        readcol, '../BL12a6.45_0.1.dat', format='D,D,D', lstar22, phistar22, arcnum22, /silent

        lstar = lstar0*3.26d*10.d^(42)
        phistar = phistar0*5.5d*10.d^(-3)
        num = fltarr(n_elements(arcnum0))

        for j=0, n_elements(arcnum0)-1 do begin
                num[j] = 12*15.47*(arcnum0[j] + arcnum1[j] +  arcnum2[j] + arcnum3[j] + arcnum4[j] + arcnum5[j] + arcnum6[j] + arcnum7[j] + arcnum8[j] + arcnum9[j] + arcnum10[j] + arcnum11[j] + arcnum12[j] + arcnum13[j] + arcnum14[j] + arcnum15[j] + arcnum16[j] + arcnum17[j] + arcnum18[j] + arcnum19[j] + arcnum20[j] + arcnum21[j] + arcnum22[j])       ;numbers are per arcmin of slit, 15.4 arcmins of slit per mask
	endfor

	deepnum = fltarr(n_elements(deeparcnum2))


	for j=0, n_elements(deeparcnum2)-1 do begin
		deepnum[j] = 83*15.47*(deeparcnum2[j] + deeparcnum3[j] + deeparcnum4[j] + deeparcnum5[j] + deeparcnum6[j]) ;+ deeparcnum0[j]+deeparcnum1[j] )	;numbers are per arcmin of slit, 15.47 arcmins of slit per mask
	endfor

	if n_elements(makedat) then begin

	openw, lun, 'totalnumbers_deep2_a12_4.4z4.9.dat', /get_lun
	
	for i=0, n_elements(arcnum0)-1 do begin
		printf, lun, deeplstar0[i], '   ', deepphistar0[i], '   ', deepnum[i]
	endfor
	free_lun, lun	

	endif

	
	if n_elements(cplot) then begin
	entry_device = !d.name
        loadcolors

	device, decomposed=0
        set_plot, 'PS'
        device, /color, xoffset=0.2, yoffset=0.2, xsize=8, ysize=8, /inches, filename='narf.ps';DEEP2LAE_4.4z4.9_a12_MS_test.ps'

        logphistar = alog10(phistar)
        loglstar = alog10(lstar)

        levels = [1,5,10,20,50,100,200,500]
        c_labels = [1,1,1,1,1,1,1,1]
        c_colors= [7,13,13,13,7,7,7,7]
	xtitle = textoidl('log(L_{*}) erg s^{-1}')
        ytitle = textoidl('log(\Phi_{*}) Mpc^{-3}')
        deepc_colors= [7,4,4,4,7,7,7,7]
	contour, deepnum, loglstar, logphistar, /irregular, /cell_fill, levels=levels, c_colors=c_colors, xrange =[41.2, 43.2], charsize=2, xthick=6, ythick=6, charthick=5, thick=4, xstyle=1, yrange =[-4.45, -2.05], ystyle=1, xtitle=xtitle, ytitle=ytitle
        contour, deepnum, loglstar, logphistar, /irregular, levels=levels, c_charsize=1.5, c_thick=5, c_charthick=4, c_labels=c_labels, /overplot
	contour, num, loglstar, logphistar, /irregular, levels=levels, c_charsize=1.5, c_thick=5, c_charthick=4, color=1, c_labels=c_labels, /overplot
	xyouts, 41.3, -4.32, textoidl('\alpha = -1.2'), charsize=2.2, charthick=4, color=0
        x = [41.3, 41.5]
	y = [-4.15, -4.15]
	oplot, x, y, thick=7, color=1
	xyouts, 41.52, -4.18, 'Cl1604', charsize=1.5, charthick=4, color=1
	oplot, x, y+0.16, thick=7, color=0
	xyouts, 41.52, -3.96, 'DEEP2', charsize=1.5, charthick=4, color=0	
	xyouts, 41.52, -4.05, '4.4<z<4.9', charsize=1.5, charthick=4, color=0

	device, /close_file
        set_plot, entry_device
        endif

	
end






