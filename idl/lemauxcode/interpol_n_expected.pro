pro interpol_n_expected, alpha

	if alpha eq 1.2 or alpha eq '1.2' or alpha eq -1.2 or alpha eq '-1.2' then begin

	readcol, '../BL12a4.65_0.1.dat', format='D,D,D', lstar0, phistar0, arcnum0, /silent
        readcol, '../BL12a4.75_0.1.dat', format='D,D,D', lstar1, phistar1, arcnum1, /silent
        readcol, '../BL12a4.85_0.1.dat', format='D,D,D', lstar2, phistar2, arcnum2, /silent
        readcol, '../BL12a4.95_0.1.dat', format='D,D,D', lstar3, phistar3, arcnum3, /silent
        readcol, '../BL12a5.05_0.1.dat', format='D,D,D', lstar4, phistar4, arcnum4, /silent
        readcol, '../BL12a5.15_0.1.dat', format='D,D,D', lstar5, phistar5, arcnum5, /silent
        readcol, '../BL12a5.25_0.1.dat', format='D,D,D', lstar6, phistar6, arcnum6, /silent
        readcol, '../BL12a5.35_0.1.dat', format='D,D,D', lstar7, phistar7, arcnum7, /silent
        readcol, '../BL12a5.45_0.1.dat', format='D,D,D', lstar8, phistar8, arcnum8, /silent
        readcol, '../BL12a5.55_0.1.dat', format='D,D,D', lstar9, phistar9, arcnum9, /silent
        readcol, '../BL12a5.65_0.1.dat', format='D,D,D', lstar10, phistar10, arcnum10, /silent
        readcol, '../BL12a5.75_0.1.dat', format='D,D,D', lstar11, phistar11, arcnum11, /silent
        readcol, '../BL12a5.85_0.1.dat', format='D,D,D', lstar12, phistar12, arcnum12, /silent
        readcol, '../BL12a5.95_0.1.dat', format='D,D,D', lstar13, phistar13, arcnum13, /silent
        readcol, '../BL12a6.05_0.1.dat', format='D,D,D', lstar14, phistar14, arcnum14, /silent
        readcol, '../BL12a6.15_0.1.dat', format='D,D,D', lstar15, phistar15, arcnum15, /silent
        readcol, '../BL12a6.25_0.1.dat', format='D,D,D', lstar16, phistar16, arcnum16, /silent
        readcol, '../BL12a6.35_0.1.dat', format='D,D,D', lstar17, phistar17, arcnum17, /silent
        readcol, '../BL12a6.45_0.1.dat', format='D,D,D', lstar18, phistar18, arcnum18, /silent
	endif 

	if alpha eq 1.6 or alpha eq '1.6' or alpha eq -1.6 or alpha eq '-1.6' then begin

	readcol, '../BL16a4.65_0.1.dat', format='D,D,D', lstar0, phistar0, arcnum0, /silent
        readcol, '../BL16a4.75_0.1.dat', format='D,D,D', lstar1, phistar1, arcnum1, /silent
        readcol, '../BL16a4.85_0.1.dat', format='D,D,D', lstar2, phistar2, arcnum2, /silent
        readcol, '../BL16a4.95_0.1.dat', format='D,D,D', lstar3, phistar3, arcnum3, /silent
        readcol, '../BL16a5.05_0.1.dat', format='D,D,D', lstar4, phistar4, arcnum4, /silent
        readcol, '../BL16a5.15_0.1.dat', format='D,D,D', lstar5, phistar5, arcnum5, /silent
        readcol, '../BL16a5.25_0.1.dat', format='D,D,D', lstar6, phistar6, arcnum6, /silent
        readcol, '../BL16a5.35_0.1.dat', format='D,D,D', lstar7, phistar7, arcnum7, /silent
        readcol, '../BL16a5.45_0.1.dat', format='D,D,D', lstar8, phistar8, arcnum8, /silent
        readcol, '../BL16a5.55_0.1.dat', format='D,D,D', lstar9, phistar9, arcnum9, /silent
        readcol, '../BL16a5.65_0.1.dat', format='D,D,D', lstar10, phistar10, arcnum10, /silent
        readcol, '../BL16a5.75_0.1.dat', format='D,D,D', lstar11, phistar11, arcnum11, /silent
        readcol, '../BL16a5.85_0.1.dat', format='D,D,D', lstar12, phistar12, arcnum12, /silent
        readcol, '../BL16a5.95_0.1.dat', format='D,D,D', lstar13, phistar13, arcnum13, /silent
        readcol, '../BL16a6.05_0.1.dat', format='D,D,D', lstar14, phistar14, arcnum14, /silent
        readcol, '../BL16a6.15_0.1.dat', format='D,D,D', lstar15, phistar15, arcnum15, /silent
        readcol, '../BL16a6.25_0.1.dat', format='D,D,D', lstar16, phistar16, arcnum16, /silent
        readcol, '../BL16a6.35_0.1.dat', format='D,D,D', lstar17, phistar17, arcnum17, /silent
        readcol, '../BL16a6.45_0.1.dat', format='D,D,D', lstar18, phistar18, arcnum18, /silent

	endif
	
	if n_elements(lstar0) le 0 then message, "Faint end slope must be set to either 1.2 or 1.6"

	z = 4.65 + 0.1*findgen(19)
	zplot = 4.05 + 0.1*findgen(30)
	N = indgen(19)
	n_exp_suppository = fltarr(19,n_elements(lstar0))	;set up array that will contain N_exepected for each value of Lstar and phistar in each row for each redshift bin
	LFfit_suppository = fltarr(6,n_elements(lstar0))	;set up array that will contain 3rd order (5th was a better fit but underestimates a lot at lower z, 3rd order is pretty good) polynomial fit parameters for each value of phistar and lstar for each redshift bin
	z415 = dblarr(n_elements(lstar0))
	z425 = dblarr(n_elements(lstar0))
	z435 = dblarr(n_elements(lstar0))
	z445 = dblarr(n_elements(lstar0))
	z455 = dblarr(n_elements(lstar0))

	n_exp_suppository[0,0:n_elements(lstar0)-1] = arcnum0	
	n_exp_suppository[1,0:n_elements(lstar0)-1] = arcnum1
	n_exp_suppository[2,0:n_elements(lstar0)-1] = arcnum2
	n_exp_suppository[3,0:n_elements(lstar0)-1] = arcnum3
	n_exp_suppository[4,0:n_elements(lstar0)-1] = arcnum4			
	n_exp_suppository[5,0:n_elements(lstar0)-1] = arcnum5
	n_exp_suppository[6,0:n_elements(lstar0)-1] = arcnum6
	n_exp_suppository[7,0:n_elements(lstar0)-1] = arcnum7
	n_exp_suppository[8,0:n_elements(lstar0)-1] = arcnum8
	n_exp_suppository[9,0:n_elements(lstar0)-1] = arcnum9
	n_exp_suppository[10,0:n_elements(lstar0)-1] = arcnum10
	n_exp_suppository[11,0:n_elements(lstar0)-1] = arcnum11
	n_exp_suppository[12,0:n_elements(lstar0)-1] = arcnum12
	n_exp_suppository[13,0:n_elements(lstar0)-1] = arcnum13
	n_exp_suppository[14,0:n_elements(lstar0)-1] = arcnum14
	n_exp_suppository[15,0:n_elements(lstar0)-1] = arcnum15
	n_exp_suppository[16,0:n_elements(lstar0)-1] = arcnum16
	n_exp_suppository[17,0:n_elements(lstar0)-1] = arcnum17
	n_exp_suppository[18,0:n_elements(lstar0)-1] = arcnum18

	openw, lun, 'n_expected_realdata_z4.65toz6.45' + string(strcompress(alpha,/REMOVE_ALL)) + '.dat', /get_lun, width=500
	
		for i=0, n_elements(lstar0)-1 do begin
			printf, lun, n_exp_suppository[*,i]
		endfor
	
	free_lun, lun

	entry_device = !d.name
        set_plot, 'PS'
	device, /color, xoffset=0.2, yoffset=0.2, xsize=8, ysize=8, /inches, filename='LFinterpolate_to4.05' + string(strcompress(alpha,/REMOVE_ALL)) + '.ps'	
	loadcolors

	for i=0, n_elements(lstar0) - 1 do begin
	
		n_expected = fltarr(19)
		for j=0, 19-1 do begin
			n_expected[j] = n_exp_suppository[j,i] 
		endfor		

		LFfit = poly_fit(z,n_expected,3)
		LFfit_suppository[0:3,i] = LFfit

		if i eq 0 then begin
			plot, zplot, LFfit[0] + LFfit[1]*zplot + LFfit[2]*zplot^2 + LFfit[3]*zplot^3, position = [0.17, 0.1, 0.93, 0.93], /normal, yrange=[1.5e-07,3e-05], ystyle=1, xrange=[4,6.5], xstyle=1, xtitle='z', ytitle='N_expected', charthick=5, thick=4, charsize=1.5, /nodata ;+ LFfit[4]*zplot^4 +  LFfit[5]*zplot^5 , yrange=[1.5e-07,5e-06], ystyle=1, xrange=[4,6.5], xstyle=1, /nodata
			oplot, zplot,  LFfit[0] + LFfit[1]*zplot + LFfit[2]*zplot^2 + LFfit[3]*zplot^3, color=1, thick=4;+ LFfit[4]*zplot^4 +  LFfit[5]*zplot^5, color=1
			oplot, z, n_expected, psym=4, color=1, symsize=1.7
		endif			

		if i ne 0 then begin
			oplot, zplot, LFfit[0] + LFfit[1]*zplot + LFfit[2]*zplot^2 + LFfit[3]*zplot^3, color=i+1, thick=3;+  LFfit[4]*zplot^4 +  LFfit[5]*zplot^5, color=i+1
	 		oplot, z, n_expected, psym=4, color=i+1, symsize=1.7
		endif	
			
		z415[i] = LFfit[0] + LFfit[1]*4.15 + LFfit[2]*4.15^2 + LFfit[3]*4.15^3 ;+ LFfit[4]*4.15^4 +  LFfit[5]*4.15^5
		z425[i] = LFfit[0] + LFfit[1]*4.25 + LFfit[2]*4.25^2 + LFfit[3]*4.25^3 ;+ LFfit[4]*4.25^4 +  LFfit[5]*4.25^5
		z435[i] = LFfit[0] + LFfit[1]*4.35 + LFfit[2]*4.35^2 + LFfit[3]*4.35^3 ;+ LFfit[4]*4.35^4 +  LFfit[5]*4.35^5
		z445[i] = LFfit[0] + LFfit[1]*4.45 + LFfit[2]*4.45^2 + LFfit[3]*4.45^3 ;+ LFfit[4]*4.45^4 +  LFfit[5]*4.45^5
	 	z455[i] = LFfit[0] + LFfit[1]*4.55 + LFfit[2]*4.55^2 + LFfit[3]*4.55^3 ;+ LFfit[4]*4.55^4 +  LFfit[5]*4.55^5
	
	endfor	

	xyouts, 4.1, 9e-07, textoidl('\alpha=') + string(strcompress(alpha,/REMOVE_ALL)), charsize=1.6, charthick=6
	device, /close_file
	set_plot, entry_device
	alpha = string(strcompress(alpha,/REMOVE_ALL))
	
	openw, lun, 'LFfits.fromz4.65toz6.45.' + alpha + '.dat'
	
	for i=0, n_elements(lstar0)-1 do begin
		printf, lun, LFfit_suppository[*,i]
	endfor

	free_lun, lun

	openw, 1, 'BL' + alpha + 'a4.15_0.1.dat'
	openw, 2, 'BL' + alpha + 'a4.25_0.1.dat'
        openw, 3, 'BL' + alpha + 'a4.35_0.1.dat'    
	openw, 4, 'BL' + alpha + 'a4.45_0.1.dat'
	openw, 5, 'BL' + alpha + 'a4.55_0.1.dat' 

	for i=0, n_elements(lstar0)-1 do begin
		printf, 1, lstar0[i], '    ', phistar0[i], '    ', z415[i], '    ', z415[i], format = '(E11.5, a4, E11.5, a4, E12.5, a4, E12.5)'
		printf, 2, lstar0[i], '    ', phistar0[i], '    ', z425[i], '    ', z425[i], format = '(E11.5, a4, E11.5, a4, E12.5, a4, E12.5)'
		printf, 3, lstar0[i], '    ', phistar0[i], '    ', z435[i], '    ', z435[i], format = '(E11.5, a4, E11.5, a4, E12.5, a4, E12.5)'
             	printf, 4, lstar0[i], '    ', phistar0[i], '    ', z445[i], '    ', z445[i], format = '(E11.5, a4, E11.5, a4, E12.5, a4, E12.5)'
		printf, 5, lstar0[i], '    ', phistar0[i], '    ', z455[i], '    ', z455[i], format = '(E11.5, a4, E11.5, a4, E12.5, a4, E12.5)' 
	endfor

	free_lun, 1,2,3,4,5
end

		
	
		
