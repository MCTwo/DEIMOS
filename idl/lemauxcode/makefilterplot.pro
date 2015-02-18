pro makefilterplot
	
	readcol, 'V.dat', vlam, v, /silent
	v = v/1.65
	readcol, 'r.dat', rlam, r, /silent
	readcol, 'i.dat', ilam, i, /silent
	readcol, 'z.dat', zlam, z, /silent
	z = z*2
	readcol, 'f606w.dat', f606wlam, f606w, /silent
	readcol, 'f814w.dat', f814wlam, f814w, /silent
	lyalowz = makesyntheticLya(4.4)
	lyamedz = makesyntheticLya(4.8)
	lyahighz = makesyntheticLya(5.6) 

	set_plot, 'PS'
        device, /color, xoffset=0.2, yoffset=0.2, xsize=8, ysize=11, /inches, filename='filterplot_wLYA.ps'
	loadcolors
	lowzrange = [5600, 7900]
	medzrange = [6100,8400]
	highzrange = [7050, 9490]
	;device, decomposed=0	
	
	;low redshift LAE
	
	plot, vlam, v, /nodata, position = [0.10, 0.55, 0.36, 0.95], color=0, xstyle=1, xtickname = replicate(' ', 5), xrange = lowzrange,$
	yrange = [0,0.8], ystyle = 4, xthick=6
	oplot, vlam, v, color=2, linestyle=1, thick=9
        oplot, rlam, r, color=8, linestyle=2, thick=4
        oplot, ilam, i, color=4, linestyle=0, thick=4
        oplot, zlam, z, color=5, linestyle=3, thick=4
	oplot, lyalowz.lambda, lyalowz.spec/4400, color=0, thick=2
	xyouts, 5700, 0.74, 'z=4.4', charthick=5, charsize=2, color=0
	axis, yaxis=0, yrange=[0,0.8], ythick=6, ystyle=1, color=0, ytitle='Total Throughput', ytickname = [' ','0.2','0.4','0.6', '0.8','1'], charsize=1.5, charthick=5
        axis, yaxis=1, yrange=[0,119], ythick=6, ystyle=1, color=0,  ytickname=replicate(' ', 7);, thick=5 ;ytitle=textoidl('Flux Density f_{\lambda} (Arbitrary Units)')
	
	plot, f606wlam, f606w, /nodata, position = [0.10, 0.17, 0.36, 0.55], color=0, xstyle=1, xtitle = textoidl('\lambda_{Obs} (\AA)'), xtickname=['6000', ' ', '7000', ' '], xrange = lowzrange, $
        ystyle = 4, yrange=[0,0.8], charsize = 1.5, xthick=6, charthick=5, /noerase
	oplot, f606wlam, f606w, color=6, thick=4
	oplot, f814wlam, f814w, color=1, linestyle=2, thick=4
	oplot, lyalowz.lambda, lyalowz.spec/4400, color=0, thick=2
	axis, yaxis=0, yrange=[0,0.8], color=0, ythick=6, ystyle=1, ytitle='Total Throughput', ytickname = ['0','0.2','0.4','0.6', '0.8','1'], charsize=1.5, charthick=5
        axis, yaxis=1, yrange=[0,119], ystyle=2, ythick=6, color=0, ytickname=replicate(' ', 7);, thick=5 ;ytitle=textoidl('Flux Density f_{\lambda} (Arbitrary Units)')
	
	;medium redshift LAE
	plot, vlam, v, /nodata, position = [0.36, 0.55, 0.62, 0.95], color=0, xstyle=1, xtickname = replicate(' ', 5), xrange = medzrange,$
        yrange = [0,0.8], ystyle = 4, xthick=6, /noerase
        oplot, vlam, v, color=2, linestyle=1, thick=9
        oplot, rlam, r, color=8, linestyle=2, thick=4
        oplot, ilam, i, color=4, linestyle=0, thick=4
        oplot, zlam, z, color=5, linestyle=3, thick=4
        oplot, lyamedz.lambda, lyamedz.spec/4400, color=0, thick=2
	;legend stuff
	x1 = [6800,7500]
	y1 = [0.71,0.71]
	y2 = [0.66,0.66]
	y3 = [0.61,0.61]
	y4 = [0.56, 0.56]
	oplot, x1, y1, color=2, linestyle=1, thick=9
	xyouts, x1[1]+50, y1[0]-0.01, 'V', color=2, charsize=1.75, charthick=4
	oplot, x1, y2, color=8, linestyle=2, thick=4
	xyouts, x1[1]+50, y2[0]-0.01, textoidl('r^\prime'), color=8, charsize=1.75, charthick=3
	oplot, x1, y3, color=4, linestyle=0, thick=4
        xyouts, x1[1]+50, y3[0]-0.01, textoidl('i^\prime'), color=4, charsize=1.75, charthick=3
	oplot, x1, y4, color=5, linestyle=3, thick=4
        xyouts, x1[1]+50, y4[0]-0.01, textoidl('z^\prime'), color=5, charsize=1.75, charthick=3
       	xyouts, 6200, 0.74, 'z=4.8', charthick=5, charsize=2, color=0
	;finish up axes 
	axis, yaxis=0, yrange=[0,0.8], ythick=6, color=0, ystyle=1, ytickname=replicate(' ', 6);, thick=5
        axis, yaxis=1, yrange=[0,119], ythick=6, ystyle =1, color=0,  ytickname=replicate(' ', 7);, thick=5 ;ytitle=textoidl('Flux Density f_{\lambda} (Arbitrary Units)')
	

        plot, f606wlam, f606w, /nodata, position = [0.36, 0.17, 0.62, 0.55], color=0, xstyle=1, xtitle = textoidl('\lambda_{Obs} (\AA)'), xtickname=[' ', '7000', ' ', '8000'],$
	xrange = medzrange, ystyle = 4, yrange=[0,0.8], charsize=1.5, xthick=6, charthick=5, /noerase
	oplot, f606wlam, f606w, color=6, thick=4
        oplot, f814wlam, f814w, color=1, linestyle=2, thick=4
        oplot, lyamedz.lambda, lyamedz.spec/4400, color=0, thick=2 
       
	;legend stuff 
	oplot, x1, y1, color=6, thick=4
        xyouts, x1[1]+50, y1[0]-0.015, 'F606W', color=6, charsize=1.5, charthick=3
        oplot, x1, y2, color=1, linestyle=2, thick=4
        xyouts, x1[1]+50, y2[0]-0.015, 'F814W', color=1, charsize=1.5, charthick=3

	;finishing up axes
	axis, yaxis=0, yrange=[0,0.8], ythick=6, color=0, ystyle=1,  ytickname=replicate(' ', 6);, thick=5
        axis, yaxis=1, yrange=[0,119], ythick=6, ystyle=1, color=0, ytickname=replicate(' ', 7);, thick=5

	;high redshift LAE
	plot, vlam, v, /nodata, position = [0.62, 0.55, 0.89, 0.95], color=0, xstyle=1, xtickname = replicate(' ', 5), xrange = highzrange,$
        yrange = [0,0.8], ystyle = 4, xthick=6, /noerase
        oplot, vlam, v, color=2, linestyle=1, thick=4
        oplot, rlam, r, color=8, linestyle=2, thick=4
        oplot, ilam, i, color=4, linestyle=0, thick=4
        oplot, zlam, z, color=5, linestyle=3, thick=4
        oplot, lyahighz.lambda, lyahighz.spec/4400, color=0, thick=2
       	xyouts, 7200, 0.74, 'z=5.6', charthick=5, charsize=2, color=0
	axis, yaxis=0, yrange=[0,0.8], ythick=6, color=0, ystyle=1, ytickname=replicate(' ', 6);, thick=5
        axis, yaxis=1, yrange=[0,119], ythick=6, ystyle =1, color=0,  ytitle=textoidl('Flux Density f_{\lambda} (Arbitrary Units)'), charsize=1.5, charthick=5

        plot, f606wlam, f606w, /nodata, position = [0.62, 0.17, 0.89, 0.55], color=0, xstyle=1,xtitle = textoidl('\lambda_{Obs} (\AA)'), xtickname=[' ', '8000', ' ', '9000'],$ 
	xrange = highzrange, ystyle = 4, yrange=[0,0.8], charsize=1.5, xthick=6, charthick=5, /noerase
        oplot, f606wlam, f606w, color=6, thick=4
        oplot, f814wlam, f814w, color=1, thick=4, linestyle=2
        oplot, lyahighz.lambda, lyahighz.spec/4400, color=0, thick=2 
        axis, yaxis=0, yrange=[0,0.8], ythick=6, ystyle=1, color=0, ytickname=replicate(' ', 6);, thick=5
        axis, yaxis=1, yrange=[0,119], ythick=6, ystyle=1, color=0, ytitle=textoidl('Flux Density f_{\lambda} (Arbitrary Units)'), charsize=1.5, charthick=5
	device, /close_file

end	
	
