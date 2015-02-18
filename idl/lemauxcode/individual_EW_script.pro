pro individual_ew_script, cat

   readcol, cat,  format='A,A,D,D,A,A', mask, slit, Iband, z, q, file, skipline=0
   numread = n_elements(mask)
   numgoo = numread+1
   EWlist = replicate({EWlist}, numgoo)
   EWOII = dblarr(numgoo)
   EWHd = dblarr(numgoo)
   result1 = dblarr(numgoo)
   result2 = dblarr(numgoo)

   for i=0, numread-1 do begin
   EWlist[i] = {EWlist, maskname:mask[i], slitname:slit[i], Iband:Iband[i], z:z[i], EWOII:EWOII[i], EWHd:EWHd[i], zquality:q[i], file:file[i]}
   endfor
   
   for i=0, numread-1 do begin
   spec = mrdfits(EWlist[i].file,1)
   result1[i] = measure_simple_ew(spec,$
                [[3696.3,3716.3],[3738.3,3758.3]],$
                        [3716.3,3738.3])
   result2[i] = measure_simple_ew(spec,$
		[[4017.,4057.],[4153.,4193.]], $
			[4083.5,4122.250])
   EWlist[i].EWOII = result1[i]
   EWlist[i].EWHd = result2[i]
   endfor		
   
   openw, lun, 'output.dat', /get_lun, width = 400
   for i=0, numread-1 do begin
   printf, lun, EWlist[i].maskname, EWlist[i].slitname, EWlist[i].Iband, EWlist[i].z, EWlist[i].EWOII, EWlist[i].EWHd, EWlist[i].zquality, EWlist[i].file, format='(a6, 2x, i3, 2x, d10.6, 2x, d10.6, 2x, d12.6, 2x, d10.6, 2x, a1, 2x, a50)'
   endfor
   free_lun, lun

end  
