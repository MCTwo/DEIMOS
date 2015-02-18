; define list structure obj_info files  

pro objlist__define

    tmp = {objlist, $
	   objno:'', $
	   slitno:0L, $
	   slitfile: '', $
	   objtype: '', $
	   color: '', $
           CAT_OBJPOS: 0.0, $ 
	   CAT_FWHM: 0.0, $
	   OBJPOS:0.0, $
	   FWHM:0.0, $ 
	   CORR_FWHM:0.0, $
	   NROWS:0.0, $
	   S2N_FWHM:0.0, $
	   S2N_WINDOW:0.0,$ 
	   RA:'',$ 
	   DEC:'',$ 
	   XMM:0.0,$ 
	   YMM:0.0,$ 
	   MAGB:0.0,$ 
	   MAGR:0.0,$ 
	   MAGI:0.0, $
	   OBJPA:0.0}
end     
