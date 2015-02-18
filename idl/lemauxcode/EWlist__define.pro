; define list structure for measuring EWs of Hd and OII

pro EWlist__define

    tmp = {EWlist, $
	   maskname:'', $
	   slitname:0, $
	   Iband: 0.0, $
	   z: 0.0, $
           EWOII: 0.0, $
	   EWHd: 0.0, $
	   zquality: '', $
	   file:''}
end     
