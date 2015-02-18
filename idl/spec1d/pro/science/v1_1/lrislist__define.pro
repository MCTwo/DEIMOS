; define list structure for DEEP2 LRIS coaddition routine

pro lrislist__define

    tmp = {lrislist, $
	   maskname:'', $
	   slitname:'', $
	   Iband: 0.0, $
	   z: 0.0, $
	   zquality: '', $
	   file:'', $
	   absIflux: 0.0}
end     
