; define zresult structure for spec1d pipeline.
pro zresult__define

  tmp = {zresult, $
         class:'', $
         subclass:'', $
         objname:'', $
         slitname:'', $
         maskname:'', $    
         date:'', $
         mjd:0.0D0, $
         z:0.0, $
         z_err:0.0, $
         rchi2:0.0, $
         dof:0L, $
         rchi2diff:0.0, $
         tfile:'', $
         tcolumn:lonarr(10) - 1L, $
         npoly:0L, $
         theta:fltarr(10), $
         vdisp:0.0, $
         vdisp_err:0.0, $
         zquality:0, $
         comment:''}

end








