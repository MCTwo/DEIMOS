pro smooth_ivar, arrayin, ivarin, sigma, arrayout, ivarout
     array = arrayin
     ivar = ivarin
     if n_elements(array) ne n_elements(ivar) then begin
         print, 'Wrong array dimensions in SMOOTH_IVAR.'
     endif
     asize = n_elements(array)
     n = floor(8. * sigma) + 1
     if n gt asize then begin
         message, 'Array too small to smooth.'
     endif
     half_n = (n-1)/2
     arraysm = dblarr(asize)
     ivarsm = dblarr(asize)
     ;window = dblarr(n) + 1d                                                                                   ;boxcar
     window = (1. / (double(sigma)) * sqrt(2.*!PI)) * exp(-(indgen(n)-double(half_n))^2./(2.*double(sigma)^2.))   ;Gaussian
     for i=0,half_n-1 do begin
         if total(ivar[0:i+half_n]*window[n-1-i-half_n:n-1]) ne 0.0 then begin
             arraysm[i] = total(array[0:i+half_n]*ivar[0:i+half_n]*window[n-1-i-half_n:n-1]) / total(ivar[0:i+half_n]*window[n-1-i-half_n:n-1])
             ivarsm[i] = (total(ivar[0:i+half_n]*window[n-1-i-half_n:n-1]))^2. / total(ivar[0:i+half_n] * (window[n-1-i-half_n:n-1])^2.)
         endif
     endfor
     for i=half_n,asize-half_n-1 do begin
         if total(ivar[i-half_n:i+half_n]*window) ne 0.0 then begin
             arraysm[i] = total(array[i-half_n:i+half_n]*ivar[i-half_n:i+half_n]*window) / total(ivar[i-half_n:i+half_n]*window)
             ivarsm[i] = (total(ivar[i-half_n:i+half_n]*window))^2. / total(ivar[i-half_n:i+half_n] * window^2.)
         endif
     endfor
     for i=asize-half_n,asize-1 do begin
         if total(ivar[i-half_n:asize-1]*window[0:asize-i-1+half_n]) ne 0.0 then begin
             arraysm[i] = total(array[i-half_n:asize-1]*ivar[i-half_n:asize-1]*window[0:asize-i-1+half_n]) / total(ivar[i-half_n:asize-1]*window[0:asize-i-1+half_n])
             ivarsm[i] = (total(ivar[i-half_n:asize-1]*window[0:asize-i-1+half_n]))^2. / total(ivar[i-half_n:asize-1] * (window[0:asize-i-1+half_n])^2.)
         endif
     endfor
     arrayout = arraysm
     ivarout = ivarsm
end

