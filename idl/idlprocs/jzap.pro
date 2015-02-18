;+
; NAME:
;      JZAP
;
; PURPOSE:
;      Kills cosmic rays dead. Designed for echelle orders, but works
;      for images as well.
;
; CALLING SEQUENCE:
;      jzap, oldorder [, neworder , nsig=, mask=, skymask=, 
;                        smoothim=, /noreplace ]
;
; INPUTS:
;      OLDRORDER - Original image with cosmic rays.
;
; KEYWORD PARAMETERS:
;      NSIG      - Number of standard deviations above it's neighbors
;                  a pixel must stand in order to be considered a 
;                  cosmic ray candidate.
;
;      MASK      - Cosmic ray mask returned through this variable. Has
;                  same dimensions as OLDORDER.
;
;      SKYMASK   - Vector containing sky line mask. Length is equal to
;                  the number of columns in OLDORDER.
;
;      SMOOTHIM  - The median-smoothed image of the order is returned
;                  through this keyword.
;
;      NOREPLACE - By default, cosmic rays are replaced by the
;                  corresponding pixel in a median smoothed image of
;                  the order. /NOREPLACE does not perform this
;                  step. Use NOREPLACE when JZAP is used in
;                  conjunction with an Optimal Extraction algorthim.
;
; OPTIONAL OUTPUTS:
;      NEWORDER  - New image sans cosmic rays.
;
; SIDE EFFECTS:
;      Replaces cosmic rays with the median of the neighborhood unless
;      the NOREPLACE keyword is set.
;
; RESTRICTIONS:
;      Uses a while loop to iterate. May be kinda slow. But I haven't 
;      found a better way to get rid of cosmics in a rectified
;      exchelle order.
;
;      JZAP is most sensitive in the dispersion direction and works
;      best with low-signal spectra. (Cosmic rays in the MIKE
;      spectrograph run preferentially in the cross-dispersion direction)
;
; PROCEDURE:
;      All peaks in the image are located in the dispersion (row-wise)
;      direction. Cosmic rays are identified as those peaks which
;      stand out NSIG standard deviations above their neighbors on the
;      left *and* right. This prevents line features from being
;      identified as cosmics. Emission lines may stand out from the spetrum,
;      but each pixel in a line will not A) be a peak and B) stand out
;      NSIG from its neighbors on both sides. This process is repeated
;      in the cross-dispersion direction. 
;
; EXAMPLE:
;
; IDL> im = randomn(seed,200,20)
; IDL> im[randomu(seed,8)*2000+2000] = randomu(seed,8)*50+20
; IDL> jzap,im,nim
; IDL> display,im,max=10
; IDL> display,nim,max=10
;
; MODIFICATION HISTORY:
; Written by JohnJohn in April 2003
; 23 May 2003  JohnJohn   Added MASK and SKYMASK keywords.
; 04 June 2003 JohnJohn   Added NOREPLACE, SMOOTHIM and fixed error on
; line 113 where D3 was supposed to be D4.
;-

pro zapper,spec,newim,mask,transpose=transpose,nsig=nsig,noreplace=norep,smoothim=sm
on_error,2
if keyword_set(transpose) then begin
    im = transpose(spec) 
    mask = transpose(mask)
    if n_elements(sm) gt 0 then sm = transpose(sm)
endif else im = spec

cool = 0
c = 0
nfix = 0
prev = 0
nel = n_elements(im)
sz = size(spec)
if keyword_set(norep) then begin
    sm = spec
endif else begin
    if sz[0] eq 1 then begin
        fat = fltarr(32)
        sm = median([fat,spec,fat],31)
        sm = sm[32:n_elements(spec)+32-1]
    endif else if n_elements(sm) eq 0 then sm = filter_image(im,median=5,/all)
endelse

while not cool do begin
    imr = reform(im,nel)
    d1 = imr - shift(imr,1)
    d2 = imr - shift(imr,-1)
    peak = where(d1 gt 0 and d2 gt 0)

    d3 = im[peak] - im[peak-1]
    sig3 = stdev(d3)
    cond3 = d3-mean(d3) gt nsig*sig3
    d4 = im[peak] - im[peak+1]
    sig4 = stdev(d4)
    cond4 = d4-mean(d4) gt nsig*sig4

    bad = where(cond3 and cond4, nbad)
    if nbad gt 0 then begin
        im[peak[bad]] = sm[peak[bad]]
        mask[peak[bad]] = 0.
    endif

    c = c + 1
    if nbad eq 0 or c gt 10 or nbad eq prev then cool = 1
    prev = nbad
endwhile
if keyword_set(transpose) then begin
    newim = transpose(im) 
    mask = transpose(mask)
    sm = transpose(sm)
endif else newim = im
end

;;;;;JZAP is the driver for ZAPPER
pro jzap,spec,newim,nsig=nsig,mask=mask,skymask=skymask,noreplace=norep,smoothim=smoothim
on_error,2
sz = size(spec)
ncol = sz[1]
if n_elements(nsig) eq 0 then nsig = 5
if (size(spec))[0] eq 1 then begin
    mask = fltarr(ncol)+1
    zapper,spec,newim,mask,nsig=nsig,norep=norep
endif else begin
    nrow = sz[2]
    mask = fltarr(ncol,nrow)+1
    if keyword_set(skymask) then begin
        sk = where(skymask,ct) 
        if ct eq 0 then sk = findgen(ncol)
    endif else begin     
        sk = findgen(ncol)
    endelse
    s = spec[sk,*]
    m = mask[sk,*]
    for i = 0,1 do begin
        zapper,s,new,m,nsig=nsig,/transpose,norep=norep,smoothim=smoothim
        zapper,new,new,m,nsig=nsig,norep=norep,smoothim=smoothim
    endfor
    newim = spec
    newim[sk,*] = new
    mask[sk,*] = m
endelse
end
