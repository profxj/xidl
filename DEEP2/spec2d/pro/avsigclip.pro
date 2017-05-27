function avsigclip, array, invvar, threshold, sigmap, $
                    inmask=inmask,niteration=niteration
;+
; NAME:
;   avsigclip
;
; PURPOSE:
;   averages 3rd dimension of a 3d array to reject CR outliers
;      in mean average of 2d array
;
; CALLING SEQUENCE:
;    result = avsigclip,  invvar, [threshold,sigmap, inmask=inmask, $
;             niteration=niteration]
;
; INPUTS:
;    array  -- 3d array to be averaged
;    invvar -- 3d array of inverse variance of array
;
; OPTIONAL INPUTS:
;    threshold -- threshold of chi^/d.f. before pixel is labeled "not
;                  OK" (default = 20)
;    sigmap -- threshold for + rejection (1 sigma default)
;	
; KEYWORDS:
;    niteration - number of iterations (default 3)
;    inmask -- input mask, dimensions as array, 1 on good pixels, 0 bad
;
; OUTPUTS:
;    result -- structure containing:
;    flux: 2d array, invvar averaged across the 3d input, with rejection
;    ivar: 2d invvar array, the sum of the invvar contributing to the pixel
;    CR: map of pixels rejected in at least one frame (0: all good,
;    bit on indicating which frame had rejection, up to 7)
;
; EXAMPLES:
;
; COMMENTS:
;   intended to eliminate CR effects within 3d data sets.  No spatial
;   averaging, only temporal averaging
;
; REVISION HISTORY:
;   MD 07jun02
;   DPF 21-Jul-2002  modified to loop over bad pixels only and reject
;                       the most deviant one. 
;   DPF 24-Jul-2002  force variance to be positive (roundoff error) -DPF
;   DPF 29-Jul-2002  niteration defaults to 3; compare to weighted av
;                         when rejecting most deviant value
;----------------------------------------------------------------------
  if n_elements(threshold) eq 0 then threshold = 20.0
  if n_elements(sigmap) eq 0 then sigmap = 1.0
  
;  threshold = 20.
  ncol = (size(array, /dimen))[0]
  nrow = (size(array, /dimen))[1]
  if n_elements(size(array, /dimen)) GT 2 then $
    nz = (size(array, /dimen))[2] $
    else nz = 0
  if n_elements(inmask) ne n_elements(array) then $
    keep = fix(array*0 +1) $     ;set mask
  else keep=inmask
;  deviation = array*0.
;  cut = 1. 

  if n_elements(array) eq float(ncol)*nrow then begin
     result =  {flux: array,  ivar: invvar,  mask: array*0}
     return, result
  endif


;  if NOT keyword_set(inmask) then inmask = bytarr(ncol, nrow)
  if NOT keyword_set(niteration) then niteration = 3

  for iter=0, niteration -1 do begin
    ztot = total(keep, 3)
    s_invvar = total(invvar*keep, 3) ;summed inverse variance
    invvar_m = s_invvar/(ztot + float(ztot eq 0)) ;mean inverse variance
    s_invvar_w = s_invvar + float(s_invvar eq 0)  ;to prevent NaN
;expected mean variance, input
    var_m = 1./(invvar_m + 1E-6*(invvar_m EQ 0))
    avg = total(array*keep*invvar, 3)/s_invvar_w
    avg2 = total(array^2*keep*invvar, 3)/s_invvar_w
    var = (avg2-avg^2) > 0 ;variance of 2-d array
    rms = sqrt(var)
    OK = (var lt threshold*var_m)*(invvar_m NE 0.) 
; simply want to know if the variance is too high in each pixel
    OK = (var*ztot lt 25*var_m) AND (invvar_m NE 0.)
    bad = where(OK eq 0, nbad)
    if iter lt niteration -1 then begin ;don't do on last iteration
       for i=0L, nbad-1 do begin ; toss most deviant pixel
          bx = bad[i] mod ncol
          by = bad[i]  /  ncol
;          av = total(array[bx,by,*]*keep[bx,by,*])/(total(keep[bx, by, *])>1)
          junk = max(abs(array[bx, by, *]-avg[bx, by])*keep[bx, by, *], ind)
          keep[bx, by, ind] = 0B
       endfor
    endif 
  endfor
 
  mask_3d = fix(1-keep) ;convert, turn on for bad pixel, off for good
  mask_2d = fix(avg*0)

  for i=0, nz-1 < 14 do $ ;no more than 14 frames preserve int mask
    mask_2d = mask_2d + mask_3d[*, *, i]*2^i  ;separate bit for each frame 

  bad = where(finite(avg) eq 0, nbad)
  if nbad gt 0 then avg[bad] =  -100. ;bad pixels
  
  result = {flux: avg, ivar: s_invvar, mask: mask_2d}
  if nbad gt 0 then result.ivar[bad] = 0 ; bad pixels

  return, result
end





