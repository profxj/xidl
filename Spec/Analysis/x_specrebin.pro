;+ 
; NAME:
; x_specrebin
;   Version 1.1
;
; PURPOSE:
;    Rebin a single data set to a new wavlength scale
;      Simple linear interpolation.
;
; CALLING SEQUENCE:
;   x_specrebin, wv, fx, nwwv, nwfx, VAR=, NWVAR=, /SILENT, /REDBLUE
;
; INPUTS:
;   wv -- Original wavelength array
;   fx -- Orignal flux array
;   newwv -- New (desired) wavelength array
;
; RETURNS:
;   nwfx  -- New flux array
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  VAR= -- Original variance array
;  /REDBLUE --  Data runs from red to blue (not blue to red)
;  /SILENT
;  SMOOTH= -- Number of pixels to smooth over
;  /PRESERVE -- An old and wrong mode to preserve flambda [AVOID]
;  /FLAMBDA -- Preserve the flux (i.e. scale and divide by delta lambda)
;
; OPTIONAL OUTPUTS:
;  NWVAR= -- New variance array
;
; COMMENTS:
;
; EXAMPLES:
;   x_specrebin, gdpix, orig_wv, orig_fx, newwv, newfx
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   12-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_specrebin, wv, in_fx, nwwv, nwfx, VAR=in_var, NWVAR=nwvar, $
                 SILENT=silent, REDBLUE=redblue, CHK=chk, $
                 SMOOTH=SMOOTH, PRESERVE=preserve, FLAMBDA=flambda

;
  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
          'x_specrebin, wv, fx, nwwv, nwfx, ' + $
          'VAR=, NWVAR=, /CHK, /REDBLUE, SMOOTH=, /PRESERVE, /FLAMBDA [V1.2]'
    return
  endif 


;  Optional Keywords
  npix = n_elements(wv)
  if not keyword_set(IN_VAR) then in_var = fltarr(npix) + 1.
  nwpix = n_elements(nwwv)

; Smooth
  if keyword_set(SMOOTH) then begin
      kernel = gauss_kernel(smooth)
      fx = convol(in_fx, kernel)
      var = convol(in_var, kernel)
  endif else begin
      fx = in_fx
      var = in_var
  endelse


; Calculate wavelength endpoints for each pixel

  ; wvl = Lower endpoint
  ; wvh = Upper endpoint
  if keyword_set(REDBLUE) then begin  ; Wavelengths go red to blue on the CCD
      wvh = (wv + shift(wv, 1))/2.
      wvl = (wv + shift(wv, -1))/2.
      wvh[0] = wv[0] + (wv[0] - wv[1])/2.
      wvl[npix-1] = wv[npix-1] + (wv[npix-2] - wv[npix-1])/2.
  endif else begin
      wvl = (wv + shift(wv, 1))/2.
      wvh = (wv + shift(wv, -1))/2.
      wvl[0] = wv[0] - (wv[1] - wv[0])/2.
      wvh[npix-1,*] = wv[npix-1] + (wv[npix-1,*] - wv[npix-2,*])/2.
  endelse

  if keyword_set(FLAMBDA) then $
     scale = wvh-wvl $ ; dlambda 
  else scale = replicate(1., npix)

; Calculate endpoints of the final array
  bwv = dblarr(nwpix+1)
  bwv[0:nwpix-1] = (nwwv + shift(nwwv,1))/2.
  bwv[0] = nwwv[0] - (nwwv[1]-nwwv[0])/2.
  bwv[nwpix] = nwwv[nwpix-1] + (nwwv[nwpix-1]-nwwv[nwpix-2])/2.

; Create tmp arrays for final array

  nwfx  = fltarr(nwpix)
  tot_dwv  = fltarr(nwpix)
  nwvar = fltarr(nwpix)

;  Loop!

;  if not keyword_set(SILENT) then print, 'x_specrebin: Big Loop!'
  ; Loop on pixels
  for q=0L, npix-1 do begin

      ; Fill up the Final array
      ; No overlap
      if wvh[q] LE bwv[0] OR wvl[q] GE bwv[nwpix] then continue

      ; Find pixel that bw is within
      if wvl[q] LT bwv[0] then i1 = 0 else $
        i1 = where(wvl[q] LE shift(bwv,-1) AND wvl[q] GT bwv)
      ; Same for hw
      if wvh[q] GT bwv[nwpix] then i2 = nwpix-1 else $
        i2 = where(wvh[q] LE shift(bwv,-1) AND wvh[q] GT bwv)
      j1 = i1[0]
      j2 = i2[0]

      if j1 eq -1 or j2 eq -1 then continue ;added by KLC

      ; Now Sum up
      for kk=j1,j2 do begin
          ; Rejected pixels do not get added in
          if var[q] GT 0. then begin
              frac = ( (wvh[q] < bwv[kk+1]) - (wvl[q] > bwv[kk]) ) / (wvh[q]-wvl[q])
              if (wvh[q]-wvl[q]) LT 1e-3 then stop
              nwfx[kk] = nwfx[kk] + frac * fx[q] * scale[q]
              tot_dwv[kk] = tot_dwv[kk] + frac*scale[q]
          endif
          ; Variance
          if arg_present( NWVAR ) then begin
              if var[q] LE 0. OR nwvar[kk] EQ -1 then nwvar[kk] = -1 else $
                nwvar[kk] = nwvar[kk] + frac * var[q] * scale[q]^2
          endif
      endfor
   endfor

  if keyword_set(PRESERVE) or keyword_set(FLAMBDA) then begin
     dwv = bwv-shift(bwv,1) 
     ;dwv[0] = dwv[1]
;     printcol, dwv[1:nwpix], tot_dwv, nwwv, nwfx
     nwfx = nwfx/dwv[1:nwpix]
;     x_splot, wv, fx, xtwo=nwwv, ytwo=nwfx, ythr=sqrt(nwvar>0), /blo
  endif

  if keyword_set(CHK) then $
      x_splot, wv, fx, xtwo=nwwv, ytwo=nwfx, ythr=sqrt(nwvar>0), /blo

  return
end
