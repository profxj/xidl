;+ 
; NAME:
; x_specrebin
;   Version 1.0
;
; PURPOSE:
;    Rebin a single data set to a new wavlength scale
;      Simple adding (no weighting by S/N)
;
; CALLING SEQUENCE:
;   
;   x_specrebin, gdpix, orig_wv, orig_fx, newwv, newfx
;
; INPUTS:
;   orig_wv
;   orig_fx
;   newwv
;
; RETURNS:
;
; OUTPUTS:
;   newfx
;
; OPTIONAL KEYWORDS:
;  VAR         
;
; OPTIONAL OUTPUTS:
;  NEWVAR      
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

pro x_specrebin, wv, fx, nwwv, nwfx, VAR=var, NWVAR=nwvar, $
                   SILENT=silent, REDBLUE=redblue

;
  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
             'x_specrebin, wv, fx, nwwv, nwfx, VAR=, NWVAR=, GDPIX= [V1.0]'
    return
  endif 


;  Optional Keywords
  npix = n_elements(wv)
  if not keyword_set(VAR) then var = fltarr(npix) + 1.

  nwpix = n_elements(nwwv)

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

; Calculate endpoints of the final array
  bwv = dblarr(nwpix+1)
  bwv[0:nwpix-1] = (nwwv + shift(nwwv,1))/2.
  bwv[0] = nwwv[0] - (nwwv[1]-nwwv[0])/2.
  bwv[nwpix] = nwwv[nwpix-1] + (nwwv[nwpix-1]-nwwv[nwpix-2])/2.

; Create tmp arrays for final array

  nwfx  = fltarr(nwpix)
  nwvar = fltarr(nwpix)

;  Loop!

  if not keyword_set(SILENT) then print, 'x_specrebin: Big Loop!'
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

      ; Now Sum up
      for kk=j1,j2 do begin
          ; Rejected pixels do not get added in
          if var[q] GT 0. then begin
              frac = ( (wvh[q] < bwv[kk+1]) - (wvl[q] > bwv[kk]) ) / (wvh[q]-wvl[q])
              nwfx[kk] = nwfx[kk] + frac * fx[q]
          endif
          ; Variance
          if arg_present( NWVAR ) then begin
              if var[q] LE 0. OR nwvar[kk] EQ -1 then nwvar[kk] = -1 else $
                nwvar[kk] = nwvar[kk] + frac * var[q]
          endif
      endfor
  endfor

  return
end
