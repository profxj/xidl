;+ 
; NAME:
; x_allvelo
;
; PURPOSE:
;  Create an array of velocity arrays for a string of transitions
;
; CALLING SEQUENCE:
;   
;   all_velo = x_allvelo(wave, zabs, wrest, vmnx, ALL_PMNX=)
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   all_velo = x_allvelo(wave, zabs, wrest, vmnx, ALL_PMNX=)
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   24-Oct-2002 Written by JXP
;-
;------------------------------------------------------------------------------

function x_allvelo, wave, zabs, wrest, vmnx, ALL_PMNX=all_pmnx, NPIX=npix

  if (N_params() LT 3) then begin 
    print,'Syntax - ' + $
             'all_velo = x_allvelo(wave, zabs, wrest, vmnx, ALL_PMNX='
    print, '              ) [v1.0]'
    return, -1
  endif 

  nwv = n_elements(wave)
  ntrans = n_elements(wrest)
; Optional keywords

  if arg_present(ALL_PMNX) AND keyword_set(VMNX) then all_pmnx = lonarr(3,ntrans)
  if not keyword_set( NPIX ) then npix = nwv

  spl=2.9979d5

;;; LOOP
  fin_array = dblarr(npix, ntrans)

  for i=0L, ntrans-1 do begin
      ;; Create VELO
      velo = (wave-wrest[i]*(1.d + zabs))*spl/( wrest[i]*(1.d + zabs) )
      ;; PIXMIN, PIXMAX
      if arg_present(ALL_PMNX) AND keyword_set(VMNX) then begin
          mn = min(abs(velo-vmnx[0]),pmn)
          all_pmnx[0,i] = pmn
          mn = min(abs(velo-vmnx[1]),pmx)
          all_pmnx[1,i] = pmx
          all_pmnx[2,i] = pmx-pmn
          fin_array[0:all_pmnx[2,i],i] = velo[pmn:pmx]
      endif else begin
          mxpx = (npix-1)<(nwv-1)
          fin_array[0:mxpx,i] = velo[0:mxpx]
      endelse
  endfor

  return, fin_array
end

