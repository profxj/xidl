;+ 
; NAME:
; wfccd_writeastrct
;    Version 1.0
;
; PURPOSE:
;    Solves arc solutions for a given mask
;
; CALLING SEQUENCE:
;   
;   wfccd_writeastrct, wfstrct, WFARC=
;
; INPUTS:
;   wfstrct     - WFCCD structure
;
; RETURNS:
;
; OUTPUTS:
;   wfarc      -  WFCCD arc structure (fits file)
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_writeastrct, wfstrct
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   15-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro wfccd_writeastrct, wfarc, outfil

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'wfccd_writeastrct, wfarc, outfil [v1.0]'
    return
  endif 

;  Optional Keywords

;  Create anonymous structure

  narc = n_elements(wfarc)
  npix = n_elements(wfarc[0].wave)
  
  tmp = { $
          flg_anly: 0, $ ; Wavelength array
          cent: 0L, $ ; Row for fit
          wave: dblarr(npix), $ ; Wavelength array
          spec: fltarr(npix), $ ; Arc spectrum
          nord: 0L, $
          rms:  0.d, $
          ffit: dblarr(10), $
          lsig: 0., $
          hsig: 0., $
          niter: 0L, $
          minpt: 0L, $
          maxrej: 0L, $
          flg_rej: 0, $
          nrm:  dblarr(2) $
         }

  anon = replicate(tmp, narc)

; Fill it up
  for i=0L,narc-1 do begin
      ; Skip over bad slits
      if wfarc[i].flg_anly EQ 1 then begin
          anon[i].flg_anly = wfarc[i].flg_anly
          anon[i].cent = wfarc[i].cent
          anon[i].wave = wfarc[i].wave
          anon[i].spec = wfarc[i].spec
          anon[i].nord = wfarc[i].fit.nord
          anon[i].nrm = wfarc[i].fit.nrm
          anon[i].rms = wfarc[i].fit.rms
          anon[i].lsig = wfarc[i].fit.lsig
          anon[i].hsig = wfarc[i].fit.hsig
          anon[i].niter = wfarc[i].fit.niter
          anon[i].maxrej = wfarc[i].fit.maxrej
          anon[i].flg_rej = wfarc[i].fit.flg_rej
          anon[i].minpt = wfarc[i].fit.minpt
          anon[i].ffit[0:anon[i].nord] = *wfarc[i].fit.ffit
      endif
  endfor

; Write
  mwrfits, anon, outfil, /create

  delvarx, anon, tmp

  return
end
