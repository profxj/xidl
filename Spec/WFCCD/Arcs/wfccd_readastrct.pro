;+ 
; NAME:
; wfccd_readastrct
;    Version 1.0
;
; PURPOSE:
;    Solves arc solutions for a given mask
;
; CALLING SEQUENCE:
;   
;   wfccd_readastrct, wfstrct, WFARC=
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
;   wfccd_readastrct, wfstrct
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   15-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro wfccd_readastrct, infil, wfarc

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'wfccd_readastrct, infil, wfarc [v1.0]'
    return
  endif 

;  Optional Keywords

;  Open up the file

  anon = mrdfits(infil, 1, head, /silent)

;  Create wfarc

  narc = n_elements(anon)
  tmp = { wfccdarcstr }
  wfarc = replicate(tmp, narc)

; Fill it up
  for i=0L,narc-1 do begin
      wfarc[i].flg_anly = anon[i].flg_anly
      if tag_exist( anon[i], 'CENT' ) EQ 1B then wfarc[i].cent = anon[i].cent
      wfarc[i].wave = anon[i].wave
      wfarc[i].spec = anon[i].spec
      wfarc[i].fit.nord = anon[i].nord
      wfarc[i].fit.nrm = anon[i].nrm
      wfarc[i].fit.ffit = ptr_new(anon[i].ffit[0:anon[i].nord])
      wfarc[i].fit.func = 'POLY'
      wfarc[i].fit.niter = anon[i].niter
      wfarc[i].fit.lsig = anon[i].lsig
      wfarc[i].fit.hsig = anon[i].hsig
      wfarc[i].fit.flg_rej = anon[i].flg_rej
      wfarc[i].fit.maxrej = anon[i].maxrej
      wfarc[i].fit.minpt = anon[i].minpt
      wfarc[i].fit.rms = anon[i].rms
  endfor

  delvarx, anon, tmp

  return
end
