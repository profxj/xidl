;+
; NAME:
; kast_crays
;   Version 1.1
;
; PURPOSE:
;   Detect cosmic rays, set affected pixels to have var = -1
;
; CALLING SEQUENCE:
;   kast_crays, kast, setup, side, obj_id, [exp_id], CHK=chk
;
; INPUTS:
;   kast  --  Kast IDL structure
;  setup  --  Setup value
;   side  --  Specific camera [blue (1) vs. red (2)]
; obj_id  --  Object value
;  [exp_id]  --  Exposure indices
; 
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /CHK  -- Print out info on CRays
;
; COMMENTS:
;
; EXAMPLES:
;    kast_crays, kast, 0, 1, 2, /chk
;
; PROCEDURES/FUCTIONS CALLED:
;
; REVISION HISTORY:
;   01-September-2003 Written by GEP
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro kast_crays, kast, setup, side, obj_id, exp_id, CHK=chk

; 
  if N_params() LT 4 then begin
    print, 'Syntax - ', + $
      'kast_crays, kast, setup, side, obj_id, [exp_id], /CHK [v1.0]'
    return
  endif

; Get exposures
  allexp = where(kast.type EQ 'OBJ' AND kast.flg_anly NE 0 AND $
                 kast.setup EQ setup AND kast.side EQ side AND $
                 kast.obj_id EQ obj_id AND kast.mode EQ 1 AND $
                 kast.exp GT 0)
  if keyword_set(exp_id) then exp = allexp[exp_id] else exp=allexp
  nexp = n_elements(exp)

; Need more than one exposure
  if nexp LT 2 then begin
    print, 'Need at least two exposures to find cosmic rays!'
    return
  endif

; Obj Structure
  for i=0L,nexp-1 do begin
    tmp = xmrdfits(kast[exp[i]].obj_fil, 1, STRUCTYP='specobjstrct',/silent)
    if i EQ 0 then kastobj = tmp else kastobj = [kastobj, tmp]
  endfor
  allobj = where(kastobj.exp NE 0)

; Get flux, variance, scales
  npix = kastobj[allobj[0]].npix
  flux = findgen(nexp,npix)
  var = findgen(nexp,npix)
  scl = findgen(nexp)
  vscl = findgen(nexp)
  if nexp GT 2 then begin
    for i=0L,nexp-1 do begin
      flux[i,0:npix-1] = kastobj[allobj[i]].fx[0:npix-1]
      var[i,0:npix-1] = kastobj[allobj[i]].var[0:npix-1]
      scl[i,*] = kastobj[allobj[0]].exp/kastobj[allobj[i]].exp
      vscl[i,*] = (kastobj[allobj[0]].exp/kastobj[allobj[i]].exp)^2
      flux[i] = flux[i]*scl[i]
      vscl[i] = flux[i]*vscl[i]
    endfor
  endif else begin
    for i=0L,nexp-1 do begin
      scl[i,*] = kastobj[allobj[0]].exp/kastobj[allobj[i]].exp
      vscl[i,*] = (kastobj[allobj[0]].exp/kastobj[allobj[i]].exp)^2
    endfor
  endelse

; Normalize variance
    
; Check for cosmic rays
  if nexp EQ 2 then begin
    for i=0L,npix-1 do begin
      if kastobj[allobj[0]].fx[i] GT kastobj[allobj[1]].fx[i]*scl[1] + $
              10*sqrt(kastobj[allobj[0]].var[i]) then begin
         kastobj[allobj[0]].var[i] = -1.
         if keyword_set(CHK) then print, 'CR in obj 1, pix: ',i
      endif else begin
        if kastobj[allobj[1]].fx[i]*scl[1] GT kastobj[allobj[0]].fx[i] + $
              10*sqrt(kastobj[allobj[1]].var[i]*vscl[1]) then begin
          kastobj[allobj[1]].var[i] = -1.
          if keyword_set(CHK) then print, 'CR in obj 2, pix: ',i
        endif
      endelse
    endfor
  endif else begin
    for i=0L,npix-1 do begin
      fxmed = median(flux[allobj,i])
      for kk=0L,nexp-1 do begin
        if flux[kk,i] GT fxmed + 10*sqrt(var[kk,i]) then begin
            kastobj[allobj[kk]].var[i] = -1.
            if keyword_set(CHK) then print, 'CR in obj ',kk,', pix: ',i
        endif
      endfor
    endfor
  endelse

; Confirm CRs if /CHK
  if keyword_set(CHK) then begin
    if nexp EQ 2 then begin
      x_splot, kastobj[allobj[0]].fx, ytwo=scl[1]*kastobj[allobj[1]].fx, /block
    endif else begin
      for i=1L,nexp-1 do begin
        x_splot, flux[0,0:npix-1], ytwo=flux[i,0:npix-1], /block
      endfor
    endelse
  endif

  for i=0L,nexp-1 do begin
    print, 'kast_crays: Updating ', kast[exp[i]].obj_fil
    mwrfits, kastobj[allobj[i]], kast[exp[i]].obj_fil, /create
    spawn, 'gzip -f '+kast[exp[i]].obj_fil
  endfor

  print, 'kast_crays: All done!'

end
