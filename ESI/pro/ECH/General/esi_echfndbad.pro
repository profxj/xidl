;+ 
; NAME:
; esi_echfndbad   
;     Version 1.1
;
; PURPOSE:
;    Set bad pixels to 0
;
; CALLING SEQUENCE:
;   
;  esi_echfndbad, img
;
; INPUTS:
;   img     - 2D IMG array
;
; RETURNS:
;
; OUTPUTS:
;  Full image with fndbadels suppressed
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Only set for 1x1 binning for now
;
; EXAMPLES:
;   esi_echfndbad, img
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   23-Oct-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echfndbad, esi, cbin, rbin

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'esi_echfndbad, esi, cbin, rbin [v1.0]'
      return
  endif 
  
;  Optional Keywords

  gddrk = where(esi.mode EQ 2 AND esi.flg_anly NE 0 AND $
                esi.cbin EQ cbin AND esi.rbin EQ rbin AND $
                esi.type EQ 'DRK', ngddrk)

;  Bias subtract
  esi_subbias, esi, gddrk, /NOBAD
  xcombine, 'OV/ov_'+esi[gddrk].img_root, img_dark, head, $
    FCOMB=2, SCALE='MED', GAIN=esi[gddrk[0]].gain, $
    RN=esi[gddrk[0]].readno

  sz = size(img_dark, /dimensions)
; Stats
  djs_iterstat, img_dark[0:1024/cbin,2484/rbin:*], sigrej=2.5, $
    maxiter=9, mask=mask, sigma=fsig, median=fmed

; Median

  gdmsk = bytarr(sz[0],sz[1])
  gdmsk[*] = 1B
  gdmsk[0:1024/cbin,2484/rbin:*] = mask

  gdmsk = median(gdmsk, 3)

;  Write to Bias directory

  case cbin of 
      4: begin
          case rbin of 
              4: begin
                  gdmsk[106:110,661L:*] = 0B
                  gdmsk[2:2,940L:*] = 0B
              end
              else: stop
          endcase
      end
      else: stop
  endcase

  outfil = 'Bias/MSK_'+strtrim(cbin,2)+'x'+strtrim(rbin,2)+'.fits'
  mwrfits, gdmsk, outfil, /create

  print, 'esi_echfndbad: Writing ',  outfil

  return
end
