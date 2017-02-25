;+ 
; NAME:
; mike_combsens   
;   Version 1.1
;
; PURPOSE:
;    Combine sensitivity functions.  Pretty manual right now
;
; CALLING SEQUENCE:
; mike_combsens
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
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-May-2004 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_combsens 

;
;  if  N_params() LT 2  then begin 
;      print,'Syntax - ' + $
;        'mike_calibstd, mike, indx, HSTFIL=, FXFIL= [v1.0]'
;    return
;  endif 

;  Optional Keywords

  ;;  Blue side
  bfeb12 = getenv('MIKE_DIR')+'pro/Std/Archive/sens_b12feb04.fits'
  bfeige = getenv('MIKE_DIR')+'pro/Std/Archive/sens_feige110b.fits'
  outfil = getenv('MIKE_DIR')+'pro/Std/Archive/sens_blue2.fits'

  ;; Copy
  spawn, 'cp '+bfeige+' '+outfil

  ;; Add from feb12
  for i=30L,35 do begin
      sens = xmrdfits(bfeb12, i, /silent)
      mwrfits, sens, outfil
  endfor

  ;; Red
  rfeb12 = getenv('MIKE_DIR')+'pro/Std/Archive/sens_r12feb04.fits'
  rfeige = getenv('MIKE_DIR')+'pro/Std/Archive/sens_feige110r.fits'
  vrfeige = getenv('MIKE_DIR')+'pro/Std/Archive/sens_feige110vr.fits'
  outfil = getenv('MIKE_DIR')+'pro/Std/Archive/sens_red2.fits'

  ;; Copy

  spawn, 'cp '+rfeige+' '+outfil  ;; Order 76--48

  ;; Add from Feige110 vr  Order 47--28
  for i=17L,36 do begin
      sens = xmrdfits(vrfeige, i, /silent)
      mwrfits, sens, outfil
  endfor

  ;; Testing
  print, 'mike_combsens: Testing...'
  for i=1L, 49 do begin
      sens = xmrdfits(outfil, i, /silent)
      print, sens.ordr
  endfor

  return
end

  

