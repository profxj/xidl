;+ 
; NAME:
; wfccd_fspeczhist
;    Version 1.0
;
; PURPOSE:
;   Converts SDSS eigenfunctions into WFCCD
;
; CALLING SEQUENCE:
;   
;   indx = wfccd_fspeczhist(wffspec, obj_nm)
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   XSIZE      - Size of gui in screen x-pixels (default = 1000)
;   YSIZE      - Size of gui in screen y-pixels (default = 600)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_fspeczhist, wfccd, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_fspeczhist, fspec_fil, BIN=bin, zmnx=zmnx, zqso=zqso, PSFIL=psfil, $
                      AONLY=aonly

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'wfccd_fspeczhist, fspec_fil, BIN=, ZMNX=, ZQSO=, PSFIL=, AONLY= [v1.0]'
    return
  endif 

; Optional Keywords
  if not keyword_set( BIN ) then bin = 0.005

; PSFIL

  if keyword_set( PSFIL ) then begin
      device, decompose=0
      !x.thick = 5
      !y.thick = 5
      !p.charthick = 4
      ps_open, file=psfil, font=1, /color, bpp=8
;      !y.margin = [5,2]
  endif

; Grab obj
  wfccd_wrfspec, wffspec, fspec_fil, /read

  if keyword_set( AONLY ) then $
    gdobj = where(wffspec.flg_anly NE 0 AND strtrim(wffspec.obj_id,2) EQ 'a', ngd) $
  else gdobj = where(wffspec.flg_anly NE 0, ngd)

  ; PLOT
  clr = getcolor(/load)
  wset, 0
  plothist, wffspec[gdobj].zans.z, bin=bin, xrange=zmnx, $
    xtitle='!17z', ymargin=[5,0], background=clr.white, color=clr.black, $
    /fill, fcolor=clr.green, ystyle=1, charsize=1.5

  if keyword_set( ZQSO ) then oplot, [zqso, zqso], $
    [0., 1e5], color=clr.red, linestyle=2

; CLOSE PSFIL

  if keyword_set(PSFIL) then begin
      ps_close, /noprint, /noid
      device, decomposed=1
      !x.thick = 1
      !y.thick = 1
      !p.charthick = 1
;      !y.margin = [4,2]
  endif

  return
end
