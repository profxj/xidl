pro setcrds, strct, FITS=fits

;  This program sets the structure cards for April 17, 2001 WFCCD data

  if not keyword_set( FITS ) then fits = 'wfccd_17apr01.fits'

  print, 'Setting the cards for 17apr01 WFCCD data'

; Translate

  nimg = n_elements(strct)

; YSTRT
  strct.ystrt = 400L

; Names
  strct[x_fmidx(strct,113):x_fmidx(strct,121)].Obj = 'HE1029-14'
  strct[x_fmidx(strct,122):x_fmidx(strct,142)].Obj = 'PG1116+22'
  strct[x_fmidx(strct,143):x_fmidx(strct,153)].Obj = '3C273'
  strct[x_fmidx(strct,154):x_fmidx(strct,163)].Obj = 'PKS1302-102'
  strct[x_fmidx(strct,164):x_fmidx(strct,172)].Obj = 'MRK1383'
  strct[x_fmidx(strct,173):x_fmidx(strct,182)].Obj = 'RJ1556+11'


; Afternoon 'tests'
  strct[x_fmidx(strct,lindgen(9)+84)].flg_anly = 0

; Bias
  strct[x_fmidx(strct,lindgen(20)+93)].type = 'ZRO'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; HE1029-14 S0

  strct[x_fmidx(strct,113)].type = 'FLT'
  strct[x_fmidx(strct,113)].masknm = 'S0'
  strct[x_fmidx(strct,114)].type = 'FLT'
  strct[x_fmidx(strct,114)].masknm = 'S0'
  strct[x_fmidx(strct,115)].type = 'ARC'
  strct[x_fmidx(strct,115)].masknm = 'S0'
  
  strct[x_fmidx(strct,116)].flg_anly = 0
  strct[x_fmidx(strct,117)].flg_anly = 0
  strct[x_fmidx(strct,118)].flg_anly = 0

  strct[x_fmidx(strct,lindgen(3)+113)].mask_id = 0
  strct[x_fmidx(strct,lindgen(3)+119)].mask_id = 0
  strct[x_fmidx(strct,119+lindgen(3))].masknm = 'S0'
  strct[x_fmidx(strct,121)].type = 'ARC'

  ; Mask output file
  strct[x_fmidx(strct,lindgen(9)+113)].msk_fil = 'Masks/H1029.small.0.stdout'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PG1116+22 B0
;    Calibs
  strct[x_fmidx(strct,122)].type = 'ARC'
  strct[x_fmidx(strct,122)].masknm = 'B0'
  strct[x_fmidx(strct,123)].type = 'FLT'
  strct[x_fmidx(strct,123)].masknm = 'B0'
  strct[x_fmidx(strct,124)].type = 'FLT'
  strct[x_fmidx(strct,124)].masknm = 'B0'

; Mask ID
  strct[x_fmidx(strct,lindgen(3)+122)].mask_id = 1
  strct[x_fmidx(strct,lindgen(3)+128)].mask_id = 1

;    Alignment
  strct[x_fmidx(strct,125+lindgen(3))].flg_anly = 0
;    Science
  strct[x_fmidx(strct,[128,129])].masknm = 'B0'
;    Calibs
  strct[x_fmidx(strct,130)].type = 'ARC'
  strct[x_fmidx(strct,130)].masknm = 'B0'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PG1116+22 B1 
;    Calibs
  strct[x_fmidx(strct,131)].type = 'ARC'
  strct[x_fmidx(strct,131)].masknm = 'B1'
;    Alignment
  strct[x_fmidx(strct,132)].flg_anly = 0
  strct[x_fmidx(strct,133)].flg_anly = 0
  strct[x_fmidx(strct,134)].flg_anly = 0
;    Science
  strct[x_fmidx(strct,135)].masknm = 'B1'
;    Alignment
  strct[x_fmidx(strct,136)].flg_anly = 0
  strct[x_fmidx(strct,137)].flg_anly = 0
  strct[x_fmidx(strct,138)].flg_anly = 0
;    Science
  strct[x_fmidx(strct,139)].masknm = 'B1'
;    Calibs
  strct[x_fmidx(strct,140)].type = 'ARC'
  strct[x_fmidx(strct,140)].masknm = 'B1'
  strct[x_fmidx(strct,141)].type = 'FLT'
  strct[x_fmidx(strct,141)].masknm = 'B1'
  strct[x_fmidx(strct,142)].type = 'FLT'
  strct[x_fmidx(strct,142)].masknm = 'B1'
; Mask ID
  strct[x_fmidx(strct,131)].mask_id = 2
  strct[x_fmidx(strct,135)].mask_id = 2
  strct[x_fmidx(strct,139+lindgen(4))].mask_id = 2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 3C273 B1 
;    Calibs
  strct[x_fmidx(strct,143)].type = 'FLT'
  strct[x_fmidx(strct,143)].masknm = 'B1'
  strct[x_fmidx(strct,146)].type = 'FLT'
  strct[x_fmidx(strct,146)].masknm = 'B1'
  strct[x_fmidx(strct,147)].type = 'ARC'
  strct[x_fmidx(strct,147)].masknm = 'B1'

;    Alignment
  strct[x_fmidx(strct,145)].flg_anly = 0
  strct[x_fmidx(strct,148)].flg_anly = 0
  strct[x_fmidx(strct,149)].flg_anly = 0
  strct[x_fmidx(strct,150)].flg_anly = 0
;    Science
  strct[x_fmidx(strct,151)].masknm = 'B1'
  strct[x_fmidx(strct,152)].masknm = 'B1'
;    Calibs
  strct[x_fmidx(strct,153)].type = 'ARC'
  strct[x_fmidx(strct,153)].masknm = 'B1'
; Mask ID
  strct[x_fmidx(strct,143)].mask_id = 3
  strct[x_fmidx(strct,[146,147])].mask_id = 3
  strct[x_fmidx(strct,151+lindgen(3))].mask_id = 3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PKS1302-102 B0 
;    Calibs
  strct[x_fmidx(strct,155)].type = 'FLT'
  strct[x_fmidx(strct,155)].masknm = 'B0'
  strct[x_fmidx(strct,156)].type = 'FLT'
  strct[x_fmidx(strct,156)].masknm = 'B0'
  strct[x_fmidx(strct,157)].type = 'ARC'
  strct[x_fmidx(strct,157)].masknm = 'B0'
;    Alignment
  strct[x_fmidx(strct,154)].flg_anly = 0
  strct[x_fmidx(strct,158)].flg_anly = 0
  strct[x_fmidx(strct,159)].flg_anly = 0
  strct[x_fmidx(strct,160)].flg_anly = 0
;    Science
  strct[x_fmidx(strct,161)].masknm = 'B0'
  strct[x_fmidx(strct,162)].masknm = 'B0'
;    Calibs
  strct[x_fmidx(strct,163)].type = 'ARC'
  strct[x_fmidx(strct,163)].masknm = 'B0'
; Mask ID
  strct[x_fmidx(strct,163+lindgen(3))].mask_id = 4
  strct[x_fmidx(strct,155+lindgen(3))].mask_id = 4

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; MRK1383 S0 
;    Calibs
  strct[x_fmidx(strct,164)].type = 'ARC'
  strct[x_fmidx(strct,164)].masknm = 'S0'
;    Alignment
  strct[x_fmidx(strct,165)].flg_anly = 0
  strct[x_fmidx(strct,166)].flg_anly = 0
  strct[x_fmidx(strct,167)].flg_anly = 0
;    Science
  strct[x_fmidx(strct,168)].masknm = 'S0'
  strct[x_fmidx(strct,169)].masknm = 'S0'
;    Calibs
  strct[x_fmidx(strct,170)].type = 'ARC'
  strct[x_fmidx(strct,170)].masknm = 'S0'
  strct[x_fmidx(strct,171)].type = 'FLT'
  strct[x_fmidx(strct,171)].masknm = 'S0'
  strct[x_fmidx(strct,172)].type = 'FLT'
  strct[x_fmidx(strct,172)].masknm = 'S0'
; Mask ID
  strct[x_fmidx(strct,164)].mask_id = 5
  strct[x_fmidx(strct,168+lindgen(5))].mask_id = 5

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; RJ1556+11 T0 
;    Calibs
  strct[x_fmidx(strct,173)].type = 'FLT'
  strct[x_fmidx(strct,173)].masknm = 'T0'
  strct[x_fmidx(strct,174)].type = 'ARC'
  strct[x_fmidx(strct,174)].masknm = 'T0'
;    Alignment
  strct[x_fmidx(strct,175)].flg_anly = 0
  strct[x_fmidx(strct,176)].flg_anly = 0
  strct[x_fmidx(strct,177)].flg_anly = 0
;    Science
  strct[x_fmidx(strct,178)].masknm = 'T0'
  strct[x_fmidx(strct,179)].masknm = 'T0'
;    Calibs
  strct[x_fmidx(strct,180)].type = 'ARC'
  strct[x_fmidx(strct,180)].masknm = 'T0'
  strct[x_fmidx(strct,181)].type = 'FLT'
  strct[x_fmidx(strct,181)].masknm = 'T0'
  strct[x_fmidx(strct,182)].type = 'FLT'
  strct[x_fmidx(strct,182)].masknm = 'T0'
; Mask ID
  strct[x_fmidx(strct,[173,174])].mask_id = 6
  strct[x_fmidx(strct,178+lindgen(5))].mask_id = 6

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; OUTPUT
  write_wfccdstr, strct, /anonly, outfil='wfccd_17apr01.lst', FITS=fits

  print, 'All done!'
  print, '----------------------------------------'

return
end
