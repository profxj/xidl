;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; civ_instantstrct.pro               
; Author: Kathy Cooksey                      Date: 20 Feb 2008
; Project: HST Metal-line System Survey (CIV Focus) with
;          Xavier Prochaska
; Description: Assign values to every element of structure
;              civcandstrct
; Input: 
;   num - integer number of civcandstrct to create
; Optional:
;   synth - add 'overdensity' and 'temp' tags to structure
; Output: 
;   returns (array) of instantiated civcandstrct
; Example:
;   civcand = civ_instantstrct(10)
; History:
;   20 Feb 2008  created by KLC
;    7 Apr 2009  added synth option, KLC
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function civ_instantstrct,num, synth=synth

  if keyword_set(synth) then begin
     tmp2 = civ_instantstrct(1)
     arr = replicate(0.d,30,3)
     tmp = create_struct(tmp2,'overdensity',arr,'temp',arr) ; name= doesn't work
  endif else tmp = { civcandstrct }
  
  if num gt 1 then civcand = replicate(tmp,num) $
  else civcand = tmp

  ;; Scalars (9)
  civcand.qso = ' '             ;string
  civcand.ra = 0.               ;double
  civcand.dec = 0.              ;double
  civcand.zqso = 0.             ;double
  civcand.search_fil = ' '      ;string
  civcand.instr_fil = ' '       ;string
  civcand.rating_eye = 0        ;int
  civcand.other_comments = ' '  ;string


  ;; 1-D arrays (12)
  civcand.comment_eye = 0b      ;byte array
  civcand.flg_sys = 0           ;integer
  civcand.ion = ' '             ;string
  civcand.wrest = 0.            ;double
  civcand.flg_colm = 0          ;integer
  civcand.instr = 0             ;integer
  civcand.ew = 0.               ;double
  civcand.sigew = 0.            ;double
  civcand.ncolm = 0.            ;double
  civcand.signcolm = 0.         ;double
  civcand.b = 0.                ;double
  civcand.sigb = 0.             ;double
  civcand.zabs = 0.             ;double
  civcand.zsig = 0.             ;double

  ;; 2-D arrays (4)
  civcand.aodm_vel = 0.         ;double
  civcand.aodm_civ = 0.         ;double
  civcand.sigaodm_civ = 0.      ;double
  civcand.wv_lim=0.             ;double

  return,civcand

end
