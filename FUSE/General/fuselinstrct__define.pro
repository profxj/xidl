;+ 
; NAME:
; fuselinstrct__define
;  V1.1
;
; PURPOSE:
;    Defines the structure for FUSE absorption lines
;
; CALLING SEQUENCE:
;  tmp = {fuselinstrct}
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
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Nov-2003   Created by JXP
;-
;------------------------------------------------------------------------------
pro fuselinstrct__define

;  This routine defines the line list structure

  tmp = {fuselinstrct, $
         ion: ' ', $
         wrest: 0.d, $
         wv_lim: dblarr(2), $
         flg: 0, $         ; 1 = Analysis, 2 = Lower limit, 4 = Upper limit
         instr: 0, $        ; Instrument (see below)
         EW: fltarr(20), $           ; Colm (log) 
         sigEW: fltarr(20), $
         Ncolm: 0., $
         sigNcolm: 0., $
         zabs: 0.d,  $       ; Redshift
         zsig: 0.d  $       ; Error
         }

end
  
; Instrument 1=1bsic, 2=2asic, 4=2blif, 8=1alif, 16=1asic,
;     32=1blif, 64=2alif, 128=STIS, 256=?
         
