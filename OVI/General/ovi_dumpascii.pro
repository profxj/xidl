;+ 
; NAME:
; ovi_dumpascii   
;     Version 1.1
;
; PURPOSE:
;    Creates a Summary ASCII file of the spectroscopic data
;
; CALLING SEQUENCE:
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
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   25-Jan-2008 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro ovi_dumpascii, struct, outfil, FLGANLY=flganly 

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'ovi_dumpascii, struct, outfil, FLGANLY=flganly [v1.0]'
      return
  endif 

  if not keyword_set(FLGANLY) then FLGANLY = 0L
  ;; 

  ;; Read in
  ovistr = xmrdfits(struct,1)

  ;; Parse
  gd = where(ovistr.flg_anly GE FLGANLY)
  ovistr = ovistr[gd]

  ;; Writeout
  
  writecol, outfil, $
            ovistr.id, $
            ovistr.obj_id, $
            ovistr.xypix[0], $
            ovistr.xypix[1], $
            ovistr.mag[0], $
            ovistr.mag[1], $
            ovistr.flg_anly, $
            ovistr.z, $
            ovistr.gal_coeff[0], $
            ovistr.gal_coeff[1], $
            ovistr.gal_coeff[4], $
            FMT='(i5,1x,a1,1x,f7.2,1x,f7.2,1x,f5.2,1x,f5.2,1x,i1,1x,f8.5,'+$
            '1x,f5.2,1x,f5.2,1x,f5.2)'
  
; All done
  print, 'ovi_dumpascii: All done!'

end
