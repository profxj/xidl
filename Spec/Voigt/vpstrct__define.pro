;+ 
; NAME:
; vpstrct__define
;   Version 1.1
;
; PURPOSE:
;   creates structure to contain data from a vpfit fort.26 output file,
;   including flux file names, region start/stop, abs. lines, N, Nerr, etc.
;   For the lines, has both string and double representations of the data,
;   so information on 'tied' or 'locked' lines is not lost.
; 
;
; CALLING SEQUENCE:
;   tmp = {vpstrct}
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
;   Limits total number of regions to 500, total number of lines to 1000
;
; EXAMPLES:
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   Written by GEP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro vpstrct__define

  tmp = {vpstrct, $
         nreg: 0L, $
         fluxfil: strarr(500), $
         reg_beg: dblarr(500),        $   
         reg_end: dblarr(500),        $   
         nion: 0L, $
         ion: strarr(1000),      $    
         z: dblarr(1000), $         
         z_str: strarr(1000), $         
         zerr: dblarr(1000), $         
         n: dblarr(1000), $         
         nerr: dblarr(1000), $         
         n_str: strarr(1000), $         
         b: dblarr(1000), $         
         berr: dblarr(1000), $         
         b_str: strarr(1000) $ 
         }

end
  
         
