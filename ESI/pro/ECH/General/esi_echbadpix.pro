;+ 
; NAME:
; esi_echbadpix   
;     Version 1.1
;
; PURPOSE:
;    Set bad pixels to 0
;
; CALLING SEQUENCE:
;   
;  esi_echbadpix, img
;
; INPUTS:
;   img     - 2D IMG array
;
; RETURNS:
;
; OUTPUTS:
;  Full image with badpixels suppressed
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Only set for 1x1 binning for now
;
; EXAMPLES:
;   esi_echbadpix, img
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   12-Aug-2002 Written by JXP
;   03-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echbadpix, esi, img

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'esi_echbadpix, img [v1.1]'
      return
  endif 
  
;  Optional Keywords
  cbin = esi.cbin
  rbin = esi.rbin

;  Set bad pix to 0
  case esi.cbin of
      1: begin
         case esi.rbin of 
             1: begin
                 ;; HOT SPOT
                 img[35:150,3850:3930] = 0.
                 img[151:209,3808:3880] = 0.
                 img[151:190,3881:3960] = 0.
                 ;; Bad CLMS
                 img[421:437,2647:4095] = 0.
                 img[875:905,3827:4095] = 0.
                 img[891:910,3814:3827] = 0.
             end
             else: stop
         endcase
     end
     4: begin
         case esi.rbin of 
             4: begin
                 mskfil = 'Bias/MSK_'+strtrim(cbin,2)+'x'+strtrim(rbin,2)+'.fits'
                 mskimg = xmrdfits(mskfil, 0, /silent)
                 msk = where(mskimg EQ 0B)
                 img[msk] = 0.
             end
             else: stop
         endcase
     end
     else: stop
 endcase
                 

  return
end
