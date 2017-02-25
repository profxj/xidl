;+ 
; NAME:
;  tvstruct
;   Version 1.1
;
; PURPOSE:
;  Creates (defines) a structure for use in image GUIs
;
; CALLING SEQUENCE:
;   
;   tmp = {tvstruct}
;
; EXAMPLES:
; 
; REVISION HISTORY:
;
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro tvstruct__define

;  This routine defines the structure for direct images

  tmp = {tvstruct, $
         pos: fltarr(4),     $    ; Maps image onto entire window  
         svxymnx: fltarr(4), $
         xymnx: fltarr(4),   $       ; Mapped values of entire draw window
         winsize: lonarr(2), $       ; Number of draw window pixels
         size: lonarr(2), $          ; Number of draw window pixels
         gridsize: lonarr(2), $      ; Number of pixels for image
         xcurs: 0., $                ; x-coordinate of mouse
         ycurs: 0.  $                ; y-coordinate of mouse
         }

end
  
         
