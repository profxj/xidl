;+ 
; NAME:
; apfkey__define
;     Version 1.1
;
; PURPOSE:
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
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   24-Jan-2005 Written by JXP
;-
;------------------------------------------------------------------------------
pro apfkey__define

;  This routine defines the structure for HIRES keywords
;  [v1.1]

  tmp = {apfkey, $
            frmcrd: ' ', $
            utcrd: ' ', $
            ccdspeed: ' ', $
            ccdgain: ' ', $
            racrd: ' ', $
            deccrd: ' ', $
            expcrd: ' ', $
            decker: ' ', $
            objcrd: ' ', $
            binning: ' ', $
            block: ' ', $
            eqxcrd: ' ', $
            xdispers: ' ', $
            xdangl: ' ', $
            echangl: ' ', $
            obstyp: ' ', $
            hatch: ' ', $
            pane: ' ', $
            lamp: ' ', $
            lampfil: ' ', $
            ampmod: ' ', $
            mosmod: ' ', $
            amcrd: ' ' $
        }
end
  
         
