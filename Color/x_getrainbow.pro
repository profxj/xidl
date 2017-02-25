;+ 
; NAME:
; x_getrainbow   
;   Version 1.1
;
; PURPOSE:
;    Passes back a structure defining a rainbow
;
; CALLING SEQUENCE:
;   colm = ew_to_colm([lambda], [EW], /RVRS)
;
; INPUTS:
;   lambda  - Rest Wavelength (Ang)
;   EW       - EW (mA) [or column density (linear)]
;
; RETURNS:
;   colm   - Column density (linear)
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /SILENT -- Suppress written output
;  /RVRS   -- Take a column density input and output the EW assuming
;             the linear COG
;  TOLER=  -- Tolerance on wavelength [default: 0.1]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   colm = ew_to_colm([1215.6701],[50.])
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-May-2008 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_getrainbow 

 clr = getcolor(/load)
 tmp = { name: 'Red', $
         color: clr.red, $
         wvmnx: dblarr(2) $
       }
 rainbow = replicate(tmp,7)

 ;; Red
 rainbow[6].name = 'Red'
 rainbow[6].color = clr.red
 rainbow[6].wvmnx = [6300., 6800]

 ;; Orange
 rainbow[5].name = 'Orange'
 rainbow[5].color = clr.orange
 rainbow[5].wvmnx = [5800., 6300]

 ;; Yellow
 rainbow[4].name = 'Yellow'
 rainbow[4].color = clr.yellow
 rainbow[4].wvmnx = [5300., 5800]

 ;; Green
 rainbow[3].name = 'Green'
 rainbow[3].color = clr.darkgreen
 rainbow[3].wvmnx = [4800., 5300]

 ;; Blue
 rainbow[2].name = 'Blue'
 rainbow[2].color = clr.blue
 rainbow[2].wvmnx = [4300., 4800]

 ;; Violet
 rainbow[1].name = 'Violet'
 rainbow[1].color = clr.violet
 rainbow[1].wvmnx = [3800., 4300]

 ;; Purple
 rainbow[0].name = 'Purple'
 rainbow[0].color = clr.purple
 rainbow[0].wvmnx = [3300., 3800]

 return, rainbow

end
     
