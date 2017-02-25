;+ 
; NAME:
; xsp1dstrct__define
;    Version 1.0
;
; PURPOSE:
;  This routine defines the structure for individual object spectra
;
; CALLING SEQUENCE:
;   tmp = {xsp1dstrct}
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
;   02-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------
pro xsp1dstrct__define

;  This routine defines the structure for individual object spectra

   tmp = create_struct( $
    name = 'xsp1dstrct', $
    'class'      ,  '', $
    'subclass'   ,  '', $
    'z'          , 0.0, $
    'z_err'      , 0.0, $
    'rchi2'      , 0.0, $
    'dof'        ,  0L, $
    'rchi2diff'  , 0.0, $
    'tfile'      ,  '', $
    'tcolumn'    , lonarr(10) - 1L, $
    'npoly'      ,  0L, $
    'theta'      , fltarr(10), $
    'vdisp'      , 0.0, $
    'vdisp_err'  , 0.0  $
   )

   return
end
