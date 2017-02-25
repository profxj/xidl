;+ 
; NAME:
; esi_tweakarc   
;     Version 1.1
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   
;  esi_fitarc, esi, setup, obj_id, [chip], /INTER, LINLIST=, /CHK, /CLOBBER,
;  SIGREJ=, /DEBUG, IORDR=, /PINTER 
;
; INPUTS:
;   esi     -  HIRES structure
;   setup    -  Integer defining setup
;   obj_id   -  Object identifier
;   [chip]   -  Blue (1), Red (2), or both [1,2L]    (Default: [1,2L])
;
; RETURNS:
;
; OUTPUTS:
;  IDL fit file (one per order)  (e.g. Arcs/ArcECH_##fit.idl)
;
; OPTIONAL KEYWORDS:
;   /PINTER   - Perform fit for pre-identified lines
;   /INTER    - Identify lines interactively and then fit
;   LINLIST   -  Arc line list (default: $XIDL_DIR/Spec/Arcs/Lists/esi_thar.lst
;   /CHK      - Manually check steps along the way
;   /DEBUG    - Debugging
;   SIGREJ=   - Rejection sigma for outliers in arc line fitting
;              (default: 2.)
;   IORDR     - Initial order for analysis
;   /CLOBBER  - Overwrite previous fits
;   SHFTPRM=  - Fit structure for shifting the orders of the arc
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_fitarc, esi, 1, 1
;
;
; PROCEDURES/FUNCTIONS CALLED:
;   x_fitarc
;
; REVISION HISTORY:
;   Summer-2005 Created by JXP 
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro esi_tweakarc, twkfil, ordrs, _EXTRA=extra, $
                    QAFIL=qafil, OSTR_FIL=ostr_fil, CUAR=cuar


;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_tweakarc, twkfil, ordrs, [templfil], QAFIL= [v1.0]'
      return
  endif 
  
;;  Optional Keywords
  if not keyword_set( LINLIST ) then begin
      if keyword_set( CUAR ) then $
        linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/ESI_CuArech.lst' $
      else linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/ESI_ech.lst'
  endif


  ;; Open templfil
  if not keyword_set(TEMPLFIL) then begin
      if keyword_set(CUAR) then $
        templfil = getenv('XIDL_DIR')+'/ESI/CALIBS/ECH_CuArarcfit.idl' $
      else templfil = getenv('XIDL_DIR')+'/ESI/CALIBS/ECH_arcfit.idl' 
  endif
      

  x_tweakarc, twkfil, ordrs, templfil, QAFIL=qafil, OSTR_FIL=ostr_fil, $
    GUESS_ORDR=lindgen(10), LINLIST=linlist, _EXTRA=extra

; All done
  print, 'esi_tweakarc: All done!'

  return
end
