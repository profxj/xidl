;+ 
; NAME:
; uves_tweakarc   
;     Version 1.1
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   
;  uves_fitarc, uves, setup, obj_id, [chip], /INTER, LINLIST=, /CHK, /CLOBBER,
;  SIGREJ=, /DEBUG, IORDR=, /PINTER 
;
; INPUTS:
;   uves     -  HIRES structure
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
;   LINLIST   -  Arc line list (default: $XIDL_DIR/Spec/Arcs/Lists/uves_thar.lst
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
;   uves_fitarc, uves, 1, 1
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
pro uves_tweakarc, twkfil, ordrs, templfil, _EXTRA=extra, $
                    QAFIL=qafil, OSTR_FIL=ostr_fil


;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'uves_tweakarc, twkfil, ordrs, [templfil], QAFIL= [v1.0]'
      return
  endif 
  
;;  Optional Keywords
  if not keyword_set( LINLIST ) then $
    linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/hires_thar.lst' 

  ;; Open line list
;  x_arclist, linlist, lines


  ;; Open templfil
  if keyword_set( TEMPLFIL ) then begin
      templfil = getenv('XIDL_DIR')+ $
                 '/VLT/UVES/pro/Arcs/Templates/'+templfil
      
;      tguess = guess_ordr
;      tfit = all_arcfit
;      tspec = sv_aspec
  endif

  x_tweakarc, twkfil, ordrs, templfil, QAFIL=qafil, OSTR_FIL=ostr_fil, $
    _EXTRA=extra

; All done
  print, 'uves_tweakarc: All done!'

  return
end
