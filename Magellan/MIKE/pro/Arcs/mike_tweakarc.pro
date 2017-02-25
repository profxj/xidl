;+ 
; NAME:
; mike_tweakarc   
;     Version 1.1
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   
;  mike_tweakarc, twkfil, 
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;  IDL fit file (one per order)  (e.g. Arcs/ArcECH_##fit.idl)
;
; OPTIONAL KEYWORDS:
;   /PINTER   - Perform fit for pre-identified lines
;   /INTER    - Identify lines interactively and then fit
;   LINLIST   -  Arc line list (default: $XIDL_DIR/Spec/Arcs/Lists/mike_thar.lst
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
;   mike_tweakarc, 'Arcs/Fits/mr1114_fit.idl', [72], /redblu
;   mike_tweakarc, 'Arcs/Fits/mr1114_fit.idl', [45],
;      'templ_arc_2x2R.idl', /redblue
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
pro mike_tweakarc, twkfil, ordrs, templfil, $
                   QAFIL=qafil, OSTR_FIL=ostr_fil, _EXTRA=extra


;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_tweakarc, twkfil, ordrs, [templfil], QAFIL= [v1.0]'
      return
  endif 
  
;;  Optional Keywords
  if not keyword_set( LINLIST ) then $
    linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/mike_thar.lst' 


  ;; Open templfil
  if keyword_set( TEMPLFIL ) then begin
      templfil=getenv('MIKE_DIR')+'/pro/Arcs/'+templfil
  endif

  x_tweakarc, twkfil, ordrs, templfil, LINLIST=linlist, $
              QAFIL=qafil, OSTR_FIL=ostr_fil, _EXTRA=extra


; All done
  print, 'mike_tweakarc: All done!'

  return
end
