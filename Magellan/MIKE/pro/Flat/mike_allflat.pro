;+ 
; NAME:
; mike_allflat   
;     Version 1.1
;
; PURPOSE:
;    Performs the 3 main flat steps:  mike_mkmflat, mike_mktflat,
;    mike_edgeflat
;
; CALLING SEQUENCE:
;   
;  mike_allflat, mike, setup, [side], /CLOBBER, _EXTRA=extra
;
; INPUTS:
;   mike     -  ESI structure
;   [slit]  -  Slit size (e.g. 0.5, 0.75, 1.0)
;
; RETURNS:
;
; OUTPUTS:
;  One normalized flat per slit width
;
; OPTIONAL KEYWORDS:
;   /CHK     - Check steps in mike_edgeflat
;   _EXTRA   - User may supply one of many extra keywords to
;              the various routines.  You should consult the specific
;              routines to decide if you wish to include any of these options.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_allflat, mike, 1, 
;
; PROCEDURES/FUNCTIONS CALLED:
;  mike_mkmflat
;  mike_mktflat
;  mike_edgeflat
;
; REVISION HISTORY:
;   14-Aug-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_allflat, mike, setup, side, CLOBBER=clobber, CHK=chk, _EXTRA=extra

;
 if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_allflat, mike, setup, [side], /CLOBBER, /TFCHK [v2.0]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( SIDE ) then side = [1L,2L]
  if not keyword_set( clobber ) then clobber = 0L
  if not keyword_set( iflat ) then iflat = 0L
 

  for i=0, n_elements(side) -1 do begin

    qq = side[i]
 
    ; Milky Flat
    mike_mkmflat, mike, setup, qq, CLOBBER=clobber, _EXTRA=EXTRA

   ; Trace flat
    ;; COMBINE
     mike_mktflat, mike, setup, qq, CLOBBER=clobber, _EXTRA=EXTRA

    ;; EDGEFLAT
    mike_edgeflat, mike, setup, qq, CLOBBER=clobber, CHK=chk, _EXTRA=EXTRA

  endfor

  ;; SLIT PROFILE  
  ;;  Call mike_sliflat manually after mike_allarc
  print, 'Call mike_slitflat manually after mike_allarc:: '
  print, 'IDL> mike_slitflat, mike, setup, CLOBBER=clobber'
  print, 'Call mike_slitflat manually after mike_allarc:: '

  print, 'mike_allflat: All done!'
  return
end
