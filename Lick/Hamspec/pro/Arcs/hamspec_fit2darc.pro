;+ 
; NAME:
; hamspec_fit2darc   
;     Version 1.1
;
; PURPOSE:
;  To fit the arc lines identified in hamspec_fitarc as a fucntion of
;  their y-centroid and order number.  The main routine is in
;  hamspec_fit2darc_work.  The fit is a simple least-squares with one
;  round of rejection.
;
; CALLING SEQUENCE:
;   
;  hamspec_fit2darc, hamspec, setup, [obj_id, chip]
;
; INPUTS:
;   hamspec    -  HIRES structure
;   setup    -  Integer defining setup
;   [obj_id] -  Object identifier
;   [chip] -  Blue (1), Green (2), Red (3), or multiple [1,2L]
;            (Default: [1,2,3L])
;
; RETURNS:
;
; OUTPUTS:
; A fits file containing the 2D solution.  Named something like
; 'Arcs/Fits/Arc_B0539_fit2D.fits' 
;
; OPTIONAL KEYWORDS:
;   NOCOEFF   - Number of coefficients to use in the x-direction
;               (default: 4)
;   NYCOEFF   - Number of coefficients to use in the y-direction
;               (default: 4)
;   /CLOBBER  - Overwrite any previous solution
;   /DEBUG    - debug
;   /CHKRES   - Plot the residuals
;   ARCFIL=   - Name of the arc file to process (Optional to using
;               setup, chip, etc.)
;   SIGREJ=   - Sigma for rejectin [default: 6]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hamspec_fit2darc, hamspec, 1.0
;
; PROCEDURES/FUNCTIONS CALLED:
;  hamspec_fit2darc_work
;
; REVISION HISTORY:
;   28-Feb-2003 Written by SB
;   18-Apr-2003 Revised by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hamspec_fit2darc, hamspec, setup, CLOBBER=clobber, $
                    nycoeff=nycoeff, nocoeff=nocoeff, DEBUG=debug, $
                    CHKRES=chkres, ARCFIL=arcfil, SIGREJ=sigrej
;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hamspec_fit2darc, hamspec, setup, nycoeff=, nocoeff=' 
      print, '     /CLOBBER, /CHKRES, /DEBUG, ARCFIL=, SIGREJ= [v1.1]'
      return
  endif

  if NOT keyword_set(nycoeff) then nycoeff = 4
  if NOT keyword_set(SIGREJ) then sigrej = 6.

  if keyword_set(nocoeff) then flg_ocoeff = 1L else flg_ocoeff = 0L

; Loop on chip
  ;; Nocoeff
  if flg_ocoeff NE 1 then begin
     ;; This is only good for the UV cross disperser for now
     nocoeff = 5
  endif

  ;; Grab all obj indices
  if not keyword_set(ARCFIL) then begin
     arcfil = hamspec_getfil('arc_fil', setup, /name)
  endif 

  ;; Dewar specific
  idx = where(hamspec.flg_anly NE 0 AND $
               hamspec.setup EQ setup) 
  if strtrim(hamspec[idx[0]].ccd,2) EQ 'e2vCCD203-824kx4kthin' then begin
     sigrej = 3.
  endif
  
  ;; ORD_STR
  ordr_str = hamspec_getfil('ordr_str', setup)

  ;; Outfil
  pos = strpos(arcfil, '.fits')
  frame = long(strmid(arcfil,pos-4,4))
  
  ;;  Check for outfil
  outfil = hamspec_getfil('arc_2Dfit', setup, $
                          /name, CHKFIL=chkf,  FRAME=frame)
  if chkf NE 0 and not keyword_set( CLOBBER ) then begin
     print, 'hamspec_fit2darc: Arc fit file exists. ' + $
            'Continuing..'
     return
  endif
     
  ;;  QAFIL
  qafil = hamspec_getfil('qa_arc2dfit', setup, $
                         /name, CHKFIL=chkf,  FRAME=frame)
  ;; Open arc 1D fit
  arc_info = hamspec_getfil('arc_fit', setup, $
                            /name, CHKFIL=chkfil,  FRAME=frame)
  ;;
  rslt = x_fit2darc( arcfil, ordr_str, arc_info, SIGREJ=sigrej, $
                     nycoeff=nycoeff, CHKRES=chkres, $
                     nocoeff=nocoeff, OUT_STR=out_str, $
                     CLOBBER=clobber, DEBUG=debug, QAFIL=qafil)
  ;;
  mwrfits, out_str, outfil, /create
  ;; 
  if size(rslt,/tname) NE 'STRING' then stop
      
  print, 'hamspec_fit2darc: All done'
  return

end

