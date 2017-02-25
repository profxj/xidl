;+ 
; NAME:
; apf_fit2darc   
;     Version 1.1
;
; PURPOSE:
;  To fit the arc lines identified in apf_fitarc as a fucntion of
;  their y-centroid and order number.  The main routine is in
;  apf_fit2darc_work.  The fit is a simple least-squares with one
;  round of rejection.
;
; CALLING SEQUENCE:
;   
;  apf_fit2darc, apf, setup, [obj_id, chip]
;
; INPUTS:
;   apf    -  HIRES structure
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
;   apf_fit2darc, apf, 1.0
;
; PROCEDURES/FUNCTIONS CALLED:
;  apf_fit2darc_work
;
; REVISION HISTORY:
;   28-Feb-2003 Written by SB
;   18-Apr-2003 Revised by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro apf_fit2darc, apf, setup, obj_id, CLOBBER=clobber, $
                    nycoeff=nycoeff, nocoeff=nocoeff, DEBUG=debug, $
                    CHKRES=chkres, ARCFIL=arcfil, SIGREJ=sigrej
;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'apf_fit2darc, apf, setup, obj_id, [chip], nycoeff=, nocoeff=' 
      print, '     /CLOBBER, /CHKRES, /DEBUG, ARCFIL=, SIGREJ= [v1.1]'
      return
  endif

  if NOT keyword_set(nycoeff) then nycoeff = 4
  if NOT keyword_set(SIGREJ) then sigrej = 6.

  if keyword_set(nocoeff) then flg_ocoeff = 1L else flg_ocoeff = 0L

; Loop on chip
  ;; Nocoeff
  nocoeff = 3

  ;; Grab all obj indices
  if not keyword_set(ARCFIL) then begin
     indx = where(apf.flg_anly NE 0 AND $
                  apf.obj_id EQ obj_id AND apf.setup EQ setup AND $
                  strtrim(apf.type,2) EQ 'OBJ', nindx)
     if nindx EQ 0 then begin
        print, 'apf_fitarc: No Obj found!  Returning' 
        return
     endif
     arcfil_all= apf[indx].arc_fil
     asort = sort(arcfil_all)
     auniq = uniq(arcfil_all[asort])
  endif else begin
     arcfil_all = [arcfil]
     auniq = [0L]
     asort = [0L]
     pos = strpos(arcfil, '.fits')
     frame = long(strmid(arcfil,pos-5,5))
     indx = where(apf.frame EQ frame) 
  endelse

  
  ;; ORD_STR
  ordr_str = apf_getfil('ordr_str', setup) 

  ;; LOOP
  for mm=0L,n_elements(auniq)-1 do begin
     idx = indx[asort[auniq[mm]]]
     
     arcfil = arcfil_all[asort[auniq[mm]]]
     
     ;; Outfil
     pos = strpos(arcfil, '.fits')
     frame = long(strmid(arcfil,pos-5,5))
     
     ;;  Check for outfil
     outfil = apf_getfil('arc_2Dfit', setup, $
                         /name, CHKFIL=chkf,  FRAME=frame)
     if chkf NE 0 and not keyword_set( CLOBBER ) then begin
        print, 'apf_fit2darc: Arc fit file exists. ' + $
               'Continuing..'
        return
     endif

     ;;  QAFIL
     qafil = apf_getfil('qa_arc2dfit', setup, $
                        /name, CHKFIL=chkf,  FRAME=frame)
     ;; Open arc 1D fit
     arc_info = apf_getfil('arc_fit', setup, $
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
  endfor
  
   print, 'apf_fit2darc: All done'
   return

end

