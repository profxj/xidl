;+ 
; NAME:
; uves_fit2darc   
;     Version 1.1
;
; PURPOSE:
;  To fit the arc lines identified in uves_fitarc as a fucntion of
;  their y-centroid and order number.  The main routine is in
;  uves_fit2darc_work.  The fit is a simple least-squares with one
;  round of rejection.
;
; CALLING SEQUENCE:
;   
;  uves_fit2darc, uves
;
; INPUTS:
;   uves     -  MIKE structure
;   setup    -  Integer defining setup
;   obj_id   -  Object identifier
;   [side]   -  Blue (1), Red (2), or both [1,2L]    (Default: [1,2L])
;
; RETURNS:
;
; OUTPUTS:
; A fits file containing the 2D solution.  Named something like
; 'Arcs/Fits/Arc_mb0539_fit2D.fits' 
;
; OPTIONAL KEYWORDS:
;   NOCOEFF   - Number of coefficients to use in the x-direction
;               (default: 6 for blue, 7 for red)
;   NYCOEFF   - Number of coefficients to use in the y-direction
;               (default: 4)
;   /CLOBBER  - Overwrite any previous solution
;   /DEBUG    - debug
;   /CHKRES   - Plot the residuals
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   uves_fit2darc, uves, 1.0
;
; PROCEDURES/FUNCTIONS CALLED:
;  uves_fit2darc_work
;
; REVISION HISTORY:
;   28-Feb-2003 Written by SB
;   18-Apr-2003 Revised by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function uves_fit2darc_work, arcfil, setup, side, nycoeff=nycoeff, $
                             nocoeff=nocoeff, CLOBBER=clobber, $
                             DEBUG=debug, CHKRES=chkres, out_str=out_str, $
                             OWGT=owgt, SIG2DREJ=sigrej, _EXTRA=extra
                   
;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'rslt = uves_fit2darc_work( arcfil, setup, side, nycoeff= ) [v2.0]'
      return, -1
  endif

   if NOT keyword_set(nycoeff) then nycoeff = 4
   if NOT keyword_set(sigrej) then sigrej = 5.
   if keyword_set(nocoeff) then flg_ocoeff = 1L else flg_ocoeff = 0L

   ;; Nocoeff
   if flg_ocoeff NE 1 then begin
       case side of 
           1: nocoeff=5         ; Blue
           2: nocoeff=5         ; Red
           else: return, -1
       endcase
   endif

   ;; wcen
   wcen = uves_getwcen(arcfil,FRAME=frame, /AFIL)

;  Check for outfil
   outfil = uves_getfil('arc_2Dfit', setup, FRAME=frame, $
                        WCEN=wcen, /name, CHKFIL=chkf)
   if chkf NE 0 and not keyword_set( CLOBBER ) then begin
       print, 'uves_fit2darc: Arc fit file exists. ' + $
         'Continuing..'
       return, outfil
   endif

; Order structure
   ordr_str = uves_getfil('ordr_str', setup, WCEN=wcen)

; Open arc 1D fit
   arc_info = uves_getfil('arc_fit', setup, WCEN=wcen, FRAME=frame, /name)
   qafil = uves_getfil('qa_arc2dfit', setup, WCEN=wcen, FRAME=frame, /name)

   rslt = x_fit2darc( arcfil, ordr_str, arc_info, $
                      nycoeff=nycoeff, CHKRES=chkres, SIGREJ=sigrej, $
                      nocoeff=nocoeff, OUT_STR=out_str, $
                      CLOBBER=clobber, DEBUG=debug, QAFIL=qafil);, /yval)

   ;; Output
   mwrfits, out_str, outfil, /create

   return, outfil

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro uves_fit2darc, uves, setup, obj_id, side, CLOBBER=clobber, ARCFIL=arcfil, $
                   nycoeff=nycoeff, nocoeff=nocoeff, DEBUG=debug, CHKRES=chkres
;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'uves_fit2darc, uves, setup, obj_id, [side], nycoeff=, nocoeff=' + $
        '/CLOBBER, /CHKRES, /DEBUG [v2.0]'
      return
  endif

   if NOT keyword_set(nycoeff) then nycoeff = 4
   if not keyword_set( SIDE ) then side = [1L,2L]
   if keyword_set(nocoeff) then flg_ocoeff = 1L else flg_ocoeff = 0L

; Label
  if setup GE 10 then c_s = strtrim(setup,2) else c_s = '0'+strtrim(setup,2) 

; Loop on side
  for kk=0L,n_elements(side)-1 do begin
      qq = side[kk]
      ;; SIDE
      if qq EQ 1 then print, 'uves_fit2darc: Fitting BLUE arc' $
      else print, 'uves_fit2darc: Fitting RED arc'

      ;; Nocoeff
      if flg_ocoeff NE 1 then begin
          case qq of 
              1: nocoeff=5 ; Blue
              2: nocoeff=5 ; Red
              else: stop
          endcase
      endif

      ;; Grab all obj indices
      indx = where(uves.flg_anly NE 0 AND uves.side EQ qq AND $
                   uves.obj_id EQ obj_id AND uves.setup EQ setup AND $
                   (strtrim(uves.type,2) EQ 'OBJ' $
                    OR strtrim(uves.type,2) EQ 'STD'), nindx)
      if nindx EQ 0 then begin
          print, 'uves_fitarc: No Obj found!  Returning' 
          continue
      endif


      arcfil_all= uves[indx].arc_fil
      asort = sort(arcfil_all)
      auniq = uniq(arcfil_all[asort])

  
      ;; LOOP
      for mm=0L,n_elements(auniq)-1 do begin
          idx = indx[asort[auniq[mm]]]

          rslt = uves_fit2darc_work( uves[idx].arc_fil, $
                                   setup, qq, nycoeff=nycoeff, CHKRES=chkres, $
                                   nocoeff=nocoeff, OUT_STR=out_str, $
                                   CLOBBER=clobber, DEBUG=debug)
          ;; Trouble?
          if size(rslt,/tname) NE 'STRING' then stop
      endfor
      
  endfor
;

   print, 'uves_fit2darc: All done'
   return

end

