;+ 
; NAME:
; mike_fit2darc   
;     Version 1.1
;
; PURPOSE:
;  To fit the arc lines identified in mike_fitarc as a fucntion of
;  their y-centroid and order number.  The main routine is in
;  mike_fit2darc_work.  The fit is a simple least-squares with one
;  round of rejection.
;
; CALLING SEQUENCE:
;   
;  mike_fit2darc, mike
;
; INPUTS:
;   mike     -  MIKE structure
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
;   mike_fit2darc, mike, 1.0
;
; PROCEDURES/FUNCTIONS CALLED:
;  mike_fit2darc_work
;
; REVISION HISTORY:
;   28-Feb-2003 Written by SB
;   18-Apr-2003 Revised by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mike_fit2darc_work, arcfil, setup, side, nycoeff=nycoeff, $
                             nocoeff=nocoeff, CLOBBER=clobber, $
                             DEBUG=debug, CHKRES=chkres, out_str=out_str, $
                             OWGT=owgt, SIG2DREJ=sigrej
                   
;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'rslt = mike_fit2darc_work( arcfil, setup, side, nycoeff= ) [v2.0]'
      return, -1
  endif

   if NOT keyword_set(nycoeff) then nycoeff = 4
   if NOT keyword_set(sigrej) then sigrej = 10.
   if keyword_set(nocoeff) then flg_ocoeff = 1L else flg_ocoeff = 0L

   print, 'mike_fit2darc: Rejecting with ', sigrej

   ;; Nocoeff
   if flg_ocoeff NE 1 then begin
       case side of 
           1: nocoeff=5         ; Blue
           2: nocoeff=5         ; Red
           else: return, -1
       endcase
   endif

;  Check for outfil
   outfil = mike_getfil('arc_2Dfit', SUBFIL=arcfil, /name, CHKFIL=chkf)
   if chkf NE 0 and not keyword_set( CLOBBER ) then begin
       print, 'mike_fit2darc: Arc fit file exists. ' + $
         'Continuing..'
       return, outfil
   endif

; Order structure
   ordr_str = mike_getfil('ordr_str', setup, SIDE=side)
;   nordr = n_elements(ordr_str)
;   ordr_str[nordr-2:*].flg_anly = -1

; Open arc 1D fit
   arc_info = mike_getfil('arc_fit', subfil=arcfil, /name)
;   restore, arc_info

   qafil = mike_getfil('qa_arc2dfit', setup, SUBFIL=arcfil) 

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

pro mike_fit2darc, mike, setup, obj_id, side, CLOBBER=clobber, $
                   nycoeff=nycoeff, nocoeff=nocoeff, DEBUG=debug, CHKRES=chkres
;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'mike_fit2darc, mike, setup, obj_id, [side], nycoeff=, nocoeff=' + $
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
      if qq EQ 1 then print, 'mike_fitarc: Fitting BLUE arc' $
      else print, 'mike_fitarc: Fitting RED arc'

      ;; Nocoeff
      if flg_ocoeff NE 1 then begin
          case qq of 
              1: nocoeff=5 ; Blue
              2: nocoeff=6 ; Red
              else: stop
          endcase
      endif

      ;; Grab all obj indices
      indx = where(mike.flg_anly NE 0 AND mike.side EQ qq AND $
                   mike.obj_id EQ obj_id AND mike.setup EQ setup AND $
                   (strtrim(mike.type,2) EQ 'OBJ' $
                    OR strtrim(mike.type,2) EQ 'STD'), nindx)
      if nindx EQ 0 then begin
          print, 'mike_fitarc: No Obj found!  Returning' 
          continue
      endif


      arcfil_all= mike[indx].arc_fil
      asort = sort(arcfil_all)
      auniq = uniq(arcfil_all[asort])

  
      ;; LOOP
      for mm=0L,n_elements(auniq)-1 do begin
          idx = indx[asort[auniq[mm]]]

          rslt = mike_fit2darc_work( mike[idx].arc_fil, $
                                   setup, qq, nycoeff=nycoeff, CHKRES=chkres, $
                                   nocoeff=nocoeff, OUT_STR=out_str, $
                                   CLOBBER=clobber, DEBUG=debug)
          ;; Trouble?
          if size(rslt,/tname) NE 'STRING' then stop
      endfor
      
  endfor
;

   print, 'mike_fit2darc: All done'
   return

end

