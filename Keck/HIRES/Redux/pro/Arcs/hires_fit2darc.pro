;+ 
; NAME:
; hires_fit2darc   
;     Version 1.1
;
; PURPOSE:
;  To fit the arc lines identified in hires_fitarc as a fucntion of
;  their y-centroid and order number.  The main routine is in
;  hires_fit2darc_work.  The fit is a simple least-squares with one
;  round of rejection.
;
; CALLING SEQUENCE:
;   
;  hires_fit2darc, hires, setup, [obj_id, chip]
;
; INPUTS:
;   hires    -  HIRES structure
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
;   hires_fit2darc, hires, 1.0
;
; PROCEDURES/FUNCTIONS CALLED:
;  hires_fit2darc_work
;
; REVISION HISTORY:
;   28-Feb-2003 Written by SB
;   18-Apr-2003 Revised by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hires_fit2darc, hires, setup, obj_id, chip, CLOBBER=clobber, $
                    nycoeff=nycoeff, nocoeff=nocoeff, DEBUG=debug, $
                    CHKRES=chkres, ARCFIL=arcfil, SIGREJ=sigrej
;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hires_fit2darc, hires, setup, obj_id, [chip], nycoeff=, nocoeff=' 
      print, '     /CLOBBER, /CHKRES, /DEBUG, ARCFIL=, SIGREJ= [v1.1]'
      return
  endif

  if NOT keyword_set(nycoeff) then nycoeff = 4
  if NOT keyword_set(SIGREJ) then sigrej = 6.

  if keyword_set(ARCFIL) then begin
      ;; Grab Chip from name
      case strmid(arcfil,9,1) of
          'S': chip = -1L
          'B': chip = 1L
          'G': chip = 2L
          'R': chip = 3L
          else: stop
      endcase
  endif
  if not keyword_set( CHIP ) then chip = [1L,2L,3L]
  if keyword_set(nocoeff) then flg_ocoeff = 1L else flg_ocoeff = 0L

; Loop on chip
  for kk=0L,n_elements(chip)-1 do begin
      qq = chip[kk]
      ;; CHIP
      case qq of 
          -1: print, 'hires_fit2darc: 2D fit for Single chip'
          1: print, 'hires_fit2darc: Tracing BLUE CCD' 
          2: print, 'hires_fit2darc: Tracing GREEN CCD' 
          3: print, 'hires_fit2darc: Tracing RED CCD' 
          else: stop
      endcase

      ;; Nocoeff
      if flg_ocoeff NE 1 then begin
          ;; This is only good for the UV cross disperser for now
          nocoeff = 4
;          case qq of 
;              1: nocoeff=4 ; Blue
;              2: nocoeff=4 ; Green
;              3: nocoeff=4 ; Red
;              else: stop
;          endcase
      endif

      ;; Grab all obj indices
      if not keyword_set(ARCFIL) then begin
          indx = where(hires.flg_anly NE 0 AND hires.chip EQ qq AND $
                       hires.obj_id EQ obj_id AND hires.setup EQ setup AND $
                       strtrim(hires.type,2) EQ 'OBJ', nindx)
          if nindx EQ 0 then begin
              print, 'hires_fitarc: No Obj found!  Returning' 
              continue
          endif
          arcfil_all= hires[indx].arc_fil
          asort = sort(arcfil_all)
          auniq = uniq(arcfil_all[asort])
      endif else begin
          arcfil_all = [arcfil]
          auniq = [0L]
          asort = [0L]
          pos = strpos(arcfil, '.fits')
          frame = long(strmid(arcfil,pos-4,4))
          indx = where(hires.frame EQ frame AND hires.chip EQ chip)
      endelse

  
      ;; ORD_STR
      ordr_str = hires_getfil('ordr_str', setup, CHIP=chip)

      ;; LOOP
      for mm=0L,n_elements(auniq)-1 do begin
          idx = indx[asort[auniq[mm]]]

          arcfil = arcfil_all[asort[auniq[mm]]]

          ;; Outfil
          pos = strpos(arcfil, '.fits')
          frame = long(strmid(arcfil,pos-4,4))

          ;;  Check for outfil
          outfil = hires_getfil('arc_2Dfit', setup, CHIP=chip, $
                                 /name, CHKFIL=chkf,  FRAME=frame)
          if chkf NE 0 and not keyword_set( CLOBBER ) then begin
              print, 'hires_fit2darc: Arc fit file exists. ' + $
                'Continuing..'
              continue
          endif

          ;;  QAFIL
          qafil = hires_getfil('qa_arc2dfit', setup, CHIP=chip, $
                                 /name, CHKFIL=chkf,  FRAME=frame)
          ;; Open arc 1D fit
          arc_info = hires_getfil('arc_fit', setup, CHIP=chip, $
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
      
  endfor
;

   print, 'hires_fit2darc: All done'
   return

end

