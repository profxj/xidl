;+ 
; NAME:
; hires_fittrcarc
;     Version 1.1
;
; PURPOSE:
;   To fit the slope of the arc lines as a function of order number
;   and y position on the CCD.  This information is then used to
;   construct a 2D wavelength image.  The fitting routine is the usual
;   least-squares algorithm with two rounds of rejection.
;
; CALLING SEQUENCE:
; hires_fittrcarc, hires, setup, obj_id, chip, /CHK, /CLOBBER, $
;                   ARCFIL=arcfil
;
; INPUTS:
;   hires    -  HIRES structure
;   setup    -  Integer defining setup
;   obj_id   -  Object identifier
;   [chip]   -  Blue (1), Green (2), Red (3), or multiple [1,2L]
;            (Default: [1,2,3L])
;
; RETURNS:
;
; OUTPUTS:
;  Fits file with the coefficients of the 2D fit.  Filename like
;  'Arcs/TRC/Arc_B0539_F.fits' 
;
; OPTIONAL KEYWORDS:
;  /CHK  -- Plots residuals
;  /CLOBBER -- Overwrite previous solution
;  ARCFIL=   - Name of the arc file to process (Optional to using
;               setup, chip, etc.)
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
;   24-Feb-2005 Written by JXP
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hires_fittrcarc, hires, setup, obj_id, chip, CHK=chk, CLOBBER=clobber, $
                     ARCFIL=arcfil
;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hires_fittrcarc, hires, setup, obj_id, [chip], /CHK, /CLOBB, ARCFIL= [v1.1]'
      return
  endif 

;  Optional Keywords
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

; Loop on chip
  for kk=0L,n_elements(chip)-1 do begin
      qq = chip[kk]

      ;; CHIP
      case qq of
          -1: print, 'hires_fittrcarc: Fitting Single arc'
          1: print, 'hires_fittrcarc: Fitting BLUE arc' 
          2: print, 'hires_fittrcarc: Fitting GREEN arc' 
          3: print, 'hires_fittrcarc: Fitting RED arc' 
      endcase

      if not keyword_set(ARCFIL) then begin
          ;; Grab all obj indices
          indx = where(hires.flg_anly NE 0 AND hires.chip EQ qq AND $
                       hires.obj_id EQ obj_id AND hires.setup EQ setup AND $
                       (strtrim(hires.type,2) EQ 'OBJ' OR $ 
                        strtrim(hires.type,2) EQ 'STD'), nindx)
          if nindx EQ 0 then begin
              print, 'hires_fittrcarc: No Obj found!  Returning' 
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
      ordr_str = hires_getfil('ordr_str', setup, CHIP=chip, fil_nm=ordr_fil)

      ;; LOOP on exposures
      for mm=0L,n_elements(auniq)-1 do begin
          arc_fil = arcfil_all[asort[auniq[mm]]]
          idx = indx[asort[auniq[mm]]]

          ;; Outfil
          pos = strpos(arc_fil, '.fits')
          frame = long(strmid(arc_fil,pos-4,4))

          ;;  Check for outfil
          out_fil = hires_getfil('arc_fittrc', setup, CHIP=chip, $
                                 /name, CHKFIL=chkf,  FRAME=frame)
          ;;  Check for outfil
          if chkf NE 0 and not keyword_set( CLOBBER ) then begin
              print, 'hires_fittrcarc: Arc fit file exists. ' + $
                'Continuing..'
              continue
          endif

          ;; TRC_FIL
          trc_fil = hires_getfil('arc_trc', setup, CHIP=chip, $
                                 /name, CHKFIL=chkf,  FRAME=frame)
          ;; QA
          qafil = hires_getfil('qa_fittrcarc', setup, CHIP=chip, $
                                 /name, CHKFIL=chkf,  FRAME=frame)

          ;; Main driver
          rslt = x_fittrcarc(arc_fil, trc_fil, ordr_str, out_fil, qafil, $
                             CHK=chk, CLOBBER=clobber, $
                             ORDR_FIL=ordr_fil) 
          if size(rslt,/tname) NE 'STRING' then stop
      endfor
  endfor

  print, 'hires_fittrcarc: All done!'

  return
end
