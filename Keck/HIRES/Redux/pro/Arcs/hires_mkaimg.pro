;+ 
; NAME:
; hires_mkaimg
;     Version 1.1
;
; PURPOSE:
;   Given the 2D solution for the slope of the lines as a function of
;   position this code creates a wavelength image (i.e. assigns a
;   unique wavelength to each pixel in each order).  The xoffset is 
;   input in order to properly determine the edges of each order.  A
;   simple spline interpolation is used to determine the values.
;
;
; CALLING SEQUENCE:
;  hires_mkaimg, hires, setup, [obj_id, chip], /CHK, /CLOBBER, ARCFIL=
;
; INPUTS:
;   hires     -  MIKE structure
;   setup    -  Integer defining setup
;   [obj_id]   -  Object identifier
;   [chip]   -  Blue (1), Green (2), Red (3), or multiple [1,2L]
;            (Default: [1,2,3L])
;
; RETURNS:
;
; OUTPUTS:
;  2D wavelength image with name like 'Arcs/Arc_mb0439I.fits'
;
; OPTIONAL KEYWORDS:
;  /CLOBBER  - Overwrite previous image
;  /CHK      - Display the final image
;  ARCFIL=   - Name of the arc file to process (Optional to using
;               setup, chip, etc.)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hires_mkaimg, hires, setup, obj_id, chip
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   15-May-2003 Written by SB
;
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hires_mkaimg, hires, setup, obj_id, chip, CHK=chk, CLOBBER=clobber, $
                  ARCFIL=arcfil, BAD_ANLY=bad_anly, _EXTRA=extra

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hires_mkaimg, hires, setup, [obj_id, chip], /CHK, /CLOBBER, ARCFIL= [v1.1]'
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
          -1: print, 'hires_mkaimg: Making Single chip'
          1: print, 'hires_mkaimg: Making BLUE arc' 
          2: print, 'hires_mkaimg: Making GREEN arc' 
          3: print, 'hires_mkaimg: Making RED arc' 
      endcase

      if not keyword_set(ARCFIL) then begin
          ;; Grab all obj indices
          indx = where(hires.flg_anly NE 0 AND hires.chip EQ qq AND $
                       hires.obj_id EQ obj_id AND hires.setup EQ setup AND $
                       (strtrim(hires.type,2) EQ 'OBJ' OR $
                        strtrim(hires.type,2) EQ 'STD'), nindx)
          if nindx EQ 0 then begin
              print, 'hires_mkaimg: No Obj found!  Returning' 
              continue
          endif
          arcfil_all= hires[indx].arc_fil
          asort = sort(arcfil_all)
          auniq = uniq(arcfil_all[asort])
      endif else begin ;; Single frame
          arcfil_all = [arcfil]
          auniq = [0L]
          asort = [0L]
          pos = strpos(arcfil, '.fits')
          frame = long(strmid(arcfil,pos-4,4))
          indx = where(hires.frame EQ frame AND hires.chip EQ chip)
      endelse

      arc_roots = strmid(hires.arc_fil, max(strpos(hires.arc_fil, '.')-5), 5)

      ;; ORD_STR
      ordr_str = hires_getfil('ordr_str', setup, CHIP=chip, fil_nm=ordr_fil)

      ;; LOOP
      for mm=0L,n_elements(auniq)-1 do begin

          arc_fil = arcfil_all[asort[auniq[mm]]]
          idx = indx[asort[auniq[mm]]]

          ;; Outfil
          pos = strpos(arc_fil, '.fits')
          frame = long(strmid(arc_fil,pos-4,4))

          ;;  Check for outfil
          out_fil = hires_getfil('arc_mkaimg', setup, CHIP=chip, $
                                 /name, CHKFIL=chkf,  FRAME=frame)
          ;;  Check for outfil
          if chkf NE 0 and not keyword_set( CLOBBER ) then begin
              print, 'hires_fittrcarc: Arc fit file exists. ' + $
                'Continuing..'
              continue
          endif

          arc2d_fil = hires_getfil('arc_2Dfit', setup, CHIP=chip, $
                                 /name, CHKFIL=chkf,  FRAME=frame)
          fil_fittrc = hires_getfil('arc_fittrc', setup, CHIP=chip, $
                                 /name, CHKFIL=chkf,  FRAME=frame)
;          stop ;; Odds are that arc_xyoff is not set!
          rslt = x_mkaimg(arc_fil, ordr_str, arc2d_fil, fil_fittrc, $
                          out_fil, CHK=chk, CLOBBER=clobber, $
;                              SHFTPRM=hires[idx].arc_xyoff, $
                          BAD_ANLY=bad_anly, _EXTRA=extra) 

          if bad_anly GT 0 then stop
;             remove = where(strpos(hires.img_root, arc_roots[idx]) GT -1)
;             if remove[0] NE -1 then hires[remove].flg_anly = 0
;          endif

          if size(rslt,/tname) NE 'STRING' then stop
      endfor
  endfor

  print, 'hires_mkaimg: All done!'

  return
end
