;+ 
; NAME:
; apf_mkaimg
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
;  apf_mkaimg, apf, setup, [obj_id], /CHK, /CLOBBER, ARCFIL=
;
; INPUTS:
;   apf     -  MIKE structure
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
;   apf_mkaimg, apf, setup, obj_id, chip
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   15-May-2003 Written by SB
;
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro apf_mkaimg, apf, setup, obj_id, CHK=chk, CLOBBER=clobber, $
                  ARCFIL=arcfil, BAD_ANLY=bad_anly, _EXTRA=extra

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'apf_mkaimg, apf, setup, [obj_id, /CHK, /CLOBBER, ARCFIL= [v1.1]'
      return
  endif 

;  Optional Keywords

  if not keyword_set(ARCFIL) then begin
     ;; Grab all obj indices
     indx = where(apf.flg_anly NE 0 AND $
                  apf.obj_id EQ obj_id AND apf.setup EQ setup AND $
                  (strtrim(apf.type,2) EQ 'OBJ' OR $
                   strtrim(apf.type,2) EQ 'STD'), nindx)
     if nindx EQ 0 then begin
        print, 'apf_mkaimg: No Obj found!  Returning' 
        return
     endif
     arcfil_all= apf[indx].arc_fil
     asort = sort(arcfil_all)
     auniq = uniq(arcfil_all[asort])
  endif else begin ;; Single frame
     arcfil_all = [arcfil]
     auniq = [0L]
     asort = [0L]
     pos = strpos(arcfil, '.fits')
     frame = long(strmid(arcfil,pos-5,5))
     indx = where(apf.frame EQ frame)
  endelse
  
  arc_roots = strmid(apf.arc_fil, max(strpos(apf.arc_fil, '.')-5), 5)

  ;; ORD_STR
  ordr_str = apf_getfil('ordr_str', setup, fil_nm=ordr_fil)
  
  ;; LOOP
  for mm=0L,n_elements(auniq)-1 do begin
     
     arc_fil = arcfil_all[asort[auniq[mm]]]
     idx = indx[asort[auniq[mm]]]
     
     ;; Outfil
     pos = strpos(arc_fil, '.fits')
     frame = long(strmid(arc_fil,pos-5,5))
     
     ;;  Check for outfil
     out_fil = apf_getfil('arc_mkaimg', setup, $
                          /name, CHKFIL=chkf,  FRAME=frame)
     ;;  Check for outfil
     if chkf NE 0 and not keyword_set( CLOBBER ) then begin
        print, 'apf_fittrcarc: Arc fit file exists. ' + $
               'Continuing..'
        continue
     endif
     
     arc2d_fil = apf_getfil('arc_2Dfit', setup, $
                            /name, CHKFIL=chkf,  FRAME=frame)
     fil_fittrc = apf_getfil('arc_fittrc', setup, $
                             /name, CHKFIL=chkf,  FRAME=frame)
;          stop ;; Odds are that arc_xyoff is not set!
     rslt = x_mkaimg(arc_fil, ordr_str, arc2d_fil, fil_fittrc, $
                     out_fil, CHK=chk, CLOBBER=clobber, $
;                              SHFTPRM=apf[idx].arc_xyoff, $
                     BAD_ANLY=bad_anly, NOEXTRA=1, _EXTRA=extra) 

     if bad_anly GT 0 then stop
;             remove = where(strpos(apf.img_root, arc_roots[idx]) GT -1)
;             if remove[0] NE -1 then apf[remove].flg_anly = 0
;          endif

     if size(rslt,/tname) NE 'STRING' then stop
  endfor
  
  print, 'apf_mkaimg: All done!'

  return
end
