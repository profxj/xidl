;+ 
; NAME:
; apf_fittrcarc
;     Version 1.1
;
; PURPOSE:
;   To fit the slope of the arc lines as a function of order number
;   and y position on the CCD.  This information is then used to
;   construct a 2D wavelength image.  The fitting routine is the usual
;   least-squares algorithm with two rounds of rejection.
;
; CALLING SEQUENCE:
; apf_fittrcarc, apf, setup, obj_id, /CHK, /CLOBBER, $
;                   ARCFIL=arcfil
;
; INPUTS:
;   apf    -  HIRES structure
;   setup    -  Integer defining setup
;   obj_id   -  Object identifier
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
;               setup, etc.)
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

pro apf_fittrcarc, apf, setup, obj_id, CHK=chk, CLOBBER=clobber, $
                     ARCFIL=arcfil
;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'apf_fittrcarc, apf, setup, obj_id, /CHK, /CLOBB, ARCFIL= [v1.1]'
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
        print, 'apf_fittrcarc: No Obj found!  Returning' 
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
  ordr_str = apf_getfil('ordr_str', setup, fil_nm=ordr_fil)
  
  ;; LOOP on exposures
  for mm=0L,n_elements(auniq)-1 do begin
     arc_fil = arcfil_all[asort[auniq[mm]]]
     idx = indx[asort[auniq[mm]]]
     
     ;; Outfil
     pos = strpos(arc_fil, '.fits')
     frame = long(strmid(arc_fil,pos-5,5))
     
     ;;  Check for outfil
     out_fil = apf_getfil('arc_fittrc', setup, $
                          /name, CHKFIL=chkf,  FRAME=frame)
     ;;  Check for outfil
     if chkf NE 0 and not keyword_set( CLOBBER ) then begin
        print, 'apf_fittrcarc: Arc fit file exists. ' + $
               'Continuing..'
        return
     endif
     
     ;; TRC_FIL
     trc_fil = apf_getfil('arc_trc', setup, $
                          /name, CHKFIL=chkf,  FRAME=frame)
     ;; QA
     qafil = apf_getfil('qa_fittrcarc', setup, $
                        /name, CHKFIL=chkf,  FRAME=frame)
     
     ;; Main driver
     rslt = x_fittrcarc(arc_fil, trc_fil, ordr_str, out_fil, qafil, $
                        CHK=chk, CLOBBER=clobber, $
                        ORDR_FIL=ordr_fil) 
     if size(rslt,/tname) NE 'STRING' then stop
  endfor

  print, 'apf_fittrcarc: All done!'

  return
end
