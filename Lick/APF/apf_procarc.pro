;+ 
; NAME:
; apf_procarc   
;     Version 1.1
;
; PURPOSE:
;    Process (bias subtract, flatten) the arc files.  In addition, the
;    code first calls apf_arcxyoff which uses a 2D FFT to determine
;    the offset between the Template Arc and the current arc due to
;    thermal expansion in the instrument.
;   
;    This file also includes the routine apf_procarc_sngl which
;    allows the processing of a single Arc image given the filename of
;    the Raw image.
;
; CALLING SEQUENCE:
;   
;  apf_procarc, apf, setup, obj_id, ATEMPL=, /CLOBBER,
;  FLATFIL=
;  print, apf_procarc_sngl(raw_fil, apf, setup, [obj_id, 
;  exten], ATEMPL=, /CLOBBER)
;
; INPUTS:
;   apf     -  MIKE structure
;   setup    -  Integer defining setup
;   obj_id   -  Object identifier
;
; RETURNS:
;
; OUTPUTS:
;  One processed, combined Arc Image (e.g. Arcs/Arc_mb0539.fits)
;
; OPTIONAL KEYWORDS:
;   /CLOBBER - Overwrite exisiting arc image
;   ATEMPL -- Index of the Template Arc image (default: 0L)
;   FLATFIL -- Filename of the milky flat (pixel to pixel correction)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   apf_procarc, apf, 1, 1, 1, /CLOBBER
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  apf_getfil
;  apf_arcxyoff
;  apf_proc
;
; REVISION HISTORY:
;   13-Aug-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function apf_procarc_sngl, raw_fil, apf, setup, $
  TEMPL_FIL=templ_fil, FRAME=frame, $
  CLOBBER=clobber, FLATFIL=flatfil, CHK=chk, $
  XYOFF=xyoff, GAIN=gain, READNO=readno, ALIGN=align,$
  FITPRM=fitprm, _EXTRA=extra

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'rslt = apf_procarc_sngl( raw_fil, setup, ' + $
        'FLATFIL=, /CLOBBER, /ALIGN ) [v1.1]'
      return, -1
  endif 
  
;  Optional Keywords
  compile_opt strictarr

  if not keyword_set( ATEMPL ) then atempl = 0L

  ;; QA
  if setup GE 10 then c_s = strtrim(setup,2) else c_s = $
    '0'+strtrim(setup,2) 
  fqa = 'QA/Arcs'+c_s
  a = findfile(fqa, count=count)
  if count EQ 0 then file_mkdir, fqa
  
  ;; Check for prior image
  if not keyword_set( FRAME ) then begin
      pos = strpos(raw_fil, '.fits')
      frame = long(strmid(raw_fil,pos-4,4))
  endif
  outfil = apf_getfil('arc_fil', FRAME=frame,  /name, CHKFIL=chkfil)
  idx = where(apf.frame EQ frame AND apf.setup EQ setup)
  if CHKFIL NE 0 AND not keyword_set( CLOBBER ) then begin
      print, 'apf_procarc: Arc ', outfil, ' exists.  Returning..'
      return, outfil
  endif

  ;; Determine xy_off with the FFT
  ;; Process
  rslt = apf_proc_sngl(raw_fil, FLATFIL=flatfil, $
                         CLOBBER=clobber, $
                         REDOOV=clobber, /ARC, IOUTFIL=outfil, $
                         GAIN=apf[idx].gain, $
                         READNO=apf[idx].readno, SETUP=setup)

  ;; ALL DONE
  print, 'apf_procarc_sngl: All Done! '
  return, outfil

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro apf_procarc, apf, setup, obj_id, ATEMPL=atempl, $
                  CLOBBER=clobber, FLATFIL=flatfil, ALIGN=align

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'apf_procarc, apf, setup, obj_id, FLATFIL=, /CLOBBER ' + $
        'ATEMPL= [v1.1]'
      return
  endif 
  
  ;; QA
  if setup GE 10 then c_s = strtrim(setup,2) else c_s = '0'+strtrim(setup,2) 
  fqa = 'QA/Arcs'+c_s
  a = findfile(fqa, count=count)
  if count EQ 0 then file_mkdir, fqa

;  Optional Keywords
  if not keyword_set( ATEMPL ) then atempl = 0L
  
  ;; Grab all Arc files
  arcs = where(apf.setup EQ setup AND apf.flg_anly NE 0 AND $
               strtrim(apf.type,2) EQ 'ARC', narc)
  if narc EQ 0 then begin
     print, 'apf_procarc: No Arcs found! Returning' 
     return
  endif

  ;; Grab all obj files
  indx = where(apf.flg_anly NE 0 AND $
               apf.obj_id EQ obj_id AND apf.setup EQ setup AND $
               (strtrim(apf.type,2) EQ 'OBJ' OR $
                strtrim(apf.type,2) EQ 'STD'), nindx)
  if nindx EQ 0 then begin
     print, 'apf_procarc: No Obj found!  Returning' 
     stop
     return
  endif
  
  ;; LOOP
  for kk=0L,nindx-1 do begin
     
     print, 'apf_procarc: Index ', strtrim(kk,2)
     
     ;; Check for arc_fil name
     if strlen(strtrim(apf[indx[kk]].arc_fil,2)) EQ 0 then begin
        ;; Identify closest image in time
        diff = abs( (apf[indx[kk]].date + $
                     apf[indx[kk]].exp/2./3600./24.) - $
                    (apf[arcs].date + apf[arcs].exp/2./3600./24.) )
        mn = min(diff, gdarc)
        
        ;; Check for prior image
        outfil = apf_getfil('arc_fil', FRAME=apf[arcs[gdarc]].frame, $
                            CHIP=qq, /name)
     endif else begin
        outfil=apf[indx[kk]].arc_fil 
        pos = strpos(outfil, '.fits')
        frame = long(strmid(outfil,pos-4,4))
        gdarc = where(apf[arcs].frame EQ frame,ngd)
        if ngd NE 1 then stop
     endelse
     
     if x_chkfil(outfil+'*',/silent) NE 0 $
        AND not keyword_set( CLOBBER ) then begin
        print, 'apf_procarc: Arc ', outfil, ' exists.  ' + $
               'Updating using previous offsets and then returning'
        apf[indx[kk]].arc_fil = outfil
        continue
     endif
     
     ;; Determine xy_off
;          apf_arcxyoff, apf, qq, arcs[atempl], arcs[gdarc], xyoff
     
     ;; Process
     apf_proc, apf, arcs[gdarc], FLATFIL=flatfil, CLOBBER=clobber, $
               REDOOV=clobber, /ARC, IOUTFIL=outfil, setup=setup
     
     
     ;; Structure
     apf[indx[kk]].arc_fil = outfil
     
  endfor

  ;; ALL DONE
  print, 'apf_procarc: All Done! '
  return

end

