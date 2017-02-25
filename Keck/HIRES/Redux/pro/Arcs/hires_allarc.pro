;+ 
; NAME:
; hires_allarc
;     Version 1.1
;
; PURPOSE:
;   Runs through all of the Arc processing steps in one go.
;   Can be used to process a set of arcs for a given object or all
;   of the arcs from a given night.
;
;   This file also contains hires_allarc_sngl which can be used to
;   fully process a single arc image.  
;
; CALLING SEQUENCE:
;   
;   hires_allarc, hires, setup, [chip], /CHK, /CLOBBER, /INDX
;  print, hires_allarc_sngl(fil, hires, setup, chip, EXTEN=, FRAME=)
;
; INPUTS:
;   hires  -  HIRES structure
;   setup  -  Integer defining setup OR an array of index values to
;            process
;   [chip] -  Blue (1), Green (2), Red (3), or multiple [1,2L]
;            (Default: [1,2,3L])
;
; RETURNS:
;
; OUTPUTS:
;  A series of files related to Arc calibration
;
; OPTIONAL KEYWORDS:
;  /CLOBBER   - Overwite any previous files
;  /CHK       - Present GUIs to show results of fitting, etc.
;  /PINTER    - Perform fit for pre-identified lines
;  /INDX      - Treat setup as an array of index values to process
;
; OPTIONAL OUTPUTS:
;  FITS=     - File to write hires structure to (default: 'strct.fits')
;
; ADDITIONAL KEYWORDS TO HIRES_ALLARC_SNGL:
;   /NOIMG     - Do not create arc image
;   /NOWAV     - Do not fit the arc lines (ie. only do tracing)
;   /NOONED    - Do not redo the 1D fitting (x_fitarc)
;   /NOTRCFIT   - Do not fit the arc line traces [not recommended]
;   TEMPL_FIL= - Identify the template file for arc fitting
;
; COMMENTS:
;
; EXAMPLES:
;   hires_allarc, hires, 1L, 1L, /CHK, /CLOBBER
;   hires_allarc, hires, [104L, 107L], /INDX
;   rslt = hires_allarc_sngl('Raw/mb0539.fits', 1, 1, /PINTER)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  hires_fitarc
;  hires_fit2darc
;  hires_tracearc
;  hires_fittrcarc
;  hires_mkaimg
;
; REVISION HISTORY:
;   15-Aug-2003 Written by JXP
;
;  Usage:
;
;-
;

function hires_allarc_sngl, fil, hires, setup, chip, CHK=chk, CLOBBER=clobber, $
                            XYOFF=xyoff, TEMPL_FIL=templ_fil, NOONED=nooned,$
                            IXY=ixy, PINTER=pinter, NOIMG=noimg, $
                            NOWAV=nowav, BAD_ANLY=bad_anly, SHFTPRM=shftprm, $
                            EXTEN=exten, _EXTRA=extra, FRAME=frame, NOTRCFIT=notrcfit

;

  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'rslt = hires_allarc_sngl(fil, hires, setup, chip, /CHK, /CLOBBER' 
      print, '/PINTER, XYOFF=, TEMPL_FIL=, IXY=, /NOIMG, /NOWAV, SHFTPRM=, ' 
      print, 'EXTEN=, /NOONED, /NOTRCFIT) [v1.1]'
      return, -1
  endif 

  if not keyword_set( IXY ) then ixy = [0., 0.]
;  if keyword_set(TEMPL_FIL) then stop

  ;; Check for arc_m
  ordr_str = hires_getfil('ordr_str', setup, CHIP=chip, fil_nm=ordr_fil)
;  if ordr_str[0].arc_m[0] EQ 0. then begin
;      ;; Run setarcm
;      print, 'hires_allarc:  Running hires_setarcm..'
;      hires_setarcm, hires, fil, setup, chip, ARC_FIL=arc_fil, XYOFF=xyoff, $
;        FITPRM=shftprm, EXTEN=exten
;  endif


  ;; Process
;  rslt = hires_procarc_sngl( fil, setup, FRAME=frame, $
  rslt = hires_procarc_sngl( fil, hires, setup, FRAME=frame, $
                             chip, exten, CLOBBER=clobber, $
                             /ALIGN, CHK=keyword_set(chk), _EXTRA=extra)
  arc_fil = rslt

  ;; Fit 1D
  if size(rslt,/tname) EQ 'STRING' and not keyword_set( NOWAV ) $
    AND not keyword_set(NOONED) then $
    hires_fitarc, hires, setup, ARCFIL=arc_fil, CLOBBER=clobber, $
                             PINTER=pinter, GUESSARC=TEMPL_FIL, _EXTRA=extra
  ;; Fit 2D
  if size(rslt,/tname) EQ 'STRING' and not keyword_set( NOWAV ) then $
     hires_fit2darc, hires, setup, ARCFIL=arc_fil, CLOBBER=clobber, _EXTRA=extra 

  ;; Trace
  if size(rslt,/tname) EQ 'STRING' then $
    hires_tracearc, hires, setup, ARCFIL=arc_fil, CLOBBER=clobber
  ;; Fit Trace
  if size(rslt,/tname) EQ 'STRING' and not keyword_set(NOTRCFIT) then $
    hires_fittrcarc, hires, setup, ARCFIL=arc_fil, CLOBBER=clobber 

  ;; Produce final image
;  if keyword_set( ALIGN ) then xoff = xyoff[0] else xoff = xyoff[0]+ixy[0]
  if size(rslt,/tname) EQ 'STRING' and not keyword_set( NOIMG ) then $
    hires_mkaimg, hires, setup, ARCFIL=arc_fil, CLOBBER=clobber, $
    CHK=chk, BAD_ANLY=bad_anly, _EXTRA=extra

  print, 'hires_allarc_sngl: All done!'

  return, rslt
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hires_allarc, hires, allsetup, chip, CHK=chk, CLOBBER=clobber, $
                  PINTER=pinter, INDX=indx, FITS=fits, _EXTRA=extra, $
                  OBJ_ID=obj_id

;

  if  N_params() LT 2 then begin 
      print,'Syntax - ' + $
        'hires_allarc, hires, setup, [chip], /CHK, /PINTER, /CLOBBER,'
      print, '     /INTER, /INDX, FITS=, OBJ_ID=, /NOONED [v1.1]'
      return
  endif 
  if not keyword_set( CHIP ) then chip = [1L,2L,3L]

  ;; Loop on setups
  for kk=0L,n_elements(allsetup)-1 do begin
      setup = allsetup[kk]

      
      if not keyword_set(INDX) then begin
          ;; Loop on chip
          sv_bad = 0
          for i=0,n_elements(chip)-1 do begin
              ii = chip[i]
              ;; OBJID
              if keyword_set(OBJ_ID) then begin
                  ;; Arcfiles
                  obj = where((hires.type EQ 'OBJ' OR hires.type EQ 'STD') $
                              AND hires.flg_anly NE 0 $
                              AND hires.setup EQ setup and $
                              hires.obj_id EQ obj_id $
                              AND hires.chip EQ ii, narc) 
                  dumfil = hires_getarcfil(hires, obj, raw_fil=arcfil)
                  idx = obj
              endif else begin
                  ;; NO objid
                  arcs = where(hires.type EQ 'ARC' AND hires.flg_anly NE 0 $
                               AND hires.setup EQ setup $
                               AND hires.chip EQ ii, narc) 
                  arcfil= strtrim(hires[arcs].rootpth,2)+ $
                    strtrim(hires[arcs].img_root,2)
                  idx = arcs
              endelse
              narc = n_elements(arcfil)
              
              for qq=0L,narc-1 do begin
                  bad_anly=0 
                  rslt = hires_allarc_sngl(arcfil[qq], $
                                           hires, _EXTRA=extra, $
                                           hires[idx[qq]].setup, $
                                           hires[idx[qq]].chip, $
                                           frame=hires[idx[qq]].frame, $
                                           CLOBBER=clobber, $
                                           BAD_ANLY=bad_anly, $
                                           PINTER=pinter, CHK=keyword_set(chk), $
                                           EXTEN=hires[idx[qq]].exten)
;              hires[arcs[qq]].arc_xyoff = shftprm
                  
;                  print, 'Bad ANLY = ', arcs[qq], bad_anly
                  if bad_anly GT 0 and keyword_set(ARCS) then begin
                      hires[arcs[qq]].flg_anly = 0
                      sv_bad = 1
                  endif
                  ;; Write-out (just in case!)
                  hires_wrstrct, hires
              endfor
          endfor
      endif else begin  ;; Single files
          stop  ;; Not supported right now  JXP  
          ;; Loop on chip
          nindx = n_elements(setup)
          for ii=0L,nindx-1 do begin
              qq = indx[ii]
              if hires[qq].type NE 'ARC' then begin
                  print, 'hires_allarc:  Index ', qq, ' not an arc!!  Fix and rerun'
                  stop
              endif
              bad_anly=0 
              rslt = hires_allarc_sngl(strtrim(hires[qq].rootpth,2)+$
                                       strtrim(hires[qq].img_root,2),$
                                       hires[qq].setup, $
                                       hires[qq].chip, $
                                       FRAME=hires[qq].frame, $
                                       CLOBBER=clobber, BAD_ANLY=bad_anly, $
                                       SHFTPRM=shftprm, $
                                       PINTER=pinter, CHK=keyword_set(chk))
              hires[qq].arc_xyoff = shftprm
              
              if bad_anly GT 0 then begin
                  print, 'Bad ANLY = ', qq, bad_anly
                  hires[qq].flg_anly = 0
                  sv_bad = 1
              endif
              ;; Write-out (just in case!)
              hires_wrstrct, hires
          endfor
          ;; Rerun hires_setup
          if sv_bad NE 0 then hires_setup, hires
      endelse
  endfor
      
  ;; Save
  hires_wrstrct, hires, FITS=fits
  print, 'hires_allarc: All done!!!!!!!!!!!'

  return
end
