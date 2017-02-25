;+ 
; NAME:
; uves_allarc
;     Version 1.1
;
; PURPOSE:
;   Runs through all of the Arc processing steps in one go.
;   Can be used to process a set of arcs for a given object or all
;   of the arcs from a given night.
;
;   This file also contains uves_allarc_sngl which can be used to
;   fully process a single arc image.  In fact, if the keyword /ALL is
;   used, then uves_allarc simply loops through all of the arcs calling
;   uves_allarc_sngl for each one.
;
; CALLING SEQUENCE:
;   
;   uves_allarc, uves, setup, [side], /CHK, /CLOBBER, /INDX
;
; INPUTS:
;   uves     -  MIKE structure
;   setup  -  Integer defining setup OR an array of index values to
;            process
;   [side]   -  Blue (1), Red (2), or both [1,2L]    (Default: [1,2L])
;
; RETURNS:
;
; OUTPUTS:
;  A series of files related to Arc calibration
;
; OPTIONAL KEYWORDS:
;  /CLOBBER  - Overwite any previous files
;  /CHK      - Present GUIs to show results of fitting, etc.
;  /PINTER   - Perform fit for pre-identified lines
;  /INDX     - Treat setup as an array of index values to process
;
; OPTIONAL OUTPUTS:
;  FITS=     - File to write uves structure to (default: 'strct.fits')
;
; ADDITIONAL KEYWORDS TO MIKE_ALLARC_SNGL:
;   /NOIMG    - Do not create arc image
;   /NOWAV    - Do not fit the arc lines (ie. only do tracing)
;
; COMMENTS:
;
; EXAMPLES:
;   uves_allarc, uves, 1L, 1L, /CHK, /CLOBBER
;   uves_allarc, uves, [104L, 107L], /INDX
;   rslt = uves_allarc_sngl('Raw/mb0539.fits', 1, 1, /PINTER)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  uves_fitarc
;  uves_fit2darc
;  uves_tracearc
;  uves_fittrcarc
;  uves_mkaimg
;
; REVISION HISTORY:
;   15-Aug-2003 Written by JXP
;
;  Usage:
;
;-
;

function uves_allarc_sngl, fil, uves, setup, side, CHK=chk, CLOBBER=clobber, $
                            XYOFF=xyoff, TEMPL_FIL=templ_fil, NOONED=nooned,$
                            IXY=ixy, PINTER=pinter, NOIMG=noimg, $
                            NOWAV=nowav, BAD_ANLY=bad_anly, SHFTPRM=shftprm, $
                            EXTEN=exten, _EXTRA=extra, FRAME=frame

;

  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'rslt = uves_allarc_sngl(fil, uves, setup, side, /CHK, /CLOBBER' 
      print, '/PINTER, XYOFF=, TEMPL_FIL=, IXY=, /NOIMG, /NOWAV, SHFTPRM=, ' 
      print, 'EXTEN=, /NOONED) [v1.1]'
      return, -1
  endif 

  if not keyword_set( IXY ) then ixy = [0., 0.]
  if side EQ 2 then exten = 1 else exten = 0
  if keyword_set(TEMPL_FIL) then stop

  ;; Check for arc_m
  wcen = uves_getwcen(fil,uves,FRAME=frame)
  ordr_str = uves_getfil('ordr_str', setup, WCEN=wcen, fil_nm=ordr_fil)
;  if ordr_str[0].arc_m[0] EQ 0. then begin
;      ;; Run setarcm
;      print, 'uves_allarc:  Running uves_setarcm..'
;      uves_setarcm, uves, fil, setup, side, ARC_FIL=arc_fil, XYOFF=xyoff, $
;        FITPRM=shftprm, EXTEN=exten
;  endif


  ;; Process
;  rslt = uves_procarc_sngl( fil, setup, FRAME=frame, $
  rslt = uves_procarc_sngl( fil, setup, side, uves, FRAME=frame, CLOBBER=clobber, $
                            /ALIGN, EXTEN=exten, $
                            CHK=keyword_set(chk), _EXTRA=extra)
  arc_fil = rslt

  ;; Fit 1D
  if size(rslt,/tname) EQ 'STRING' and not keyword_set( NOWAV ) $
    AND not keyword_set(NOONED) then $
    uves_fitarc, uves, setup, side, ARCFIL=arc_fil, CLOBBER=clobber, $
                             PINTER=pinter, _EXTRA=extra
  ;; Fit 2D
  if size(rslt,/tname) EQ 'STRING' and not keyword_set( NOWAV ) then $
     print, uves_fit2darc_work(arc_fil, setup, side, CLOBBER=clobber )

  ;; Trace
  if size(rslt,/tname) EQ 'STRING' then $
    print, uves_tracearc_work(arc_fil, setup, side, CLOBBER=clobber )
  ;; Fit Trace
  if size(rslt,/tname) EQ 'STRING' then $
    print, uves_fittrcarc_work(arc_fil, setup, side, CLOBBER=clobber )

  ;; Produce final image
;  if keyword_set( ALIGN ) then xoff = xyoff[0] else xoff = xyoff[0]+ixy[0]
  if size(rslt,/tname) EQ 'STRING' and not keyword_set( NOIMG ) then $
    print, uves_mkaimg_work(arc_fil, setup, side, CLOBBER=clobber, $
                            CHK=chk, BAD_ANLY=bad_anly)

  print, 'uves_allarc_sngl: All done!'

  return, rslt
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro uves_allarc, uves, allsetup, side, CHK=chk, CLOBBER=clobber, $
                 PINTER=pinter, INDX=indx, FITS=fits, _EXTRA=extra, $
                  OBJ_ID=obj_id

;

  if  N_params() LT 2 then begin 
      print,'Syntax - ' + $
        'uves_allarc, uves, setup, [side], /CHK, /PINTER, /CLOBBER,'
      print, '     /INTER, /INDX, FITS=, OBJ_ID= [v1.1]'
      return
  endif 
  if not keyword_set( SIDE ) then side = [1L,2L]

  ;; Loop on setups
  for kk=0L,n_elements(allsetup)-1 do begin
      setup = allsetup[kk]

      
      if not keyword_set(INDX) then begin
          ;; Loop on side
          sv_bad = 0
          for i=0,n_elements(side)-1 do begin
              ii = side[i]
              ;; OBJID
              if keyword_set(OBJ_ID) then stop $
                  ;; Arcfiles
;                  obj = where((uves.type EQ 'OBJ' OR uves.type EQ 'STD') $
;                              AND uves.flg_anly NE 0 $
;                              AND uves.setup EQ setup and $
;                              uves.obj_id EQ obj_id $
;                              AND uves.side EQ ii, narc) 
;                  dumfil = uves_getarcfil(uves, obj, raw_fil=arcfil)
;                  idx = obj
              else begin
                  ;; NO objid
                  arcs = where(uves.type EQ 'ARC' AND uves.flg_anly NE 0 $
                               AND uves.setup EQ setup $
                               AND uves.side EQ ii, narc) 
                  arcfil= strtrim(uves[arcs].rootpth,2)+ $
                    strtrim(uves[arcs].img_root,2)
                  idx = arcs
              endelse
              narc = n_elements(arcfil)
              
              for qq=0L,narc-1 do begin
                  bad_anly=0 
                  rslt = uves_allarc_sngl(arcfil[qq], $
                                           uves, _EXTRA=extra, $
                                           uves[idx[qq]].setup, $
                                           uves[idx[qq]].side, $
                                           frame=uves[idx[qq]].frame, $
                                           CLOBBER=clobber, $
                                           BAD_ANLY=bad_anly, $
                                           PINTER=pinter, CHK=keyword_set(chk), $
                                           EXTEN=uves[idx[qq]].exten)
;              uves[arcs[qq]].arc_xyoff = shftprm
                  
;                  print, 'Bad ANLY = ', arcs[qq], bad_anly
                  if bad_anly GT 0 and keyword_set(ARCS) then begin
                      uves[arcs[qq]].flg_anly = 0
                      sv_bad = 1
                  endif
                  ;; Write-out (just in case!)
                  uves_wrstrct, uves
              endfor
          endfor
      endif else begin  ;; Single files
          stop  ;; Not supported right now  JXP  
          ;; Loop on side
          nindx = n_elements(setup)
          for ii=0L,nindx-1 do begin
              qq = indx[ii]
              if uves[qq].type NE 'ARC' then begin
                  print, 'uves_allarc:  Index ', qq, ' not an arc!!  Fix and rerun'
                  stop
              endif
              bad_anly=0 
              rslt = uves_allarc_sngl(strtrim(uves[qq].rootpth,2)+$
                                       strtrim(uves[qq].img_root,2),$
                                       uves[qq].setup, $
                                       uves[qq].side, $
                                       FRAME=uves[qq].frame, $
                                       CLOBBER=clobber, BAD_ANLY=bad_anly, $
                                       SHFTPRM=shftprm, $
                                       PINTER=pinter, CHK=keyword_set(chk))
              uves[qq].arc_xyoff = shftprm
              
              if bad_anly GT 0 then begin
                  print, 'Bad ANLY = ', qq, bad_anly
                  uves[qq].flg_anly = 0
                  sv_bad = 1
              endif
              ;; Write-out (just in case!)
              uves_wrstrct, uves
          endfor
          ;; Rerun uves_setup
          if sv_bad NE 0 then uves_setup, uves
      endelse
  endfor
      
  ;; Save
  uves_wrstrct, uves, FITS=fits
  print, 'uves_allarc: All done!!!!!!!!!!!'

  return
end
