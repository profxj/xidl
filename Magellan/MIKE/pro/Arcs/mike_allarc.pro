;+ 
; NAME:
; mike_allarc
;     Version 1.1
;
; PURPOSE:
;   Runs through all of the Arc processing steps in one go.
;   Can be used to process a set of arcs for a given object or all
;   of the arcs from a given night.
;
;   This file also contains mike_allarc_sngl which can be used to
;   fully process a single arc image.  In fact, if the keyword /ALL is
;   used, then mike_allarc simply loops through all of the arcs calling
;   mike_allarc_sngl for each one.
;
; CALLING SEQUENCE:
;   
;   mike_allarc, mike, setup, [side], /CHK, /CLOBBER, /INDX
;
; INPUTS:
;   mike     -  MIKE structure
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
;  It is best to consult the various routines to see the full set of
;  KEYWORD choices
;  /CLOBBER  - Overwite any previous files
;
; OPTIONAL OUTPUTS:
;  FITS=     - File to write mike structure to (default: 'strct.fits')
;
; ADDITIONAL KEYWORDS TO MIKE_ALLARC_SNGL:
;   /NOIMG    - Do not create arc image
;   /NOWAV    - Do not fit the arc lines (ie. only do tracing)
;   /NOONED    - Do not redo the 1D fitting (x_fitarc)
;
; COMMENTS:
;
; EXAMPLES:
;   mike_allarc, mike, 1L, 1L, /CHK, /CLOBBER
;   mike_allarc, mike, [104L, 107L], /INDX
;   rslt = mike_allarc_sngl('Raw/mb0539.fits', 1, 1, /PINTER)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  mike_fitarc
;  mike_fit2darc
;  mike_tracearc
;  mike_fittrcarc
;  mike_mkaimg
;
; REVISION HISTORY:
;   15-Aug-2003 Written by JXP
;
;  Usage:
;
;-
;

function mike_allarc_sngl, fil, setup, side, CHK=chk, CLOBBER=clobber, $
                           XYOFF=xyoff, TEMPL_FIL=templ_fil, $
                           IXY=ixy, NOONED=nooned, PINTER=pinter, NOIMG=noimg, $
                           NOWAV=nowav, BAD_ANLY=bad_anly, SHFTPRM=shftprm,$
                           _EXTRA=extra

;

  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'rslt = mike_allarc_sngl(fil, setup, side, /CHK, /CLOBBER' + $
        '/PINTER, XYOFF=, TEMPL_FIL=, IXY=, /NOIMG, /NOWAV, /NOONED, SHFTPRM=) [v1.1]'
      return, -1
  endif 

  if not keyword_set( IXY ) then ixy = [0., 0.]
;  if keyword_set(TEMPL_FIL) then stop

  ;; Check for arc_m
  ordr_str = mike_getfil('ordr_str', setup, SIDE=side)
  if ordr_str[0].arc_m[0] EQ 0. then begin
      ;; Run setarcm
      print, 'mike_allarc:  Running mike_setarcm..'
      mike_setarcm, fil, setup, side, ARC_FIL=arc_fil, XYOFF=xyoff, $
        FITPRM=shftprm, _EXTRA=extra
  endif

  ;; Process
  rslt = mike_procarc_sngl( fil, setup, side, CLOBBER=clobber, $
                          /ALIGN, CHK=keyword_set(chk), FITPRM=shftprm, $
                          _EXTRA=extra)

  arc_fil = rslt

  ;; Fit 1D
  if size(rslt,/tname) EQ 'STRING' and not keyword_set( NOWAV ) $
    AND not keyword_set(NOONED) then $
    rslt = mike_fitarc_work( arc_fil, setup, side, CLOBBER=clobber, $
                             GUESSARC=templ_fil,$
                             PINTER=pinter, _EXTRA=extra)
  ;; Fit 2D
  if size(rslt,/tname) EQ 'STRING' and not keyword_set( NOWAV ) then $
    rslt = mike_fit2darc_work( arc_fil, setup, side, CLOBBER=clobber, $
                               _EXTRA=extra)
  ;; Trace
  if size(rslt,/tname) EQ 'STRING' then $
    rslt = mike_tracearc_work( arc_fil, setup, side, CLOBBER=clobber, $
                             SHFTPRM=shftprm)
  ;; Fit Trace
  if size(rslt,/tname) EQ 'STRING' then $
    rslt = mike_fittrcarc_work( arc_fil, setup, side, CLOBBER=clobber )
  ;; Produce final image
;  if keyword_set( ALIGN ) then xoff = xyoff[0] else xoff = xyoff[0]+ixy[0]
  if size(rslt,/tname) EQ 'STRING' and not keyword_set( NOIMG ) then $
    rslt = mike_mkaimg_work( arc_fil, setup, side, CLOBBER=clobber, $
                             SHFTPRM=shftprm, CHK=chk, BAD_ANLY=bad_anly)

  print, 'mike_allarc_sngl: All done!'

  return, rslt
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_allarc, mike, setup, side, CHK=chk, CLOBBER=clobber, EXP=exp, $
                 PINTER=pinter, INDX=indx, FITS=fits, _EXTRA=extra


  if  N_params() LT 2 then begin 
      print,'Syntax - ' + $
        'mike_allarc, mike, setup, [side], /CHK, /PINTER, /CLOBBER,'
      print, '     /INTER, /INDX, FITS= [v1.1]'
      return
  endif 

  if keyword_set( OBJ_ID) then stop
  if not keyword_set( SIDE ) then side = [1L,2L]

  if (n_elements(indx) eq 0L) then begin
      ;; Loop on side
      sv_bad = 0
      for i=0,n_elements(side)-1 do begin
          ii = side[i]
          ;; NO objid
          arcs = where(mike.type EQ 'ARC' AND mike.flg_anly NE 0 $
                       AND mike.setup EQ setup $
                       AND mike.side EQ ii, narc) 
          ;; Exposure #
          if n_elements(EXP) NE 0 then begin
              arcs = arcs[exp]
              narc = n_elements(arcs)
          endif

          for qq=0L,narc-1 do begin
              bad_anly=0 
              rslt = mike_allarc_sngl(strtrim(mike[arcs[qq]].rootpth,2)+ $
                                      strtrim(mike[arcs[qq]].img_root,2),$
                                      mike[arcs[qq]].setup, $
                                      mike[arcs[qq]].side, $
                                      CLOBBER=clobber, BAD_ANLY=bad_anly, $
                                      SHFTPRM=shftprm, _EXTRA=extra, $
                                      PINTER=pinter, CHK=keyword_set(chk))
              
              if keyword_set(clobber) then mike[arcs[qq]].arc_xyoff = shftprm
              check_obj_xyoff = where(mike.arc_img EQ rslt AND $
                                      mike.obj_id GT 0, ncheck)
              if ncheck GT 0 then begin
                mike[check_obj_xyoff].arc_xyoff = $
                            mike[arcs[qq]].arc_xyoff # replicate(1,ncheck)
                print, 'Updating arc_xyoff for object numbers: ', $
                         string(mike[check_obj_xyoff].obj_id, format='(i4)')
              endif

              print, 'Bad ANLY = ', arcs[qq], bad_anly
              if bad_anly GT 0 then begin
                  mike[arcs[qq]].flg_anly = 0
                  sv_bad = 1
              endif
              ;; Write-out (just in case!)
              mike_wrstrct, mike
          endfor
      endfor
  endif else begin  ;; Single files
      sv_bad = 0
      ;; Loop on side
      nindx = n_elements(setup)
      for ii=0L,nindx-1 do begin
         ;print, ii, indx[ii]
          qq = setup[ii]
          if mike[qq].type NE 'ARC' then begin
              print, 'mike_allarc:  Index ', qq, ' not an arc!!  Fix and rerun'
              stop
          endif
          bad_anly=0 
          rslt = mike_allarc_sngl(strtrim(mike[qq].rootpth,2)+$
                                  strtrim(mike[qq].img_root,2),$
                                  mike[qq].setup, $
                                  mike[qq].side, $
                                  CLOBBER=clobber, BAD_ANLY=bad_anly, $
                                  SHFTPRM=shftprm, _EXTRA=extra, $
                                  PINTER=pinter, CHK=keyword_set(chk))
          if keyword_set(clobber) then mike[qq].arc_xyoff = shftprm

          check_obj_xyoff = where(mike.arc_img EQ rslt AND $
                                      mike.obj_id GT 0, ncheck)
          if ncheck GT 0 then begin
                mike[check_obj_xyoff].arc_xyoff = mike[qq].arc_xyoff # replicate(1,ncheck)
                print, 'Updating arc_xyoff for object numbers: ', $
                         string(mike[check_obj_xyoff].obj_id, format='(i4)')
          endif
          
          if bad_anly GT 0 then begin
              print, 'Bad ANLY = ', qq, bad_anly
              mike[qq].flg_anly = 0
              sv_bad = 1
          endif
          ;; Write-out (just in case!)
          mike_wrstrct, mike
      endfor
      ;; Rerun mike_setup
      if sv_bad NE 0 then mike_setup, mike
  endelse
                                                             
  ;; Save
  mike_wrstrct, mike, FITS=fits
  print, 'mike_allarc: All done!!!!!!!!!!!'

  return
end
