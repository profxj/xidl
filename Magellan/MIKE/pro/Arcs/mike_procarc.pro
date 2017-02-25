;+ 
; NAME:
; mike_procarc   
;     Version 1.1
;
; PURPOSE:
;    Process (bias subtract, flatten) the arc files.  In addition, the
;    code first calls mike_arcxyoff which uses a 2D FFT to determine
;    the offset between the Template Arc and the current arc due to
;    thermal expansion in the instrument.
;   
;    This file also includes the routine mike_procarc_sngl which
;    allows the processing of a single Arc image given the filename of
;    the Raw image.
;
; CALLING SEQUENCE:
;   
;  mike_procarc, mike, setup, obj_id, [side], ATEMPL=, /CLOBBER,
;  fil = mike_procarc_sngl( raw_fil, setup, side )
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
;  One processed, combined Arc Image (e.g. Arcs/Arc_mb0539.fits)
;
; OPTIONAL KEYWORDS:
;   /CLOBBER - Overwrite exisiting arc image
;   /ALIGN  -- Deterimne offset between the Arc and the trace flat.
;              This deals with thermal gradients in the instrument.
;              It is HIGHLY recommended.
;   FLATFIL -- Filename of the milky flat (pixel to pixel correction)
;   READNO,GAIN  -- Input values of readnoise and gain
;
; OPTIONAL OUTPUTS:
;   XYOFF -- Offset values between the Arc and the Trace Flat
;   FITPRM -- Fit values to the x-offset as a function of order#
;
; COMMENTS:
;
; EXAMPLES:
;   mike_procarc, mike, 1, 1, 1, /CLOBBER
;   print, mike_procarc_sngl('Raw/mb0133.fits', /ALIGN, FITPRM=fprm)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  mike_getfil
;  mike_arcxyoff
;  mike_proc
;
; REVISION HISTORY:
;   13-Aug-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mike_procarc_sngl, raw_fil, setup, side, $
                            CLOBBER=clobber, FLATFIL=flatfil, CHK=chk, $
                            XYOFF=xyoff, GAIN=gain, READNO=readno, $
                            ALIGN=align, FITPRM=fitprm, _EXTRA=extra

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'rslt = mike_procarc_sngl( raw_fil, setup, side, FLATFIL=, /CLOBBER, /ALIGN ) [v1.1]'
      return, -1
  endif 
  
;  Optional Keywords
  compile_opt strictarr

  if not keyword_set( ATEMPL ) then atempl = 0L

  ;; QA
  if setup GE 10 then c_s = strtrim(setup,2) else c_s = '0'+strtrim(setup,2) 
  fqa = 'QA/Arcs'+c_s
  a = findfile(fqa, count=count)
  if count EQ 0 then file_mkdir, fqa
  
  ;; Check for prior image
  outfil = mike_getfil('arc_fil', SUBFIL=raw_fil, CHKFIL=chkfil, /name)
  if CHKFIL NE 0 AND not keyword_set( CLOBBER ) then begin
      print, 'mike_procarc: Arc ', outfil, ' exists.  Returning..'
      return, outfil
  endif

  ;; Determine xy_off with the FFT
  if arg_present( XYOFF ) AND not keyword_set( ALIGN ) then begin
      stop ;; Probably shouldnt be here
      xyoff = mike_arcxyoff_work( templ_fil, raw_fil, side )
  endif

  ;; Process
  rslt = mike_proc_sngl(raw_fil, side, FLATFIL=flatfil, CLOBBER=clobber, $
                        REDOOV=clobber, /ARC, IOUTFIL=outfil, GAIN=gain, $
                        READNO=readno, SETUP=setup, _EXTRA=extra)

  if keyword_set( ALIGN ) then begin
      ;; Grab the order structure
      ordr_str = mike_getfil('ordr_str', setup, SIDE=side)
      ;; Open the images
      arc = xmrdfits(outfil, 0, /silent)
      arci = xmrdfits(outfil, 1, /silent)
      ;; Here we go
      if keyword_set(chk) then window, 0, title='mike_arcalign'
      xyoff = mike_arcalign_work(arc, arci, ordr_str, $
                                 CHK=keyword_set(chk), FITPRM=fitprm)
  endif
      
      

  ;; ALL DONE
  print, 'mike_procarc_sngl: All Done! '
  return, outfil

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_procarc, mike, setup, obj_id, side, ATEMPL=atempl, $
                  CLOBBER=clobber, FLATFIL=flatfil, ALIGN=align, $
                  FLIP=flip

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'mike_procarc, mike, setup, obj_id, [side], FLATFIL=, /CLOBBER ' + $
        'ATEMPL= [v1.1]'
      return
  endif 
  
  ;; QA
  if setup GE 10 then c_s = strtrim(setup,2) else c_s = '0'+strtrim(setup,2) 
  fqa = 'QA/Arcs'+c_s
  a = findfile(fqa, count=count)
  if count EQ 0 then file_mkdir, fqa

;  Optional Keywords
  if not keyword_set( SIDE ) then side = [1L,2L]
  if not keyword_set( ATEMPL ) then atempl = 0L
  
; Loop on side
  nside = n_elements(side)
  for ii=0L,nside-1 do begin

      qq = side[ii]
      ;; SIDE
      if qq EQ 1 then print, 'mike_procarc: Processing BLUE arc(s)' $
      else print, 'mike_procarc: Processing RED arc(s)'

      ;; Grab all Arc files
      arcs = where(mike.setup EQ setup AND mike.flg_anly NE 0 AND $
                   mike.side EQ side[ii] AND $
                   strtrim(mike.type,2) EQ 'ARC', narc)
      if narc EQ 0 then begin
          print, 'mike_procarc: No Arcs found! Returning' 
          return
      endif

      ;; Grab all obj files
      indx = where(mike.flg_anly NE 0 AND mike.side EQ qq AND $
                   mike.obj_id EQ obj_id AND mike.setup EQ setup AND $
                   (strtrim(mike.type,2) EQ 'OBJ' OR $
                    strtrim(mike.type,2) EQ 'STD'), nindx)
      if nindx EQ 0 then begin
          print, 'mike_procarc: No Obj found!  Returning' 
          continue
      endif

      ;; LOOP
      for kk=0L,nindx-1 do begin

          print, 'mike_procarc: Index ', strtrim(kk,2)

          ;; Identify closest image in time
          diff = abs( (mike[indx[kk]].date + $
                                   mike[indx[kk]].exp/2./3600./24.) - $
                                  (mike[arcs].date + mike[arcs].exp/2./3600./24.) )
          mn = min(diff, gdarc)

;;          diff = abs( (mike[indx[kk]].ut + mike[indx[kk]].exp/2.) - $
;;            (mike[arcs].ut + mike[arcs].exp/2.) )
;;         a = where( diff GT 21.*3600, na)
;;         if na NE 0 then diff[a] = abs(diff[a] - 24.*3600.)
;;         mn = min(diff, gdarc)

          ;; Check for prior image
          outfil = mike_getfil('arc_fil', SUBFIL=mike[arcs[gdarc]].img_root, $
                               /name, CHKFIL=chkf)
          if CHKF NE 0 AND not keyword_set( CLOBBER ) then begin
              print, 'mike_procarc: Arc ', outfil, ' exists.  ' + $
                'Updating using previous offsets and then returning'
              mike[indx[kk]].arc_fil = outfil
              mike[indx[kk]].arc_xyoff = mike[arcs[gdarc]].arc_xyoff
              continue
          endif

          ;; Determine xy_off
          mike_arcxyoff, mike, qq, arcs[atempl], arcs[gdarc], xyoff

          ;; Process
          mike_proc, mike, arcs[gdarc], FLATFIL=flatfil, CLOBBER=clobber, $
            REDOOV=clobber, /ARC, IOUTFIL=outfil, FLIP=flip

          if keyword_set( ALIGN ) then begin
             ;; Grab the order structure
             ordr_str = mike_getfil('ordr_str', setup, SIDE=side)
             ;; Open the images
             arc = xmrdfits(outfil, 0, /silent)
             arci = xmrdfits(outfil, 1, /silent)
             ;; Here we go
             xyoff = mike_arcalign_work(arc, arci, ordr_str, FITPRM=fitprm, /CHK)
;             mike[arcs[gdarc]].arc_xyoff = fitprm
             mike[indx[kk]].arc_xyoff = fitprm
           endif
      
          ;; Structure
          mike[indx[kk]].arc_fil = outfil

      endfor
  endfor

  ;; ALL DONE
  print, 'mike_procarc: All Done! '
  return

end

