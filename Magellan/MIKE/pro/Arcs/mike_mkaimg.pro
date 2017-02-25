;+ 
; NAME:
; mike_mkaimg
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
;   
;  mike_mkaimg, mike, setup, [side], /CHK, /CLOBBER
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
;  2D wavelength image with name like 'Arcs/Arc_mb0439I.fits'
;
; OPTIONAL KEYWORDS:
;   /CLOBBER  - Overwrite previous image
;   /CHK      - Display the final image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_mkaimg, mike, setup, obj_id
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  mike_getfil
;  mike_ordermask
;  mike_mkaimg_work
;
; REVISION HISTORY:
;   15-May-2003 Written by SB
;
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function mike_mkaimg_work, arc_fil, setup, side, CHK=chk, CLOBBER=clobber, $
                           BAD_ANLY=bad_anly, SHFTPRM=shftprm
;

  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'rslt = mike_mkaimg_work( arc_fil, setup, side, /CHK, /CLOBBER, SHFTPRM= ) [v1.1]'
      return, -1
  endif 

  ;; Optional keywords
  if not keyword_set( SHFTPRM ) then SHFTPRM = [0., 0.]
  bad_anly = 0
 
  ;; Check for outfil
  out_fil = mike_getfil('arc_img', subfil=arc_fil, /name, CHKFIL=chkf)
  if CHKF NE 0 then begin
      print, 'mike_mkaimg: File exists: ', out_fil
      if not keyword_set( CLOBBER ) then begin
          print, 'mike_mkaimg: Use /CLOBBER to overwrite.  Skipping...'
          return, out_fil
      endif else print, 'mike_mkaimg: Clobbering -- ', out_fil
  endif

  ;; Read in order structure
  ordr_str = mike_getfil('ordr_str', setup, SIDE=side)
  nordr = n_elements(ordr_str)

  ;; Grab file names
  arc2d_fil = mike_getfil('arc_2Dfit', subfil=arc_fil, /name)
  fil_fittrc = mike_getfil('arc_fittrc', subfil=arc_fil, /name, CHKFIL=chkf)

  ;; Check arc_fil for size
  head = xheadfits(arc_fil)
  sz = lonarr(2)
  sz[0] = sxpar(head, 'NAXIS1')
  sz[1] = sxpar(head, 'NAXIS2')
  ximage = lindgen(sz[0]) # replicate(1,sz[1])
  yimage = lindgen(sz[1]) ## replicate(1,sz[0])

  ;; Grab Slope fit
  mkaimg_str = mike_getfil('arc_fittrc', subfil=arc_fil)

  ;; Create Final image
  aimg = dblarr(sz[0],sz[1])

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; LOOPING

  ;; Offset the order edges by xoff
  ordr_shift = ordr_str
  shft = mike_shifti(shftprm, OSTR=ordr_shift)

  rslt = x_mkaimg(arc_fil, ordr_shift, arc2d_fil, fil_fittrc, $
                  out_fil, CHK=chk, CLOBBER=clobber, $
                  BAD_ANLY=bad_anly) 

  if bad_anly GT 0 then begin
    print, 'mike_mkaimg: Bad wavelength solution, not writing to disk'
    return, out_fil
  endif

  return, out_fil
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_mkaimg, mike, setup, obj_id, side, CHK=chk, CLOBBER=clobber

;

  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'mike_mkaimg, mike, setup, obj_id, [side], /CHK, /CLOBBER [v1.0]'
      return
  endif 

;  Optional Keywords
  if not keyword_set( SIDE ) then side = [1L,2L]

; Loop on side
  for kk=0L,n_elements(side)-1 do begin
      qq = side[kk]
      ;; SIDE
      if qq EQ 1 then print, 'mike_mkaimg: Making BLUE Arcimages' $
      else print, 'mike_mkaimg: Making RED Arcimages'

      ;; Grab all obj indices
      indx = where(mike.flg_anly NE 0 AND mike.side EQ qq AND $
                   mike.obj_id EQ obj_id AND mike.setup EQ setup AND $
                   (strtrim(mike.type,2) EQ 'OBJ' OR $
                   strtrim(mike.type,2) EQ 'STD'), nindx)
      if nindx EQ 0 then begin
          print, 'mike_mkaimg: No Obj found!  Returning' 
          continue
      endif

      arcfil_all= mike[indx].arc_fil
      asort = sort(arcfil_all)
      auniq = uniq(arcfil_all[asort])

      arc_roots = strmid(mike.arc_fil, max(rstrpos(mike.arc_fil, '.')-5), 5)

      ;; LOOP
      for mm=0L,n_elements(auniq)-1 do begin

          arc_fil = arcfil_all[asort[auniq[mm]]]
          idx = indx[asort[auniq[mm]]]
;          xoff = mike[idx].arc_xyoff[0]

          stop ;; Odds are that arc_xyoff is not set!
          rslt = mike_mkaimg_work(arc_fil, setup, qq,$
                                  CHK=chk, CLOBBER=clobber, $
                                  SHFTPRM=mike[idx].arc_xyoff, $
                                  BAD_ANLY=bad_anly) 

          if bad_anly GT 0 then begin
             remove = where(strpos(mike.img_root, arc_roots[idx]) GT -1)
             if remove[0] NE -1 then mike[remove].flg_anly = 0
          endif

          if size(rslt,/tname) NE 'STRING' then stop
      endfor
  endfor

  print, 'mike_mkaimg: All done!'

  return
end
