;+ 
; NAME:
; esi_echtrcholes   
;     Version 1.1
;
; PURPOSE:
;    Traces the pinhole images in each order to create a curvature
;    map useful for data reduction.  
;
; CALLING SEQUENCE:
;   
;  esi_echtrcholes, esi, IMG=, XERR=, FITFIL=, FIT=,
;        GDENDS=, HOLEFIL=
;
; INPUTS:
;   esi     -  ESI structure
;
; RETURNS:
;
; OUTPUTS:
;  fitfil -- Name of IDL fit to pinhole traces (default:
;            Maps/hole_fit.idl)
;  holefil --  Name of combined pinhole image  (default:
;             Maps/img_hole.fits)
;
; OPTIONAL KEYWORDS:
;  IMG -- Name of pinhole image file  (default is to combine all
;         pinhole frames)
;  CLOBBER -- Overwrite existing map file
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   1) The current version works only with data taken from April 2000 -
;   present.  
;   2) Only setup for standard 1x1 binning
;
; EXAMPLES:
;   esi_echtrcholes, esi
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   07-Aug-2002 Written by JXP
;   01-Feb-2003 Polished by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echtrcholes, esi, IMG=img, FITFIL=fitfil, HOLEFIL=holefil, $
                     CLOBBER=clobber

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'esi_echtrcholes, esi, [xcen_mid], IMG=, XERR=, GDENDS=, FITFIL=  [v1.0]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set(REFMID) then refmid = 2048L
  if not keyword_set(REFTOP) then reftop = 4040L
  if not keyword_set(REFBOT) then refbot = 55L
  if not keyword_set( FITFIL ) then fitfil = 'Maps/hole_fit.idl'
  if not keyword_set( HOLEFIL ) then holefil = 'Maps/img_hole.fits'

; Check for fitfil
  a = findfile(fitfil, count=na)
  if na NE 0 and not keyword_set( CLOBBER) then begin
      print, 'esi_echtrcholes: Hole fit file exists! Returning... ', fitfil
      return
  endif
;;;;;; GET HOLE IMG ;;;;;;;;;;
  if not keyword_set( IMG ) then begin
      ;; Hole frames
      holes = where(esi.slit EQ 9.99 AND esi.mode EQ 2 AND $
                    esi.flg_anly NE 0, nhole)
      if nhole EQ 0 then begin
          print, 'esi_echtrcholes: No Hole images!'
          return
      endif

      ;; Bias Subtract
      esi_subbias, esi, holes, /force

      ;; Combine
      if nhole GT 1 then begin
          print, 'esi_echtrcholes: Combining hole images'
          xcombine, 'OV/ov_'+esi[holes].img_root, img_hole, head, $
            FCOMB=2, SCALE=esi[holes].exp, $
            GAIN=esi[holes[0]].gain, RN=esi[holes[0]].readno
      endif else img_hole = mrdfits('OV/ov_'+esi[holes[0]].img_root,/silent)
      ;; WRITE
      mwrfits, img_hole, holefil, /create, /silent
      esi_delov, esi, holes
  endif else img_hole = x_readimg(img)

; IVAR
  ivar = 1./ (img_hole > 1)
  ;; Avoid bad columns
  ivar[422:423,2651:4095] = -1.
  ivar[425,2648:4095] = -1.
  ivar[432:437,2648:4095] = -1.
  ;; Set flux to 0.
  img_hole[422:423,2651:4095] = 0.
  img_hole[425,2648:4095] = 0.
  img_hole[432:437,2648:4095] = 0.

; Choose lines to trace
  smsh = djs_median(img_hole[*,refmid-3:refmid+3],2)
  x_fndpeaks, smsh, center, NSIG=10.
  xstart = center
  ystrt = replicate(refmid, n_elements(xstart))
  for j=0L,9 do $
    xstart = trace_fweight(img_hole, xstart, ystrt, radius=3.5, $
                           INVVAR=ivar)

  ;; CHK
  if keyword_set( CHK ) then $
    x_splot, smsh, XTWO=xstart, YTWO=smsh[round(xstart)], PSYM_Y2=1, /block

; Trace
  print, 'esi_echtrcholes: Tracing mid with trace_crude...'
;  restore, 'tmp.idl'
  xcen_mid = trace_crude(img_hole, ivar, yset=ycen_pos, XSTART=xstart, $
                         radius=3.5, ystart=refmid, xerr=xerr_mid, $
                         MAXSHIFTE=0.8, NMED=5, NAVE=5)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Edit top
  smsh = djs_median(img_hole[*,reftop-2:reftop+2],2)
  x_fndpeaks, smsh, center, NSIG=7.
  lft = where(center LT 419)
  xstart = center[lft]

  ;; CHK
  if keyword_set( CHK ) then $
    x_splot, smsh, XTWO=xstart, YTWO=smsh[round(xstart)], PSYM_Y2=1, /block

; Trace
  print, 'esi_echtrcholes: Tracing top with trace_crude...'
  xcen_top = trace_crude(img_hole, ivar, yset=ytop_pos, XSTART=xstart, $
                         radius=3.5, ystart=reftop, xerr=xerr_top, $
                         MAXSHIFTE=0.8, NMED=5, NAVE=5)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PATCH Top
;;;;;;;;;;;;;;

  print, 'esi_echtrcholes: Patching'

  ;; Defect
  for i=0L,2 do xcen_mid[3775:4095,7+i] = xcen_top[3775:4095,i]
  xerr_mid[3860:3920,7:9] = 99.99
  
  ;; Bad column(s)
  for i=6L,16 do begin
      for j=4095L,0,-1 do begin
          if(xcen_top[j,i] GT 417.) then begin
              xcen_mid[j:4095,7+i] = xcen_top[j:4095,i]
              xerr_mid[j:4095,7+i] = xerr_top[j:4095,i]
              xerr_mid[j-250:j-1,7+i] = 99.99
              break
          endif
      endfor
  endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Parse out bad lines!

  print, 'esi_echtrcholes: Parsing bad ends!'
  sz_xcen = size(xcen_mid, /dimensions)
  ntrc = sz_xcen[1]
  gdends = lonarr(ntrc,2)
  gdends[*,1] = 4095L
  ;; Tops
  for q=25L,ntrc-1 do begin
      dx = shift(xcen_mid[2800L:sz_xcen[0]-1,q],-30) - $
        xcen_mid[2800L:sz_xcen[0]-1,q] 
      a = where(dx GE -2, na)
      if (a[0]+2800L) GT 4050L then gdends[q,1] = 4095L $
      else gdends[q,1] = a[0] + 2700L
  endfor
  ;; Bottom
  duml = lindgen(1201)
  for q=0L,30 do begin
      dx = shift(xcen_mid[0L:1200L,q],-30) - xcen_mid[0L:1200L,q] 
      a = where(dx LT 0.1 AND duml LT 1165L, na)
      if na EQ 0 then gdends[q,0] = 230L else begin
          if a[na-1] LT 30L then gdends[q,0] = 0L $
          else gdends[q,0] = a[na-1] + 100L
      endelse
  endfor
  ;; Bottom
  duml = lindgen(901)
  for q=31L,ntrc-1 do begin
      dx = shift(xcen_mid[0L:900,q],-30) - xcen_mid[0L:900,q] 
      a = where(dx LT 0.1 AND duml LT 865L, na)
      if na EQ 0 then gdends[q,0] = 0L else begin
          if a[na-1] LT 30L then gdends[q,0] = 0L $
          else gdends[q,0] = a[na-1] + 100L
      endelse
  endfor

;;;;;;;;;;;;;;
;  DEAL with Blending Holes

  ;; 54
  indx = ntrc-37L
  y_end = 3100L
  sep1 = xcen_mid[y_end,indx] - xcen_mid[y_end,indx-1]
  sep2 = xcen_mid[y_end,indx] - xcen_mid[y_end,indx-2]
  frac = sep1/sep2
  xcen_mid[y_end+1:sz_xcen[0]-1,indx]= $
    (xcen_mid[y_end+1:sz_xcen[0]-1,indx-1] - $
     frac*xcen_mid[y_end+1:sz_xcen[0]-1,indx-2])/$
    (1.-frac)
  ;; 55
  indx = ntrc-36L
  y_end = 3100L
  sep1 = xcen_mid[y_end,indx+1] - xcen_mid[y_end,indx]
  sep2 = xcen_mid[y_end,indx+2] - xcen_mid[y_end,indx]
  frac = sep1/sep2
  xcen_mid[y_end+1:sz_xcen[0]-1,indx]=(xcen_mid[y_end+1:sz_xcen[0]-1,indx+1] - $
                                       frac*xcen_mid[y_end+1:sz_xcen[0]-1,indx+2])/$
    (1.-frac)
  ;; 63
  indx = ntrc-28L
  y_end = 2700L
  sep1 = xcen_mid[y_end,indx] - xcen_mid[y_end,indx-1]
  sep2 = xcen_mid[y_end,indx] - xcen_mid[y_end,indx-2]
  frac = sep1/sep2
  xcen_mid[y_end+1:sz_xcen[0]-1,indx]=(xcen_mid[y_end+1:sz_xcen[0]-1,indx-1] - $
                                       frac*xcen_mid[y_end+1:sz_xcen[0]-1,indx-2])/$
    (1.-frac)
  ;; 64
  indx = ntrc-27L
  y_end = 2700L
  sep1 = xcen_mid[y_end,indx+1] - xcen_mid[y_end,indx]
  sep2 = xcen_mid[y_end,indx+2] - xcen_mid[y_end,indx]
  frac = sep1/sep2
  xcen_mid[y_end+1:sz_xcen[0]-1,indx]=(xcen_mid[y_end+1:sz_xcen[0]-1,indx+1] - $
                                       frac*xcen_mid[y_end+1:sz_xcen[0]-1,indx+2])/$
    (1.-frac)
  ;; 72
  indx = ntrc-19L
  y_end = 3000L
  sep1 = xcen_mid[y_end,indx] - xcen_mid[y_end,indx-1]
  sep2 = xcen_mid[y_end,indx] - xcen_mid[y_end,indx-2]
  frac = sep1/sep2
  xcen_mid[y_end+1:sz_xcen[0]-1,indx]=(xcen_mid[y_end+1:sz_xcen[0]-1,indx-1] - $
                                       frac*xcen_mid[y_end+1:sz_xcen[0]-1,indx-2])/$
    (1.-frac)
  ;; 73
  indx = ntrc-18L
  y_end = 3000L
  sep1 = xcen_mid[y_end,indx+1] - xcen_mid[y_end,indx]
  sep2 = xcen_mid[y_end,indx+2] - xcen_mid[y_end,indx]
  frac = sep1/sep2
  xcen_mid[y_end+1:sz_xcen[0]-1,indx]=(xcen_mid[y_end+1:sz_xcen[0]-1,indx+1] - $
                                       frac*xcen_mid[y_end+1:sz_xcen[0]-1,indx+2])/$
    (1.-frac)

;;;;;;;;;;;;;
;  FIT

  tmp_fit = { fitstrct }
  tmp_fit.func = 'POLY'
  tmp_fit.nord = 6
  tmp_fit.hsig = 3.
  tmp_fit.lsig = 3.
  tmp_fit.niter = 2L
  tmp_fit.maxrej = 100L
  tmp_fit.flg_rej = 1

  fin_fit = replicate(tmp_fit, ntrc)
  msk = lindgen(sz_xcen[0])

  clr = getcolor(/load)
  for q=0,ntrc-1 do begin
      ;; Create mask
      msk[*] = 0L
      msk[gdends[q,0]:gdends[q,1]] = 1L
      a = where(xerr_mid[*,q] GT 0.1, na)
      if na NE 0 then msk[a] = 0L
      b = where(msk EQ 1, nb)
      ;; Fit
      if nb LT 1000L then tmp_fit.nord = 3 else tmp_fit.nord=5
      fit = x_fitrej(findgen(sz_xcen[0]), xcen_mid[*,q], $
                     MSK=msk, FITSTR=tmp_fit)
      ;; SV
      xcen_mid[*,q] = fit
      ;; Print rms
      print, 'esi_echtrcholes: ', q, ' RMS = ', tmp_fit.rms
      ;; Save
      fin_fit[q] = tmp_fit
  endfor

  ;; WRITE
  save, fin_fit, xcen_mid, filename=fitfil

  print, 'esi_echtrcholes: All done! Fit file is ', fitfil
  return
end
