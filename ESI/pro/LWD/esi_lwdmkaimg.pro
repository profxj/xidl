;+ 
; NAME:
; esi_lwdmkaimg   
;     Version 1.0
;
; PURPOSE:
;    Trace an Arc Image
;
; CALLING SEQUENCE:
;   
;  esi_lwdmkaimg, esi, slit
;
; INPUTS:
;   esi   -  ESI structure
;   slit  -  Slit size
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_lwdmkaimg, esi, slit
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro esi_lwdmkaimg, esi, slit, LINTRC=lintrc, SZ_ARC=sz_arc, CHK=chk

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_lwdmkaimg, esi, slit, LINLIST=, /INTER [v1.0]'
      return
  endif 

;  Optional Keywords
  
  if not keyword_set( SZ_ARC ) then sz_arc = [2232L, 3300L]

; Get TRC Structure
  trcfil = 'Arcs/ATRC_LWD'+esi_slitnm(slit)+'.fits'
  a = findfile(trcfil, count=na)
  if na EQ 0 then begin
      print, 'esi_lwdmkaimg: No Arc TRC found!'
      stop
  endif
  trcstr = mrdfits(trcfil, 1, /silent)
  ntrc = n_elements(trcstr.lintrc)
  sz_trc = size(trcstr.xcen, /dimensions)


; Setup Fit Structure
  fitstr = { fit2dstrct }
  fitstr.func = 'POLY'
  if keyword_set( NX ) then fitstr.nx = nx else fitstr.nx = 8
  if keyword_set( NY ) then fitstr.ny = ny else fitstr.ny = 8
  fitstr.niter = 1
  fitstr.lsig = 3.
  fitstr.hsig = 3.
  fitstr.flg_rej = 1


; Setup Arrays
  ycen = lindgen(sz_arc[1]) # replicate(1., ntrc) 
  msk = bytarr(sz_trc[0], sz_trc[1])
  for i=0,ntrc-1 do msk[trcstr.ends[i,0]:trcstr.ends[i,1], i] = 1
  gd_cen = where(msk EQ 1, ngd)
  xydat = fltarr(ngd,2)
  xydat[*,0] = trcstr.xcen[gd_cen]
  xydat[*,1] = ycen[gd_cen]

  wvdat = float(replicate(1, sz_arc[1]) # trcstr.lintrc)
  
; Fit
  if not keyword_set( SILENT ) then $
    print, 'esi_lwdmkaimg: Fitting in 2D with SVDFIT  ', systime()
  fit = x_fit2dsurf(xydat, wvdat[gd_cen], $
                    FITSTR=fitstr, /svdft)
  if arg_present( RES ) then res = fit - wvdat[gd_cen]
  print, 'esi_lwdmkaimg: RMS = ', fitstr.rms
  
; Map over the entire image
  
  adum = dindgen(sz_arc[0]) # replicate(1., sz_arc[1])
  adum2 = replicate(1., sz_arc[0]) # dindgen(sz_arc[1])
  
  ;; Funny format
  xydat = dblarr(sz_arc[0]*sz_arc[1], 2)
  xydat[*,0] = adum
  xydat[*,1] = adum2
  delvarx, adum, adum2

  ;; Mapping
  if not keyword_set( SILENT ) then $
    print, 'esi_lwdkmaimg: Mapping over the whole image  ', systime()
  tot_fit = x_calc2dfit(xydat, FITSTR=fitstr)
  map = reform(tot_fit,sz_arc[0], sz_arc[1])
  delvarx, tot_fit

; Zero out bad pix
  print, 'esi_lwdkmaimg: Elminating uncalibrated pixels  ', systime()
	  a = where(abs(trcstr.lintrc - 4046.) LT 2., na)
  if na EQ 0 then begin
      print, 'esi_lwdmkaimg: No 4046 line!'
      stop
  endif
  msk = bytarr(sz_arc[0], sz_arc[1])
  msk[*,0:230L] = 1
  msk[*,3075:sz_arc[1]-1] = 1
  ;; Trace along the 4046 line
  for i=230L,3075L do begin
      edg = (round(trcstr.xcen[i,a[0]])-100L) > (-1L)
      if edg GT -1 then msk[0:edg,i] = 1
  endfor
  a = where(msk EQ 1)
  map[a] = 0.
  fin_arc = float(map)
  delvarx, map
  
; Vacuum wavelengths
  print, 'esi_echmkaimg: Converting to vacuum wavelengths'
  a = where(fin_arc GT 0.)
  tmpaimg = fin_arc[a]
  airtovac, tmpaimg
  fin_arc[a] = temporary(tmpaimg)

; CHK

  if keyword_set( CHK ) then begin
      xatv, fin_arc, /block
      stop
  endif

; Output
  outfil = 'Arcs/AIMG_LWD'+esi_slitnm(slit)+'.fits'
  mwrfits, fin_arc, outfil, /create, /silent


; Cards
  objstd = where(esi.mode EQ 1 AND esi.flg_anly NE 0 AND $
                 esi.slit EQ slit AND $
                 (esi.type EQ 'STD' OR esi.type EQ 'OBJ'), nobj)
  if nobj NE 0 then esi[objstd].arc_fil = outfil

  ;; DONE
  print, 'esi_lwdmkaimg: All done!  ' 

  return
end
