;+ 
; NAME:
; x_tracearc   
;    Version 1.0
;
; PURPOSE:
;    Traces a series of arc lines
;
; CALLING SEQUENCE:
;   
;   x_tracearc, img, [ymnx], STRTY=, /ROT
;
; INPUTS:
;   img       - Input arc image
;   [ymnx]    - Slit edges (default to image edges)
;
; RETURNS:
;   tracestr
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   YSTRT - COLUMN/ROW to start at
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_tracearc, img, map
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   28-Jan-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_tracearc, img, ymnx, VAR=var, GAIN=gain, RN=rn, $
                     ystrt=ystrt, ROT=rot, NSIG=nsig, RADIUS=radius, $
                     MIN_XERR=min_xerr, SILENT=silent, INITROW=initrow



;  Error catching
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x_tracearc( img, [ymnx], VAR=, GAIN=, RN=, '
    print,   '           ystrt=, /ROT, MIN_XERR=, /SILENT) [v1.0]'
    return, -1
  endif 

  sz = size(img)
  if sz[0] NE 2 then begin
      print, 'x_tracearc: need a 2-D image'
      return, -1
  endif

;  Optional Keywords

  if not keyword_set( YMNX ) then begin
      ymnx = lonarr(2) 
      ymnx[0] = 0
      ymnx[1] = sz[2]-1
  endif

  if not keyword_set( MIN_XERR ) then min_xerr = 0.01
  if not keyword_set( NSIG ) then nsig = 5.
  if not keyword_set( RADIUS ) then radius = 1.5
  if not keyword_set( YSTRT ) then YSTRT = round(total(ymnx)/2.)

  ; Offset ystrt to account for subimg
  ystrt = ystrt - ymnx[0]

  ; Set tmpimg and multiply in gain
  if keyword_set( GAIN ) then tmpimg = img[*,ymnx[0]:ymnx[1]]*gain $
  else tmpimg=img[*,ymnx[0]:ymnx[1]]

  if not keyword_set( RN ) then rn = 6.  ; Readnoise
  if keyword_set( ROT ) then begin
      print, 'x_tracearc: Not set up for this'
      stop
      tmpimg = transpose(tmpimg) 
  endif

;  Create inv variance as necessary

  if not keyword_set( VAR ) then tmpvar = (tmpimg + (tmpimg EQ 0) + RN^2) $
  else tmpvar = var[*,ymnx[0]:ymnx[1]]
  

;;;;;;;;;;;;;
;  ivar and initial peaks

  sigimg = sqrt( tmpvar > 0 )
  ivar = 1./tmpvar
  ; Set bad pixels to 0 in ivar
  ivar[where(tmpimg LE 0)] = 0.

  ; Take median around ystrt (5 rows)

  if not keyword_set( INITROW ) then $
    initrow = djs_median(tmpimg[*,ystrt-2:ystrt+2], 2)

  x_fndpeaks, initrow, center, NSIG=nsig, SILENT=silent

  ; Center on the peaks with fweight (twice)
  
  xdat = findgen(sz[1])
  npk = n_elements(center)
  xstrt = fltarr(npk)
  for i=0L,npk-1 do begin
      xstrt[i] = trace_fweight(tmpimg, center[i], ystrt, $
                               radius=1.5, invvar=ivar)
      for j=0L,10-1 do $
        xstrt[i] = trace_fweight(tmpimg, xstrt[i], ystrt, $
                                 radius=3., invvar=ivar)
      for j=0L,10-1 do $
        xstrt[i] = trace_fweight(tmpimg, xstrt[i], ystrt, $
                                 radius=radius, invvar=ivar)
  endfor
      

;;;;;;;;;;;;;
;  Trace crude 
  if not keyword_set( SILENT ) then print, 'x_tracearc: Tracing'

  xcen = trace_crude(tmpimg, ivar, yset=ycen, xstart=xstrt, $
                     radius=radius, ystart=YSTRT, xerr=xerr, nave=3, $
                     maxshift0=1.0)
;;;;;;;;;;;;;;;;
;  Fit poly nomials to each line
  
  xfitstr = { fitstrct }
  xfitstr.func = 'POLY'
  xfitstr.nord = 2
  xfitstr.lsig = 3.
  xfitstr.hsig = 3.
  xfitstr.niter = 3
  xfitstr.minpt = 5
  xfitstr.maxrej = 5
  
  sz_trc = size(xcen, /dimensions)
  xfit = xcen
  for i=0L,sz_trc[1]-1 do $
      xfit[*,i] = x_fitrej(ycen[*,i], xcen[*,i], $
                           SIG=(xerr[*,i]>min_xerr), FITSTR=xfitstr)

;;;;;;;;;;;;;;
; Create structure and transpose x,y by the variable name

  tracestrct = { $
                 xstrt: long(ystrt), $
                 leny: sz[2], $
                 xcen: temporary(xcen), $
                 xfit: temporary(xfit), $
                 xerr: temporary(xerr), $
                 ycen: temporary(ycen) $
;                 msk_cen: temporary(msk_cen) $
               }

; Free memory
  delvarx, tmpimg, tmpvar, ivar

  return, tracestrct
end
