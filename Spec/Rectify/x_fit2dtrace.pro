;+ 
; NAME:
; x_fit2dtrace   
;    Version 1.1
;
; PURPOSE:
;    Fits the traced y-distortion of a flat in 2D given a set of
;    traces for at various spots in the image.
;
; CALLING SEQUENCE:
;   
;   map = x_fit2dtrace( trcstr, FITSTR=, RES=, NX=, NY=, /SILENT )
;
; INPUTS:
;   trcstr   - Trace structure
;
; RETURNS:
;   map - Map of the y-distortion
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  nx -  Number of orders along x axis [default: 4L]
;  ny -  Number of orders along y axis [default: 4L]
;
; OPTIONAL OUTPUTS:
;   FITSTR  - 2D fitting function structure of the trace
;   RES   - Residual image of the fit
;
; COMMENTS:
;   Need to watch for memory issues with 2d surfaces
;
; EXAMPLES:
;   map = x_fit2dtrace( trcstr )
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   07-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_fit2dtrace, trcstr, RES=res, FITSTR=fitstr, NX=nx, NY=ny, $
                       SILENT=silent



;  Error catching
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'map = x_fit2dtrace( trcstr, RES=, FITSTR=, NX=, NY= [v1.1]'
    return, -1
  endif 

;  Optional Keywords

  ; This is temporary
  if tag_exist( trcstr, 'xstrt' ) EQ 0 then xstrt = 400 else xstrt=trcstr.xstrt
  if tag_exist( trcstr, 'leny' ) EQ 0 then leny = 2046 else leny=trcstr.leny

  sz_cen = size(trcstr.xcen, /dimensions)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  2D Fit!

  if not keyword_set( SILENT ) then print, 'x_fit2dtrace: 2D Fitting!'

  ; Take mean of 11 columns at xstrt
  yxSTRT = replicate(1., sz_cen[0]) # $
    (total(trcstr.ycen[XSTRT-5:xSTRT+5,*],1)/11.) 

  ; Delta y as a function of x  
  dyfit = yXSTRT - trcstr.ycen
  delvarx, yxSTRT


; Reform the arrays for x_fit2dsurf

  gd_cen = where(trcstr.msk_cen EQ 1, ngd)

  xydat = fltarr(ngd,2)
  xydat[*,0] = trcstr.xcen[gd_cen]
  xydat[*,1] = trcstr.ycen[gd_cen]

  ; Fit structure
  fitstr = { fit2dstrct }
  fitstr.func = 'POLY'
  if keyword_set( NX ) then fitstr.nx = nx else fitstr.nx = 4
  if keyword_set( NY ) then fitstr.ny = ny else fitstr.ny = 4
  fitstr.niter = 1
  fitstr.lsig = 5.
  fitstr.hsig = 5.

  ; Fit
  if not keyword_set( SILENT ) then $
    print, 'x_fit2dtrace: Fitting in 2D with SVDFIT'
  xfit = x_fit2dsurf(xydat, dyfit[gd_cen], $
                     (trcstr.yerr[gd_cen]>0.05), FITSTR=fitstr, /svdft)

  if arg_present( RES ) then res = xfit - dyfit[gd_cen]
  delvarx, xfit, dyfit

  ; Map over the entire image
  
  adum = dindgen(sz_cen[0]) # replicate(1., leny)
  adum2 = replicate(1., sz_cen[0]) # dindgen(leny)
  
  xydat = dblarr(leny * sz_cen[0], 2)
  xydat[*,0] = adum
  xydat[*,1] = adum2
  delvarx, adum, adum2

  if not keyword_set( SILENT ) then $
    print, 'x_fit2dtrace: Mapping over the whole image'
  tot_fit = x_calc2dfit(xydat, FITSTR=fitstr)

  delvarx, xydat
  map = reform(tot_fit,sz_cen[0], leny)

  delvarx, tot_fit

  return, map
end
