;+ 
; NAME:
; mike_getwvarr   
;     Version 1.1
;
; PURPOSE:
;    Extract flux from 2D image to create ten 1D spectra (1 per order)
;    Output is written to the object structure (e.g. Extract/Obj_mike0024.fits)
;    The code only does boxcar extraction for now.
;
; CALLING SEQUENCE:
;   
;  mike_boxextrct, mike, obj_id, [exp], /DEBUG, /CHK, /STD, APER=,
;  RADIUS=
;
; INPUTS:
;   mike   -  ESI structure
;   indx  -  Indices of objects to process
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   /CHK    - Show final trace
;   /STD    - Extraction should be set for a standard star
;   /DEBUG  - Stop within extraction routine to check stuff
;   APER=   - Set aperture by hand (e.g. [5., 7.] )
;   RADIUS= - Size of window for setting aperture size (default: 20L)
;   ORDRS=  - Orders to extract (default: [0L,9L])
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;  1)  The program begins extracting in order 9L (physical 6) and will
;  automatically calculate an aperture for that order.  If there is
;  insufficient flux in the following orders, it will adopt the value
;  from the next higher order. 
;
; EXAMPLES:
;   mike_getwvarr, mike, 1L, [0L]
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   02-Jan-2004 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_getwvarr, ywv, iarc, nrow, side, WVARR=wvarr, NPIX=npix, YVAL=yval
;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'mike_getwvarr, ywv, iarc, nrow, side , ' + $
        'WVARR=, NPIX=, YVAL=  [v1.0]'
      return
  endif 
  
;  Optional Keywords

  bin = round(4096./nrow)

;  Spline
  wave_spline = spl_init(ywv, iarc, /double)

;  Set full wave array
  spl = 299792.458d
  if side EQ 1 then begin
      fullwv = 10^(alog10(3000.d) + $
                   dindgen(150000L/bin)*alog10(1.+ bin*1.56/spl)) 
      ;; Endpoints
      wvb = 10^spl_interp(ywv, iarc, wave_spline, 0., /double)
      a = where(fullwv GT wvb)
      i1 = a[0]
      wve = 10^spl_interp(ywv, iarc, wave_spline, nrow-1., /double)
      a = where(fullwv LT wve,na)
      i2 = a[na-1]
      
  endif else begin
      fullwv = 10^(alog10(3000.d) + $
                   dindgen(150000L/bin)*alog10(1.+ bin*2.24/spl))
      ;; Endpoints
      wvb = 10^spl_interp(ywv, iarc, wave_spline, nrow-1., /double)
      a = where(fullwv GT wvb)
      i1 = a[0]
      wve = 10^spl_interp(ywv, iarc, wave_spline, 0., /double)
      a = where(fullwv LT wve,na)
      i2 = a[na-1]

      ;; Neg
      neg = 1
  endelse

  ;; Wave array
  wvarr = fullwv[i1:i2]
  npix = i2-i1+1

; Find y values
  yval = dblarr(npix)
  for q=0L,npix-1 do begin
      yval[q] = x_fndspln(ywv, iarc, alog10(wvarr[q]), wave_spline, TOLER=1d-5)
  endfor

;  DONE
;  print, 'mike_boxextrct: All done! '
  return
end
