;+ 
; NAME:
; lris_calibstd   
;   Version 1.1
;
; PURPOSE:
;    CURRENTLY UNDER CONSTRUCTION
;
; CALLING SEQUENCE:
;   
;   spec = x_apall(ydat, [head])
;
; INPUTS:
;   ydat       - Values 
;   [head]     - Header
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   wave       - wavelength array
;   DISPLAY    - Display the sky subtracted image with xatv
;   OVR        - String array for ov region:  '[2050:2100, *]'
;   ERROR      - Variance array
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   kast_calibstd, kast, 0, 1, 0
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   03-Mar-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro lris_calibstd, kast, indx, outfil, HSTFIL=hstfil

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'kast_calibstd, kast, indx, HSTFIL=  [v1.0]'
    return
  endif 

;  Optional Keywords
  if not keyword_set( EVERYN ) then everyn = 5

; Open objfil
  objfil = kast[indx[0]].obj_fil 
  if x_chkfil(objfil+'*') EQ 0 then begin
      print, 'kast_calibstd: No Obj file ', objfil
      return
  endif
  objstr = xmrdfits(objfil, 1, STRUCTYP='specobjstrct', /silent)
  npix = objstr[0].npix

; HST
  if keyword_set(HSTFIL) then begin
      ;; Read file
      hst = xmrdfits(hstfil, 1, /silent)

      ;; BSpline
      bset = bspline_iterfit(hst.wavelength,hst.flux, yfit=yfit, everyn=everyn)

      ;; Calculate at wavelength
      sens = bspline_valu(objstr[0].wave[0:npix-1], bset) / $
        (objstr[0].fx[0:npix-1] / kast[indx[0]].exp)

      ;; BSpline
      bset = bspline_iterfit(objstr[0].wave[0:npix-1], sens, $
                             yfit=yfit, everyn=everyn)

      x_splot, sens, ytwo=yfit, /block

      ;; Output
      mwrfits, bset, outfil, /create
      
  endif

  print, 'kast_calibstd: All Done!'
  return
end
