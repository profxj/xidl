;+ 
; NAME:
; esi_echcalibstd   
;   Version 1.1
;
; PURPOSE:
;    Create a sensitivity function given a standard star and its
;    appropriate calibration file.  
;
; CALLING SEQUENCE:
; esi_echcalibstd, esi, indx, HSTFIL=, CHKFIT=, ESOFIL=
;   
; INPUTS:
;   esi   -  MIKE structure
;   indx   -  Index of the Standard star in the MIKE structure to process
;
; RETURNS:
;
; OUTPUTS:
;  Sensitivity file in the 'Extract/' directory
;
; OPTIONAL KEYWORDS:
;  OUTFIL - Name of the output file (default:
;           'Extract/Sens_{esiroot}.fits' )
;  HSTFIL=  -- CALSPEC fits file for calibration
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_echcalibstd, esi, 1, ESOFIL='fhr4469.dat'
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_echfitstd
;
; REVISION HISTORY:
;   04-Jan-2004 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echcalibstd, esi, indx, HSTFIL=hstfil, CHKFIT=chkfit, $
                     BSPLIN=bsplin, SWV=swv, SFX=sfx, ESOFIL=esofil, $
                     BSET=bset, YFIT=yfit, SENS=sens, EVERYN=everyn, $
                     OBJFIL=objfil, indx2=indx2, INTER=INTER, NCOEFF=ncoeff, $
                     SV_ALLFIT=sv_allfit, ABMAG=abmag, OTHER=other, $
                     INSAVE = INSAVE

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echcalibstd, esi, indx, HSTFIL=, FXFIL=, NCOEFF=, /ABMAG, OTHER= [v1.1]'
    return
  endif 

;  Optional Keywords
  if n_elements(indx) NE 1 then stop
;  if not keyword_set( EVERYN ) then everyn = 36/esi[indx].rowbin
  if not keyword_set( OBJFIL ) then objfil = esi[indx].obj_fil
  if not keyword_set( OUTFIL ) then $
    outfil = esi_getfil('sens_fil', SUBFIL=esi[indx].img_root, /name)
  IF not KEYWORD_SET(SAVEFIL) THEN $
    savefil =  esi_getfil('sens_fil_save', SUBFIL = esi[indx].img_root, /name)
  if not keyword_set( NCOEFF ) then ncoeff = 7

; HST
  ;; Read file
  if keyword_set( HSTFIL ) then begin
     if strlen(HSTFIL) LT 30 then $
         hstfil = getenv('XIDL_DIR')+'/Spec/Flux/HSTFIL/'+hstfil
         readcol, hstfil, swv, sfx, format = 'D,D'
  endif

  if keyword_set( ESOFIL ) then begin
      if strlen(ESOFIL) LT 30 then $
        esofil = getenv('XIDL_DIR')+'/Spec/Flux/ESOFIL/'+esofil
        readcol, esofil, swv, sfx, format = 'D,D'
;
;     let's keep it at 1d-16
;      sfx = sfx*1d-16
  endif

  if keyword_set(OTHER) then begin
      if strlen(OTHER) LT 30 then $
        othfil = getenv('XIDL_DIR')+'/Spec/Flux/OTHER/'+other
      readcol, othfil, swv, sfx, format='F,F'
      ;; AB Mag?
      if keyword_set(ABMAG) then begin
          c = x_constants()
          fnu = 10.^((sfx+48.6)/(-2.5))
          flam = fnu * c.c / (swv)^2 / 1e-8  ; erg/s/cm^2/Ang
          sfx = flam * 1e16
      endif
  endif
          

  objstr = xmrdfits(objfil, 1, /silent)
  nordr = n_elements(objstr)
  objstr = objstr[0:9<nordr]

  ;; Normalize exposure
  objstr.box_fx = objstr.box_fx / esi[indx].exp

  ;; Call the fit program
  objstr.order = 15-objstr.order
  ;; If the standard does not go red enough, extrapolate it as a power law
  maxwv = max(objstr[nordr-1L].WAVE)
  IF maxwv GT max(swv) THEN BEGIN
      ipix = WHERE(swv GT 5000.0 AND sfx GT 0.0)
      dswv = djs_median((abs(swv - shift(swv, 1)))[ipix])
      nadd = ceil((maxwv - max(swv))/dswv)
      swv_add = max(swv) + (dindgen(nadd) + 1.0d)*dswv
      ;; fit 
      sfx_coeff = ladfit(alog10(swv[ipix]), alog10(sfx[ipix]))
      logsfx_add   = poly(alog10(swv_add), sfx_coeff)
      sfx_add = 10.0d^logsfx_add
      swv = [swv, swv_add]
      sfx = [sfx, sfx_add]
  ENDIF
  x_echfitstd, objstr, swv, sfx, outfil, savefil = savefil, INSAVE = INSAVE

  print, 'esi_echcalibstd:  Writing ', outfil


  print, 'esi_echcalibstd: All Done!'
  return
end
