; NAME:
; x_initesi
;    Version 1.0
;
; PURPOSE:
;    Initialize an instrument structure for ESI (echellette)
;
; CALLING SEQUENCE:
;  tmp = {dlastruct}
;
; INPUTS:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   20-Oct-2005 Written by JXP
;-
;------------------------------------------------------------------------------
pro x_initesi, keckinstr, SLIT=slit, INFIL=infil, STR_TEL=str_tel

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'x_initesi, slit,  STR_TEL= [v1.0]'
      return
  endif 

  ;; Telescope
  x_initkeck, str_tel
  str_tel.name = 'KeckII'

  ;; Instrument
  keckinstr = {instrstruct}
  keckinstr.MAG_PERP = 6.05  ; Modified to give 0.15" pixels
  keckinstr.MAG_PARA = 6.05  ; Same
  keckinstr.PIXEL_SIZE = 18.0 ; in microns

  keckinstr.SCALE_PERP = str_tel.PLATE_SCALE* $
                         keckinstr.MAG_PERP*(keckinstr.PIXEL_SIZE/1000.) ; Arcsec
  keckinstr.SCALE_PARA = str_tel.PLATE_SCALE* $
                         keckinstr.MAG_PARA*(keckinstr.PIXEL_SIZE/1000.) ; Arcsec
  keckinstr.R     = 26000    ; 1 pixel (native) dispersion
          
  ;; These could be modified
;  keckinstr.MLAMBDA = 356385.3016 ; = 2*sigma*sin(delta)*cos(theta)
                                ;   52.676 grooves/mm 
                                ;   70.43  degress blaze 
;  keckinstr.DELY    = 0.154693 ; = f2 * Ac 

  ;; Detector
  keckinstr.readno = 2.7
  keckinstr.dark = 0.0
  keckinstr.bind = 1
  keckinstr.bins = 1

  ;; Wavelength range
  keckinstr.wvmnx = [4000., 10000.]

  ;; Slit
  if keyword_set(SLIT) then keckinstr.swidth = slit $
  else keckinstr.swidth = 0.75  ; arcsec
  keckinstr.sheight = 20.0  ; arcsec

  return
end
