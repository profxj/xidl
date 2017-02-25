; NAME:
; x_initdeimos
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
pro x_initdeimos, keckinstr, SLIT=slit, INFIL=infil, STR_TEL=str_tel, GRATING=grating

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'x_initdeimos, slit,  STR_TEL= [v1.0]'
      return
  endif 

  ;; Telescope
  x_initkeck, str_tel
  str_tel.name = 'KeckII'

  ;; Instrument
  keckinstr = {instrstruct}
  keckinstr.name = 'DEIMOS'
  keckinstr.MAG_PERP = 8.03  ; Modified to give observed resolution 0.75" maps to 4.5 pixels
  keckinstr.MAG_PARA = 8.03  ; Same
;  keckinstr.MAG_PERP = 5.735 ; Modified to give 0.1185" pixels
;  keckinstr.MAG_PARA = 5.735 ; Same
  keckinstr.PIXEL_SIZE = 15.0 ; in microns

  keckinstr.SCALE_PERP = str_tel.PLATE_SCALE* $
                         keckinstr.MAG_PERP*(keckinstr.PIXEL_SIZE/1000.) ; Arcsec
  keckinstr.SCALE_PARA = str_tel.PLATE_SCALE* $
                         keckinstr.MAG_PARA*(keckinstr.PIXEL_SIZE/1000.) ; Arcsec
          
  if not keyword_set(GRATING) then grating = '1200'
  keckinstr.grating = grating
  case grating of
     '600': keckinstr.R     = 11538.5    ; 1 pixel (native) dispersion
     '900': keckinstr.R     = 17307.8    ; 1 pixel (native) dispersion
     '1200': keckinstr.R    = 22727.3    ; 1 pixel (native) dispersion
     else: stop
  endcase

  ;; Detector
  keckinstr.readno = 2.6
  keckinstr.dark = 4.  ; electrons/pix/hr
  keckinstr.bind = 1
  keckinstr.bins = 1

  ;; Wavelength range
  keckinstr.wvmnx = [4000., 10000.]

  ;; Slit
  if keyword_set(SLIT) then keckinstr.swidth = slit $
  else keckinstr.swidth = 1.00  ; arcsec
  keckinstr.sheight = 10.0  ; arcsec

  return
end
