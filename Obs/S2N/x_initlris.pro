; NAME:
; x_initlris
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
pro x_initlris, lrisinstr, SLIT=slit, INFIL=infil, STR_TEL=str_tel, $
                GRISM=grism, GRATING=grating

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'x_initlris, slit,  STR_TEL= [v1.0]'
      return
  endif 

  ;; Telescope
  x_initkeck, str_tel
  str_tel.name = 'KeckI'

  ;; Instrument
  lrisinstr = replicate({instrstruct}, 2) ;; 2 Sides
  lrisinstr.name = 'LRIS'
  lrisinstr.MAG_PERP = 6.5  ; Modified to give 0.135" pixels
  lrisinstr.MAG_PARA = 6.5  ; Same

  lrisinstr[0].PIXEL_SIZE = 15.0 ; in microns 
  lrisinstr[1].PIXEL_SIZE = 15.0 ; in microns 

  lrisinstr.SCALE_PERP = str_tel.PLATE_SCALE* $
                         lrisinstr.MAG_PERP*(lrisinstr.PIXEL_SIZE/1000.) ; Arcsec
  lrisinstr.SCALE_PARA = str_tel.PLATE_SCALE* $
                         lrisinstr.MAG_PARA*(lrisinstr.PIXEL_SIZE/1000.) ; Arcsec

  ;; Resolution
  if not keyword_set(GRISM) then grism = 'B600'
  lrisinstr[0].grating = grism
  case grism of
     'B600': lrisinstr[0].R     = 7500    ; 1 pixel (native) dispersion
     'B300': lrisinstr[0].R     = 3304    ; 1 pixel (native) dispersion
     else: stop
  endcase

  if not keyword_set(GRATING) then grating = '600/7500'
  lrisinstr[1].grating = grating
  case grating of
     '600/7500': lrisinstr[1].R     = 11820    ; 1 pixel (native) dispersion
     '600/10000': lrisinstr[1].R     = 11820    ; 1 pixel (native) dispersion
     '400/8500': lrisinstr[1].R     = 8151    ; 1 pixel (native) dispersion
     '831/8200': lrisinstr[1].R     = 16303    ; 1 pixel (native) dispersion
     '1200/9000': lrisinstr[1].R     = 23640    ; 1 pixel (native) dispersion
     else: stop
  endcase
          
  ;; Blue Detector
  lrisinstr.readno = 3.7
  lrisinstr.dark = 0.001

  ;; Red Detector
  lrisinstr.readno = 4.5
  lrisinstr.dark = 0.001
  
  ;;
  lrisinstr.bind = 1
  lrisinstr.bins = 1

  ;; Wavelength range
  lrisinstr[0].wvmnx = [3000., 6000.]
  lrisinstr[1].wvmnx = [5000., 10000.]

  ;; Slit
  if keyword_set(SLIT) then lrisinstr.swidth = slit $
  else lrisinstr.swidth = 1.00  ; arcsec
  lrisinstr.sheight = 120.0  ; arcsec

  return
end
