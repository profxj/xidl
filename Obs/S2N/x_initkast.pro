; NAME:
; x_initkast
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
pro x_initkast, kastinstr, SLIT=slit, INFIL=infil, STR_TEL=str_tel, $
                GRISM=grism, GRATING=grating

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'x_initkast, slit,  STR_TEL= [v1.0]'
      return
  endif 

  ;; Telescope
  x_initlick, str_tel

  ;; Instrument
  kastinstr = replicate({instrstruct}, 2) ;; 2 Sides
  kastinstr.MAG_PERP = 20.9  ; Modified to give 0.15" pixels
  kastinstr.MAG_PARA = 20.9  ; Same

  kastinstr[0].PIXEL_SIZE = 15.0 ; in microns (0.43")
  kastinstr[1].PIXEL_SIZE = 27.0 ; in microns (0.78")

  kastinstr.SCALE_PERP = str_tel.PLATE_SCALE* $
                         kastinstr.MAG_PERP*(kastinstr.PIXEL_SIZE/1000.) ; Arcsec
  kastinstr.SCALE_PARA = str_tel.PLATE_SCALE* $
                         kastinstr.MAG_PARA*(kastinstr.PIXEL_SIZE/1000.) ; Arcsec

  ;; Resolution
  if not keyword_set(GRISM) then grism = 'G1'
  kastinstr[0].grating = grism
  case grism of
     'G1': kastinstr[0].R     = 2344    ; 1 pixel (native) dispersion
     'G2': kastinstr[0].R     = 4254    ; 1 pixel (native) dispersion
     'G3': kastinstr[0].R     = 5492    ; 1 pixel (native) dispersion
     else: stop
  endcase

  if not keyword_set(GRATING) then grating = '600/7500'
  kastinstr[1].grating = grating
  case grating of
     '600/7500': kastinstr[1].R     = 3164    ; 1 pixel (native) dispersion
     else: stop
  endcase
          
  ;; Blue Detector
  kastinstr.readno = 3.7
  kastinstr.dark = 0.001

  ;; Red Detector
  kastinstr.readno = 12.5
  kastinstr.dark = 0.001
  
  ;;
  kastinstr.bind = 1
  kastinstr.bins = 1

  ;; Wavelength range
  kastinstr[0].wvmnx = [3000., 6000.]
  kastinstr[1].wvmnx = [5000., 10000.]

  ;; Slit
  if keyword_set(SLIT) then kastinstr.swidth = slit $
  else kastinstr.swidth = 0.75  ; arcsec
  kastinstr.sheight = 120.0  ; arcsec

  return
end
