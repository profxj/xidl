; NAME:
; x_initkeck
;    Version 1.0
;
; PURPOSE:
;    Initialize an instrument structure for HIRES
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
pro hires_setdecker, keckinstr, decker

  case strtrim(DECKER,2) of 
      'B2': begin
          keckinstr.swidth = 0.574
          keckinstr.sheight = 7.0
      end
      'C1': begin
          keckinstr.swidth = 0.861
          keckinstr.sheight = 7.0
      end
      'C5': begin
          keckinstr.swidth = 1.148
          keckinstr.sheight = 7.0
      end
      'B5': begin
          keckinstr.swidth = 0.861
          keckinstr.sheight = 7.0
      end
      'D1': begin
          keckinstr.swidth = 1.148
          keckinstr.sheight = 14.0
      end
      'D3': begin
          keckinstr.swidth = 1.72
          keckinstr.sheight = 7.0
      end
      'E4': begin
          keckinstr.swidth = 0.400
          keckinstr.sheight = 7.0
      end
      'E5': begin
          keckinstr.swidth = 0.800
          keckinstr.sheight = 1.0
      end
      else: stop
  endcase
  return
end


pro x_inithires, keckinstr, flg, INFIL=infil, STR_TEL=str_tel, DECKER=decker

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'x_initkeck, keckinstr, [flg], KECKTEL= [v1.0]'
      return
  endif 

  if not keyword_set(FLG) then flg = 1
  if not keyword_set(STR_TEL) then begin
      case flg of 
          1: x_initkeck, str_tel
          2: x_initkeck, str_tel
          3: x_inittmt, str_tel
          else: stop
      endcase
      str_tel.name = 'KeckI'
   endif

  keckinstr = {instrstruct}
  case flg of 
      1: begin ;; Old HIRES
          keckinstr.MAG_PERP = 5.776 ; Magnification perpendicular to dispesion
          keckinstr.MAG_PARA = 8.6713 ; Magnification parallel to dispesion
          keckinstr.PIXEL_SIZE = 24.0 ; in microns

          keckinstr.SCALE_PERP = str_tel.PLATE_SCALE* $
            keckinstr.MAG_PERP*(keckinstr.PIXEL_SIZE/1000.)
          keckinstr.SCALE_PARA = str_tel.PLATE_SCALE* $
            keckinstr.MAG_PARA*(keckinstr.PIXEL_SIZE/1000.)
          keckinstr.R     = 135000      
          
          keckinstr.MLAMBDA = 356385.3016  ; = 2*sigma*sin(delta)*cos(theta)
                                        ;   52.676 grooves/mm 
                                        ;   70.43  degress blaze 
          keckinstr.DELY    = 0.154693  ; = f2 * Ac 
          ;;
          keckinstr.readno = 4.3
          keckinstr.dark = 2.0
          keckinstr.bind = 1
          keckinstr.bins = 2
          ;; Wavelength range
          keckinstr.wvmnx = [3000., 9000.]
      end
      2:  begin
          keckinstr.MAG_PERP = 5.776 ; Magnification perpendicular to dispesion
          keckinstr.MAG_PARA = 8.6713 ; Magnification parallel to dispesion
          keckinstr.PIXEL_SIZE = 15.0 ; in microns

          keckinstr.SCALE_PERP = str_tel.PLATE_SCALE* $
            keckinstr.MAG_PERP*(keckinstr.PIXEL_SIZE/1000.)
          keckinstr.SCALE_PARA = str_tel.PLATE_SCALE* $
            keckinstr.MAG_PARA*(keckinstr.PIXEL_SIZE/1000.)
          keckinstr.R     = 135000  * 24. / keckinstr.pixel_size
          
          keckinstr.MLAMBDA = 356385.3016  ; = 2*sigma*sin(delta)*cos(theta)
                                        ;   52.676 grooves/mm 
                                        ;   70.43  degress blaze 
          keckinstr.DELY    = 0.154693  ; = f2 * Ac 
          ;;
          keckinstr.readno = 2.2
          keckinstr.dark = 1.0
          keckinstr.bind = 1
          keckinstr.bins = 2
          ;; Wavelength range
          keckinstr.wvmnx = [3000., 9000.]
      end
      3: begin ;; MTHR
          keckinstr.MAG_PERP = 9.14 ; Magnification perpendicular to dispersion
          keckinstr.MAG_PARA = 9.46 ; Magnification parallel to dispersion
          keckinstr.PIXEL_SIZE = 15.0 ; in microns

          keckinstr.SCALE_PERP = str_tel.PLATE_SCALE* $
            keckinstr.MAG_PERP*(keckinstr.PIXEL_SIZE/1000.)
          keckinstr.SCALE_PARA = str_tel.PLATE_SCALE* $
            keckinstr.MAG_PARA*(keckinstr.PIXEL_SIZE/1000.)
          keckinstr.R     = 709379.   ; 0.00902 A at 6398.6
          
          keckinstr.MLAMBDA = 614073.6  ; = 2*sigma*sin(delta)*cos(theta)
                                        ;   52.676 grooves/mm 
                                        ;   70.43  degress blaze 
;          keckinstr.DELY    = 0.154693  ; Old value,  not using it.  Probably the cross-disperser arcsec/Ang.
          ;;
          keckinstr.readno = 2.2
          keckinstr.dark = 1.0
          keckinstr.bind = 4
          keckinstr.bins = 4
          ;; Wavelength range
          keckinstr.wvmnx = [3000., 9000.]
      end
      else: stop
  endcase

  ;; INFIL?
  if keyword_set(INFIL) then begin
      readcol, infil, card, val, FORMAT='A,A' 
      ;; SWIDTH
      mtch = where(strtrim(card,2) EQ 'SWIDTH',nmt)
      if nmt NE 0 then keckinstr.swidth = float(val[mtch[0]])
      ;; Decker
      mtch = where(strtrim(card,2) EQ 'DECKER',nmt)
      if nmt NE 0 then decker = val[mtch[0]]
      ;; BINC
      mtch = where(strtrim(card,2) EQ 'BINC',nmt)
      if nmt NE 0 then keckinstr.bind = long(val[mtch[0]])
      ;; BINR
      mtch = where(strtrim(card,2) EQ 'BINR',nmt)
      if nmt NE 0 then keckinstr.bins = long(val[mtch[0]])
  endif

  ;; Slit
  if keyword_set(DECKER) then begin
      hires_setdecker, keckinstr, decker
  endif else begin ;; Defaults to C5
      decker = 'C5'
      keckinstr.swidth = 1.1
      keckinstr.sheight = 7.0
  endelse

  return
end
