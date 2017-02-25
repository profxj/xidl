; NAME:
; x_initkeck
;    Version 1.0
;
; PURPOSE:
;    Initialize an instrument structure for APFSPEC
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
pro apfspec_setdecker, apfinstr, decker

  case strtrim(DECKER,2) of 
     'N': begin
          apfinstr.swidth = 0.5
          apfinstr.sheight = 8.0
     end
     'S': begin
          apfinstr.swidth = 0.75
          apfinstr.sheight = 8.0
     end
     'M': begin
          apfinstr.swidth = 1.0
          apfinstr.sheight = 8.0
     end
     'W': begin
          apfinstr.swidth = 1.0
          apfinstr.sheight = 3.0
     end
     'O': begin
          apfinstr.swidth = 8.0
          apfinstr.sheight = 8.0
     end
     'T': begin
          apfinstr.swidth = 2.0
          apfinstr.sheight = 3.0
     end
     'B': begin
          apfinstr.swidth = 2.0
          apfinstr.sheight = 8.0
     end
     else: stop
  endcase
  return
end


pro x_initapfspec, apfinstr, flg, INFIL=infil, STR_TEL=str_tel, DECKER=decker

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'x_initapfspec, apfinstr, [flg], KECKTEL= [v1.0]'
      return
  endif 

  if not keyword_set(FLG) then flg = 1
  if not keyword_set(STR_TEL) then begin
     x_initapftel, str_tel
     str_tel.name = 'APF'
   endif

  apfinstr = {instrstruct}
  case flg of 
      1: begin ;;  APFSPEC
          apfinstr.MAG_PERP = 5.5136 ; Magnification perpendicular to dispesion
          apfinstr.MAG_PARA = 4.9913 ; Magnification parallel to dispesion
          apfinstr.PIXEL_SIZE = 13.5 ; in microns

          apfinstr.SCALE_PERP = str_tel.PLATE_SCALE* $
            apfinstr.MAG_PERP*(apfinstr.PIXEL_SIZE/1000.)
          apfinstr.SCALE_PARA = str_tel.PLATE_SCALE* $
            apfinstr.MAG_PARA*(apfinstr.PIXEL_SIZE/1000.)
          apfinstr.R     = 282160.8 
          
          apfinstr.MLAMBDA = 465980.24  ; measured from order central wavelengths

          apfinstr.DELY    = 0.154693  ; = f2 * Ac 
          ;;
          apfinstr.readno = 3.75
          apfinstr.dark = 7.3
          apfinstr.bind = 1
          apfinstr.bins = 2
          ;; Wavelength range
          apfinstr.wvmnx = [3742., 7700.]
      end
      else: stop
  endcase

  ;; INFIL?
  if keyword_set(INFIL) then begin
      readcol, infil, card, val, FORMAT='A,A' 
      ;; SWIDTH
      mtch = where(strtrim(card,2) EQ 'SWIDTH',nmt)
      if nmt NE 0 then apfinstr.swidth = float(val[mtch[0]])
      ;; Decker
      mtch = where(strtrim(card,2) EQ 'DECKER',nmt)
      if nmt NE 0 then decker = val[mtch[0]]
      ;; BINC
      mtch = where(strtrim(card,2) EQ 'BINC',nmt)
      if nmt NE 0 then apfinstr.bind = long(val[mtch[0]])
      ;; BINR
      mtch = where(strtrim(card,2) EQ 'BINR',nmt)
      if nmt NE 0 then apfinstr.bins = long(val[mtch[0]])
  endif

  ;; Slit
  if keyword_set(DECKER) then begin
      apfspec_setdecker, apfinstr, decker
  endif else begin ;; Defaults to C5
      decker = 'big'
      apfinstr.swidth = 1.0
      apfinstr.sheight = 3.0
  endelse

  return
end
