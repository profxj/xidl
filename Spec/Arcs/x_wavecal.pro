;+ 
; NAME:
; x_wavecal   
;   Version 1.0
;
; PURPOSE:
;    Wavelength calibrate a given input img/spectrum
;
; CALLING SEQUENCE:
;   
;   calib = x_wavecal(img/spec, [extrct], LINELIST=)
;
; INPUTS:
;   img/spec   - 2D image or 1D spectrum
;   [extrct]   - Extraction structure
;
; RETURNS:
;   FSTRCT - Structure defining the fit
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   LINELIST   - List of lines
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_wavecal, spec
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   07-Dec-2001 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_wavecal, arc, extrct, LINELIST=linelist, DEBUG=debug, $
                    REDBLUE=redblue, DISP=disp, CALIB=calib, FLUX=flux


;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'calib = x_wavecal(arc, [extrct], LINELIST=,/DEBUG, /REDBLUE'
    print, '        DISP=, CALIB=, FLUX=) [v1.0]'
    return, -1
  endif 

;  Optional Keywords

;  Read in image if necessary

  dat = x_readimg(arc, /fscale)

;  Extract 1D if necessary 
  
  sz = size(dat)
  if sz[0] EQ 2 then begin
      if not keyword_set( extrct ) then begin
          print, 'Must provide extraction structure'
          return, -1
      endif
      ; EXTRACT
      spec = x_apall(dat, CLINE=extrct.cline, APER=extrct.aper, $
                     TRACE=*extrct.trace, /NOOV, /NOSKY)
  endif else spec = dat
  delvarx, dat

; Begin identify

  x_identify, spec, calib, LINELIST=linelist, DEBUG=debug, WAVE=wave, $
    REDBLUE=redblue, DISP=disp, FLUX=flux

  return, wave
end
