;+ 
; NAME:
; x_wavecal   
;   Version 1.1
;
; PURPOSE:
;    Wavelength calibrate a given input Arc img (or spectrum).  This
;    routine simply extracts the Arc spectrum (if necessary) and then
;    calls x_identify
;
; CALLING SEQUENCE:
;  wave = x_wavecal( arc, [extrct], LINELIST=, /DEBUG, $
;                   /REDBLUE, DISP=, CALIB=, FLUX=)
;
; INPUTS:
;   img/spec   - 2D image or 1D spectrum
;   [extrct]   - Extraction structure (required for an input 2D image)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   LINELIST   - Reference list (user can set interactively)
;   DISP       - Guess at dispersion (A or km/s per pix [pos/neg
;                value])
;   /FLUX       - Assumes wide slit and therefore wide arc lines
;   /REDBLUE    - Spectrum runs from red to blue (not blue to red)
;   /DEBUG
;   /ROT       -  Transpose the Arc image (code requires orders
;                parallel to rows)
;
; OPTIONAL OUTPUTS:
;   CALIB= - Structure defining the fit
;   SPEC=  - 1D extracted spectrum
;
; COMMENTS:
;
; EXAMPLES:
;   wave = x_wavecal(spec)
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_readimg
;  x_identify
;
; REVISION HISTORY:
;   07-Dec-2001 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_wavecal, arc, extrct, LINELIST=linelist, DEBUG=debug, $
                    REDBLUE=redblue, DISP=disp, CALIB=calib, FLUX=flux, $
                    ROT=rot, SPEC=spec


;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'calib = x_wavecal(arc, [extrct], LINELIST=,/DEBUG, /REDBLUE'
    print, '        DISP=, CALIB=, FLUX=) [v1.1]'
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
                     TRACE=*extrct.trace, /NOOV, /NOSKY, ROT=rot)
  endif else spec = dat
  delvarx, dat

; Begin identify

  x_identify, spec, calib, LINELIST=linelist, DEBUG=debug, WAVE=wave, $
    REDBLUE=redblue, DISP=disp, FLUX=flux

  return, wave
end
