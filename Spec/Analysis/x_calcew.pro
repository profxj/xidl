;+ 
; NAME:
; x_calcew   
;   Version 1.0
;
; PURPOSE:
;    Fits a continuum to spectroscopic data interactively
;
; CALLING SEQUENCE:
;   
;   ew = x_calcew(wav, fx, lmts, [sig], sigew, WAV= )
;
; INPUTS:
;   wav      - Wavelength 
;   fx       - Flux (presumed normalized)
;   lmts     - Limits of integration :: Assumed pix unless wav set
;   [sig]    - Error in fx
;
; RETURNS:
;   ew   - Equivalent width
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;   sigew  - Error in ew
;
; COMMENTS:
;
; EXAMPLES:
;   ew = x_calcew( fx, sig, [50, 100], ERR=sigew) 
;
;
; PROCEDURES/FUNCTIONS CALLED:
; REVISION HISTORY:
;   05-Mar-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_calcew, wav, fx, lmts, sig, sigew, REDUCED=reduced, FPIX=fpix

;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'ew = x_calcew(fx, lmts, [sig, sigew], /FPIX, /REDUCED) [v1.0]'
    return, -1
  endif 

; Optional Keywords

; Pixel

  if keyword_set( FPIX ) then pix = lindgen(lmts[1]-lmts[0]+1)+lmts[0] $
  else begin
      pix = where(wav GE lmts[0] AND wav LE lmts[1], npix)
      if npix EQ 0 then begin
          print, 'x_calcew: No pixels in the defined limits!'
          return, -1
      endif
  endelse

  ; Shift
  pixp = pix+1L
  ; Add
  ew = total( (1.-fx[pix])*(wav[pixp]-wav[pix]) )
  ; Error
  if keyword_set( SIG ) then $
    sigew = sqrt(total( (sig[pix]*(wav[pixp]-wav[pix]))^2 ))

  ; Reduced
  if keyword_set( REDUCED ) then begin
      ew = ew/ wav[pix[n_elements(pix)/2]]
      if keyword_set( SIG ) then sigew = sigew/wav[pix[n_elements(pix)/2]]
  endif

  return, ew
end
      
