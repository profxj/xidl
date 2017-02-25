;+ 
; NAME:
; x_calcew   
;   Version 1.1
;
; PURPOSE:
;    Calculates the EW given a wavelength and flux array and limits to
;    integrate over.  Can also calculate the error in the EW as well
;    as the reduced values.
;
; CALLING SEQUENCE:
;   ew = x_calcew(wav, fx, lmts, [sig, sigew], /REDUCED, FPIX=)
;
; INPUTS:
;   wav      - Rest wavelength [Angstromgs]
;   fx       - Flux (presumed normalized)
;   lmts     - Limits of integration (wavelength unless /FPIX is set)
;   [sig]    - Error in fx
;
; RETURNS:
;   ew   - Equivalent width
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
; /REDUCED -- Divide EW value by lambda
; /FPIX  -- If flagged, this indicates the lmts array contains the
;           pixel endpoints for the EW calculation.
;
; OPTIONAL OUTPUTS:
;   sigew  - Error in ew
;
; COMMENTS:
;  The flux array should be normalized to unity 
;
; EXAMPLES:
;   ew = x_calcew( fx, sig, [50, 100], ERR=sigew, /FPIX) 
;
; PROCEDURES/FUNCTIONS CALLED:
;
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


  ;; Delta wave
  pixp = pix+1L
  dwv = abs(wav[pixp]-wav[pix])

  ;; Cut out bad pixels
  gd = where(sig[pix] GT 0., ngd)
  if ngd EQ 0 then return, -1.
  gdpix = pix[gd]

  ;; Add
  ew = total( (1.-fx[gdpix])*dwv[gdpix]) 
  ;; Error
  if keyword_set( SIG ) then $
    sigew = sqrt(total( (sig[gdpix]*dwv[gdpix])^2 ))
    ;sigew = sqrt(total( (sig[pix]*(wav[pixp]-wav[pix]))^2 ))

  ;; Reduced
  if keyword_set( REDUCED ) then begin
      ew = ew/ wav[gdpix[n_elements(gdpix)/2]]
      if keyword_set( SIG ) then sigew = sigew/wav[gdpix[n_elements(gdpix)/2]]
  endif

  return, ew
end
      
