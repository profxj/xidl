;+ 
; NAME:
; x_setline
;
; PURPOSE:
;    Given an atomic wavelength, return the fvalue and name
;
; CALLING SEQUENCE:
;   
;   x_mkabslin
;
; INPUTS:
;   wave       - ionic transition
;
; RETURNS:
;   f          - oscillator strength
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;   nam        - Name of transition
;
; COMMENTS:
;
; EXAMPLES:
;   x_setline, 1215.6701, fval, name
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   17-Oct-2002 Written by JXP
;-
;------------------------------------------------------------------------------

function x_setline, wave, LINFIL=linfil

 ;
  if  N_params() NE 1  then begin 
    print,'Syntax - ' + $
             'line = x_setline(wave, LINFIL) [V1.1]'
    return, -1
  endif 

; Optional keywords

  if not keyword_set( LINFIL ) then $
    linfil = getenv('XIDL_DIR')+'/Spec/Lines/all_lin.fits'

; Grab the data
  all_lin = xmrdfits( linfil, 1, structyp='abslinstrct', /silent)

; Search for a match
  
  a = where(abs(all_lin.wrest-wave) LT 0.01, na)
  case na of 
      0: begin
          print, 'x_setline:  No match in database!'
          return, -1
      end
      1: return, all_lin[a]
      else: begin
          print, 'x_setline:  Multiple hits.  Returning: ', all_lin[a[0]].ion
          return, all_lin[a[0]]
      end
  endcase
          
  return, -1
end
