;+ 
; NAME:
; x_fitswave
;   V1.0
;
; PURPOSE:
;    Returns the wavelength array after passed a header
;
; CALLING SEQUENCE:
;   
;   wave = x_fitswave(head)
;
; INPUTS:
;   head    - Header
;
; RETURNS:
;   wave  - Double array of wavelength values
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wave = x_fitswave(head)
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   13-Sep-2001 Written by JXP
;-
;------------------------------------------------------------------------------

function x_fitswave, head

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'wave = x_fitswave(head) [v1.0]'
    return, -1
  endif 

; Keywords

; Make Wavelength array

  ctype = strtrim(sxpar(head, 'CTYPE1'),2)
  dcflag = sxpar(head, 'DC-FLAG')
  npix = sxpar(head, 'NAXIS1')
  w0 = sxpar(head, 'CRVAL1')
  cdelt = sxpar(head, 'CDELT1')
  crpix = sxpar(head, 'CRPIX1')

  case ctype of
      'LINEAR' : begin
          case dcflag of 
              0 : begin
                  wave = dindgen(npix)+1
                  wave = w0 + cdelt*(wave-double(crpix))
              end
              1 : begin
                  wave = dindgen(npix)+1   ; CRPIX1=1 for FORTRAN arrays
                  wave = 10.0D^(w0 + cdelt*(wave-double(crpix)))
              end
              else : return, -1
          endcase
      end
      else : return, -1
  endcase

  return, wave

end

