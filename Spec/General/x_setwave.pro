;+ 
; NAME:
;  x_setwave   
;   Version 1.1
;
; PURPOSE:
;    Sets a wavelength array given a header from a FITS file
;
; CALLING SEQUENCE:
;   
;   wave = x_setwave(head, ntot)
;
; INPUTS:
;   head     - Header
;   ntot     - number of pixels
;
; RETURNS:
;   wave     - Wavelength array
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
;   wave = x_setwave(head, 1000L)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Nov-2001 Written by JXP 
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_setwave, head, ntot, CRVAL1=crval1, CDELT1=cdelt1, CRPIX1=crpix1, SILENT = SILENT

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'wave = x_setwave(head, ntot) [v1.1]'
    return, -1
  endif 

;  Optional Keywords

;  Parse the header
  
  crpix1 = sxpar(head, 'CRPIX1', count=count, silent = silent)
  if count EQ 0 then crpix1 = 1.0d
  crval1 = sxpar(head, 'CRVAL1', count=count, silent = silent)
  if count EQ 0 then crval1 = 1.0d
  cdelt1 = sxpar(head, 'CDELT1', count=count, silent = silent)
  if count EQ 0 then cdelt1 = 1.0d
;
  cd1_1 = sxpar(head, 'CDELT1', count=count, silent = silent)
  if count EQ 0 or abs(cd1_1-1.) LT 1e-5 then begin
      cdelt1 = sxpar(head, 'CD1_1', count=count, silent = silent)
      if count EQ 0 then cdelt1 = 1.0d
  endif

; ctype
  ctype1 = sxpar(head, 'CTYPE1', count=count, silent = silent)
  dcflag = sxpar(head, 'DC-FLAG', count = dc_count, silent = silent)
  case strmid(strtrim(ctype1, 2), 0, 4) of 
      'POLY': begin
          print, 'Not able to handle poly_lambda'
          return, -1
      end
      'LINE': begin
          if dcflag EQ 1 then stype = 3 else stype = 2
          ctype2 = sxpar(head, 'CTYPE2', count = count, silent = silent)
          if strmid(strtrim(ctype2, 2), 0, 4) EQ 'PIXE' then begin
             flg_lp  = 1
             stype = 3
          endif 
       end
      ;; Added by JFH 09/11 to handle Murphy's UVES spectra
      'LOGL': begin
         if dcflag EQ 1 then stype = 3 else stype = 2
      end
      'LAMB': begin
          if dcflag EQ 1 then stype = 3 else stype = 2
      end
      'WAVE': begin  ; UVES
          if dcflag EQ 1 then stype = 3 else stype = 2
       end
      'PIXE': begin ; X-shooter?
         linear = sxpar(head, 'CTYPE2', count = count, silent = silent)
         IF strmatch(linear, '*LINEAR*') THEN stype = 2
      end
      '0': BEGIN
         IF dc_count EQ 0 THEN BEGIN
            IF crval1 LE 0 THEN stype = 3 ELSE stype = 2
         ENDIF ELSE BEGIN
            if dcflag EQ 1 then stype = 3 else stype = 2
         ENDELSE
      END
      else: BEGIN 
         if dcflag EQ 1 then stype = 3 else stype = 2
                                ;splog, 'Warning: No CTYPE1 in header'
         ;;stype = -2
      end
   endcase
; Make wave array
;  if strmatch(sxpar(head,'INSTRUME'), 'XSHOOTER') and (not keyword_set(flg_lp)) $
;  then SCL = 10. else SCL = 1.
  SCL=1.
  
  case stype of
      -2: return, temporary(dindgen(ntot))
      2: return, SCL* (crval1 + ( cdelt1 * (temporary(dindgen(ntot)) $
                                      + 1.0d -crpix1) ) ); Linear
      3: return, SCL * 10.d^(crval1 + ( cdelt1 * (temporary(dindgen(ntot)) $
                                          + 1.0d - crpix1) )) ; Log
      else: return, -1
  endcase
end
