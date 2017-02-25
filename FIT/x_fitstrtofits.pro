;+ 
; NAME:
; x_fitstrtofits
;   Version 1.1
;
; PURPOSE:
;    Writes a fitstr as a binary fits file (or read in a fit struct
;    from a binary FITS file)
;
; CALLING SEQUENCE:
;   x_fitstrtofits, fit_str, fits_fil, /REVERSE
;
; INPUTS:
;    fit_str   - Fit structure
;    fits_fil  - Fitsfile
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   /REVERSE  -- Default is to write the file.  If /REVERSE is set
;                then the FITS file is read into a FIT structure
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_fitstrtofits, fit_str, fits_fil
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   15-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_fitstrtofits, fit_str, fits_fil, REVERSE=reverse

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'x_fitstrtofits, fit_str, fits_fil, /REVERSE [v1.1]'
    return
  endif 

; Create an anonymous structure

  if not keyword_set( REVERSE ) then begin
      dumstr = { $
                 func: fit_str.func, $
                 nord: fit_str.nord, $
                 nrm: fit_str.nrm, $
                 lsig: fit_str.lsig, $
                 hsig: fit_str.hsig, $
                 niter: fit_str.niter, $
                 minpt: fit_str.minpt, $
                 maxrej: fit_str.maxrej, $
                 flg_rej: fit_str.flg_rej, $
                 rms: fit_str.rms $
               }
      case strtrim(dumstr.func,2) of
          'BSPLIN': dumstr2 = { fullbkpt: (*fit_str.ffit).fullbkpt,$
                                bkmask: (*fit_str.ffit).bkmask,$
                                nord: (*fit_str.ffit).nord,$
                                coeff: (*fit_str.ffit).coeff,$
                                icoeff: (*fit_str.ffit).icoeff $
                              }
          else: dumstr2 = { ffit: *fit_str.ffit }
      endcase

      mwrfits, dumstr, fits_fil, /create
      mwrfits, dumstr2, fits_fil
  endif else begin  ; READ FROM FITS FILE
      dumstr = mrdfits(fits_fil, 1, /silent)
      fit_str = { fitstrct }
      copy_struct, dumstr, fit_str, EXCEPT_TAGS=['ffit']
      if tag_exist(dumstr, 'ffit') EQ 1 then $
        fit_str.ffit = ptr_new(dumstr.ffit) $
      else begin
          dumstr2 = mrdfits(fits_fil, 2, /silent)
          case strtrim(fit_str.func,2) of 
              'BSPLIN': fit_str.ffit = ptr_new(dumstr2)
              else: fit_str.ffit = ptr_new(dumstr2.ffit)
          endcase
      endelse
  endelse

  return
end
