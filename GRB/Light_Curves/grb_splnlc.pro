;+ 
; NAME:
;  grb_splnlc
;   Version 1.1
;
; PURPOSE:
;    Given a grb structure (magnitude at a time t,  z, spectral slope)
;    calculates the luminsoity at any time 
;
; CALLING SEQUENCE:
;   
;   lum_nu = grb_calclum(grb, nu, [t])
;
; INPUTS:
;     grb -- GRB structure
;     nu  -- Frequency (GRB frame)
;     [t] -- Time (observer frame). Default = t0 in GRB structure
;
; RETURNS:
;   lum_nu= 
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
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Feb-2006 Written by JXP 
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;  Initializer ;;
pro grb_splnlc, grb

common grb_splinelc

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'grb_splnlc, grb [v1.0]'
    return
  endif 

  ;; Open the file
  case round(grb.mag_wv) of
      6588: filter = 'Rc'
      else: stop
  endcase
  fil = getenv('XIDL_DIR')+'/GRB/Light_Curves/grb'+ $
        strtrim(grb.nam,2) + '_' + filter + '.dat'

  readcol, fil, src, tdays, filt, mag, sig_mag, format='A,D,A,F,F'
  gd = where(strtrim(filt,2) EQ filter)
  tdays = tdays[gd]
  mag = mag[gd]
  sig_mag = sig_mag[gd]

  ;; Sort and save
  srt = sort(tdays)
  grbsp_t = tdays[srt] * 24 * 3600.  ; seconds
  grbsp_mag = mag[srt]

  ;; Spline
  grbsp_splin = spl_init(grbsp_t, grbsp_mag, /double)

  return
end
