;+ 
; NAME:
; sdss_getqmag   
;    Version 1.1
;
; PURPOSE:
;  Return the magnitude of a QSO with plate+fiber
;
; CALLING SEQUENCE:
; mag = sdss_getqmag(qalstr, platfib, /ALL)
;
; INPUTS:
;
; RETURNS:
;  mag -- Array of u-z magnitudes for SDSS QSO
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
; QSOFIL=  -- File containing the QSO summary FITS table [default:
;             dr5_qso.fits]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  mag = sdss_getqmag([424,124])
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   13-Oct-2004 Written by JXP
;-
;------------------------------------------------------------------------------

function sdss_getqmag, platfib, QSOSTR=qsostr, QSOFIL=qsofil

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'indx = sdss_fndqso, qalstr, platfib, /ALL [v1.1]'
      return, -1
  endif 
  
  if not keyword_set(QSOFIL) then qsofil = $
    getenv('SDSSPATH')+'/DR5_QSO/dr5_qso.fits'

  if not keyword_set(QSOSTR) then qsostr = xmrdfits(qsofil, 1, /silent)

  ;; 

  sz = size(platfib, /dimensions)

  if n_elements(sz) EQ 1 then begin
      mtch = where(qsostr.plate EQ platfib[0] AND $
                   qsostr.fiberid EQ platfib[1], nm)
      if nm EQ 0 then begin
          print, 'sdss_getqmag: No match for ', platfib
          return, replicate(-1,5L)
      endif else return, [ qsostr[mtch].psf_u, $
                           qsostr[mtch].psf_g, $
                           qsostr[mtch].psf_r, $
                           qsostr[mtch].psf_i, $
                           qsostr[mtch].psf_z]
  endif else begin
      finr = fltarr(5,sz[1])
      for qq=0L,sz[1]-1 do $
        finr[*,qq] = sdss_getqmag(platfib[*,qq], QSOSTR=qsostr)
      return, finr
  endelse
                   
  return, -1
end
