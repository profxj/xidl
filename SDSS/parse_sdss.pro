;+ 
; NAME:
; parse_sdss   
;   Version 1.0
;
; PURPOSE:
;    Fits a continuum to spectroscopic data interactively
;
; CALLING SEQUENCE:
;   
;   dla_sdssrich, fil
;
; INPUTS:
;
; RETURNS:
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
;   parse_sdss, fil
;
;
; PROCEDURES/FUNCTIONS CALLED:
; REVISION HISTORY:
;   05-Sep-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro parse_sdss, fil, flux, wave, conti, SIG=sig, IVAR=ivar, ZQSO=zqso, $
                   NPIX=npix, HEAD=head, CFIL=cfil, CDIR=cdir

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'parse_sdssqso, fil, flux [v1.0]'
    return
  endif 

; Optional Keywords

;  Open the file
  dat = xmrdfits(strtrim(fil,2), 0, head, /silent)
  flux = dat[*,0]
  npix = n_elements(flux)


; Wave
  iwave = sxpar(head, 'CRVAL1')
  wave= 10^(iwave + findgen(npix)*(0.0001))-10^(iwave)+10^(iwave)

; SIG
  sig = dat[*,2]

; IVAR
  if arg_present( IVAR ) then begin
      gds = where(sig GT 0.)
      ivar = fltarr(npix)
      ivar[gds] = 1. / sig[gds]^2
  endif
      
; Zqso
  if arg_present( ZQSO ) then zqso = sxpar(head, 'Z')

; Continuum
  if keyword_set( CDIR ) then begin
      sdss_name = strtrim( sxpar(head, 'NAME'), 2)+'-'+ $
        strtrim( sxpar(head, 'FIBERID'), 2)
      cfil = CDIR+sdss_name+'.fits'
  endif
  if keyword_set( CFIL ) then conti = xmrdfits(cfil, /silent)

  return
end

