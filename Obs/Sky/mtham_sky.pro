;+ 
; NAME:
; mtham_sky
;    Version 1.1
;
; PURPOSE:
;    Gives an estimate of the sky brightness for Mauna Kea
;
; CALLING SEQUENCE:
;  mag = maunakea_sky( wave, phase )
;
; INPUTS:
;  wave=  -- Wavelength array
;  phase= -- Phase of the moon [days]
;
; RETURNS:
;  mag=  -- Sky brightness in AB mags
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
;   23-Aug-2011 Written by JXP based on HIRES S2N code
;-
;------------------------------------------------------------------------------
function mtham_sky, wave, iphase

common mtham_sky_cmmn, mtham_swv, mtham_sfx, mtham_sfil
;     Sky brightness at Mt Hamilton
;..............................................................................
;

  ;; Using MaunaKea for now
;  msky = maunakea_sky(wave, iphase)

  ;; Using an empirical 'dark' sky measurement
  tfil = getenv('XIDL_DIR')+'/Obs/Sky/Empirical/lick_sky_d55_2011aug29.fits'
  if not keyword_set(MTHAM_SFIL) then mtham_sfil = ''
  if not strmatch(mtham_sfil,tfil) then begin
     mtham_swv = xmrdfits(tfil, 0)
     mtham_sfx = xmrdfits(tfil, 1)
     mtham_sfil = tfil
  endif

  ;; Interpolate
  flam = interpol(mtham_sfx, mtham_swv, wave)
  a = where(wave LT mtham_swv, na)
  if na GT 0 then flam[a] = mtham_sfx[0]

  fnu = flam / (3e10) * wave * (wave * 1e-8)
  msky = -1. * alog10(fnu) / 0.4 - 48.6
;        
  return, msky
end
