;+ 
; NAME:
; dla_esitosdss   
;   Version 1.0
;
; PURPOSE:
;    Fits a continuum to spectroscopic data interactively
;
; CALLING SEQUENCE:
;   
;   dla_esitosdss, fil, NHI, Z, NPIX=, SNR=, OUTFIL=, NSIG=
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
;   dla_esitosdss, fil, NHI, Z
;
;
; PROCEDURES/FUNCTIONS CALLED:
; REVISION HISTORY:
;   09-Dec-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro dla_esitosdss, fil, wave, fx, sig

;
  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
             'dla_esitosdss, fil, wave, fx, sig [v1.0]'
    return
  endif 

; Optional Keywords

; Read data file
  sig_fil = strmid(fil,0,strlen(fil)-6)+'e.fits'
  esi_fx = x_readspec(fil, FIL_SIG=sig_fil, WAV=esi_wv, sig=esi_sig)

; Smooth

  kernel = gauss_kernel(7.)
  smooth_fx = convol(esi_fx, kernel)
  smooth_sig = convol(esi_sig, kernel)

; Rebin
  wave = rebin(esi_wv[0:27299], 3900L)
  fx = rebin(smooth_fx[0:27299], 3900L)
  sig = rebin(smooth_sig[0:27299], 3900L)

; Check EW
  a = where(esi_wv GT 7900. AND esi_wv LT 7924.)
  dwv = esi_wv[a[1]]-esi_wv[a[0]]
  print, 'No smooth: ',total( 1-esi_fx[a] )*dwv

  a = where(wave GT 7900. AND wave LT 7924.)
  dwv = wave[a[1]]-wave[a[0]]
  print, 'Smooth: ', total( 1-fx[a] )*dwv

; Plot
;  x_specplot, fx, sig, wave=wave, inflg=4, /qal, /block

  return
end
