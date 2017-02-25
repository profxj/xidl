;+ 
; NAME:
; dla_guessew   
;   Version 1.1
;
; PURPOSE:
;    Calcualtes the EW limits for a series of input transitions for a
;    DLA with given NHI and metallicity Z.
;
; CALLING SEQUENCE:
;   dla_guessew, fil, NHI, Z, NPIX=, SNR=, NSIG=, DWV=
;
; INPUTS:
;  fil -- File with list of wavlengths
;  NHI -- N(HI) value of the DLA (logarithmic)
;  Z   -- Metallcitiy  (logarithmic)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  NPIX=  -- Number of pixels making up the feature [default: 4]
;  NSIG=  -- Number of sigma significant [default: 3.]
;  SNR=   -- Signal-to-noise per pixel [default: 10.]
;  DWV=   -- Width of pixel in Angstromgs [default: 0.1 Ang]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_guessew, fil, NHI, Z
;
;
; PROCEDURES/FUNCTIONS CALLED:
; REVISION HISTORY:
;   14-Nov-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro dla_guessew, fil, NHI, Z, NPIX=npix, SNR=snr, NSIG=nsig, $
               DWV=dwv

;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'dla_guessew, fil, NHI, Z, NPIX=, SNR=, NSIG=, ZABS=, DWV= [v1.1]'
    return
  endif 

; Optional Keywords
  if not keyword_set( NSIG ) then nsig = 3.
  if not keyword_set( NPIX ) then npix = 4
  if not keyword_set( SNR ) then snr = 10.
  if not keyword_set( DWV ) then dwv = 0.1  ;; Ang

; Grab the transitions
  readcol, fil, wv
  nwv = n_elements(wv)

; Calculate M column
  M = NHI + Z

; Calculate limiting EW
  EWlmt = sqrt(npix) * dwv / snr / nsig

  print, 'The limiting EW for S/N = ', SNR, ', NPIX= ', NPIX, ', NSIG= ', NSIG, $
    ', and ', DWV, 'A pixels is ', EWlmt*1000., 'mA'
  
  print, 'The following are rest EW (mA)!'
; Loop
  for q=0L,nwv-1 do begin

      ;; Grab name
      getfnam, wv[q], f, nm

      ;; Get Z
      getabnd, nm, Znum, abnd

      ;; Calculate expected N
      N = M - 12 + abnd

      ;; Calculate rest EW
      EW = 8.85249d-13 * 10^N * f * wv[q] * 1d-8 * wv[q]  ;; Ang

      print, strtrim(nm,2), wv[q], EW*1000., FORMAT='(a10, 1x,f9.3,1x,f7.2)'  
  endfor

end
