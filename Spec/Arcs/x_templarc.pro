;+ 
; NAME: 
; x_templarc   
;    Version 1.0
;
; PURPOSE:
;    ID's lines from an initial fit + a lines structure
;
; CALLING SEQUENCE:
;   
;   x_templarc, spec, lines, guessfit, /FFT
;
; INPUTS:
;   spec       - Input arc image or spectrum
;   lines      - Arc line structure
;
; RETURNS:
;
; OUTPUTS:
;   lines     - Sets flg_plt to 1 those lines which are ID'd
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_templarc, spec, lines, guessfit
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   26-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_templarc_gauss, x, x0, DISP=disp
  ;; Sigma
  if not keyword_set( DISP ) then disp = 3.

  ;; Function
  return, 10.*exp( -(x-x0)^2 / disp^2 )
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_templarc, spec, lines, guessfit, FFT=fft, MSK=msk, MXSHFT=mxshft, $
                MOCK_FFT=mock_fft, SHFT=shft, ALL_PK=all_pk, PKWDTH=pkwdth, $
                THIN=thin, FORDR=fordr, PKSIG=pksig, SKIPLIN=skiplin, $
                MXOFF=mxoff, LOG=log

;  Error catching
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'x_templarc, spec, lines, newfit, /FFT [v1.0]'
    	return
  endif 

; Optional Keywords
  npix = n_elements(spec)
  if not keyword_set( MXOFF ) then mxoff = 3L
  if not keyword_set( SHFT ) then shft = 0L
  if not keyword_set( PKSIG ) then pksig = 7.
  if not keyword_set( MSK ) then msk = bytarr(npix) + 1B
  if not keyword_set( MXSHFT ) then mxshft = 10L

; Make wave array
  wave = x_calcfit(dindgen(npix), FITSTR=guessfit)
  if keyword_set( LOG ) then wave = 10^wave
  wvmn = min(wave,max=wvmx)

; Get shift from FFT
  if keyword_set( FFT ) then begin
      if not keyword_set( MOCK_FFT ) then begin
          shft_lin = where(lines.flg_qual GT 3, nshft)
          if nshft LE 3 then begin
              print, 'x_templarc: Not enough lines for shift!'
              shift = 0L
          endif
          ;; Create mock spec
          pix = fltarr(nshft)
          mock = fltarr(npix)
          stop  ;; NEED TO UPDATE FOR LOG WAVELENGTHS
          for i=0L,nshft-1 do begin
              pix = x_fndfitval(lines[shft_lin[i]].wave, guessfit, $
                                dindgen(npix), wave)
              if msk[round(pix)] NE 1B then continue
              lhs = 0L > round(pix-10.)
              rhs = round(pix+10.) < (npix-1L)
              mock[lhs:rhs] = x_templarc_gauss(float(lhs)+findgen(rhs-lhs+1),pix)
          endfor
          ;; FFT
          fft_mock = fft(mock)
      endif else fft_mock = mock_fft
      fft_spec = fft(spec)
      corr = fft( fft_spec * conj(fft_mock), /inverse)
      ans = double(corr)
      tmp = shift(ans,mxshft)
      ;; Find max within 10pix
      mx = max(tmp[0:mxshft*2-1],imx)
      shft = imx-mxshft
      if abs(shft) EQ mxshft then stop
  endif else shft = shft[0]

  wave = shift(wave, round(shft))
  ;; Zero out
  wave[0:abs(shft)] = 0.
  wave[npix-1-abs(shft):npix-1] = 0.

; Find peaks
  x_fndpeaks, spec, peak, NSIG=pksig, /silent, PKWDTH=pkwdth, THIN=thin, $
    MSK=msk, NORDB=fordr

  ;; Mask
  rnd_pk = round(peak)
  peak = peak(where(msk[rnd_pk] EQ 1))
  npk = n_elements(peak)

  ;; Save all peaks
  if arg_present(ALL_PK) then all_pk = peak

  if keyword_set( SKIPLIN ) then return
; Setup LINES
  gdlin = where(lines.wave GT wvmn AND $
                lines.wave LT wvmx AND lines.flg_qual GT 0)

; ID lines
  ;; Loop on peaks
  for i=0L,npk-1 do begin
      ;; Require that a calibration line lies within +/- 3 pixels
      pk = round(peak[i])
      gdln = where( (lines[gdlin].wave - wave[(pk-mxoff)>0])* $
                    (lines[gdlin].wave - wave[(pk+mxoff)<(npix-1)]) LE 0., $
                    nmatch)
      if nmatch NE 0 then begin
          gdln = gdlin[gdln]
          ;; Find the closest calibration line
          sep = abs(lines[gdln].wave - wave[pk])
          minsep = min(sep, imin)
          lines[gdln[imin]].flg_plt = 1 
          ;; Center up on the peak
          lines[gdln[imin]].pix = peak[i] 
;          lines[gdln[imin]].wave_fit = x_calcfit(peak[i], FITSTR=guess_fit[qq])
      endif
  endfor

  return
end
