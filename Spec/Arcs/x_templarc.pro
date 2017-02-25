;+ 
; NAME: 
; x_templarc   
;    Version 1.1
;
; PURPOSE:
;    Given an archived Arc spectrum and solution and a new arc
;    spectrum, find the pixel offset and auto-id the arc lines in new
;    arc spectrum.  One must turn /FFT on to calculate the shift.
;
; CALLING SEQUENCE:
;   
; x_templarc, spec, lines, guessfit, /FFT, MSK=, MXSHFT=, $
;               MOCK_FFT=, SHFT=, ALL_PK=, PKWDTH=, $
;               /THIN, FORDR=, PKSIG=, SKIPLIN=, $
;               MXOFF=, /LOG
;
; INPUTS:
;   spec       - Input arc image or spectrum
;   lines      - Arc line structure
;  guessfit    - FIT to archived arc
;
; RETURNS:
;
; OUTPUTS:
;   lines     - Sets flg_plt to 1 those lines which are ID'd
;
; OPTIONAL KEYWORDS:
;  /FFT    -- Calculate the most likely offset between the archived
;             arc and the new one using the FFT procedure
;  MOCKFFT -- Archived FFT [default: create one from the line list]
;  MXOFF=  -- Maximum offset between ID line and a line in the
;             linelist  [default: 3 pixels]
;  SHFT=   -- Shift between archived and new arc (Input and/or Output)
;  PKSIG=  -- Sigma significant of arc lines (as found by x_fndpeak)
;             [default: 7.]
;  MSK=    -- Mask of the lines
;  MXSHFT= -- Maximum shift between archived and new solution
;  /THIN   -- Allow for very narrow arc lines
;  /LOG    -- Do the algorithm with logarithmic wavelengths
;  FORDR   -- Order for the fit for the continuum of the arc 
;  /SKIPLIN -- Do not bother to identify good lines for subsequent
;              fitting (algorithm has only calucated the shift)
;
; OPTIONAL OUTPUTS:
;  ALL_PK= --  Pix values of the arc line peaks
;
; COMMENTS:
;
; EXAMPLES:
;   x_templarc, spec, lines, guessfit, SHFT=shft
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_fndpeaks
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
                MXOFF=mxoff, LOG=log, FLG=flg, BEST=best, toler = toler, $
                MAXQUAL = MAXQUAL, FWEIGHT = FWEIGHT, ICLSE = ICLSE, FGAUSS = fgauss, $
                AUTOFIT = AUTOFIT, STRETCH=stretch, DEBUG=debug

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
  IF NOT KEYWORD_SET(TOLER) THEN TOLER = 2.0D
  IF NOT KEYWORD_SET(MAXQUAL) THEN MAXQUAL = 99999L
  flg = 1

; Make wave array
  if keyword_set (STRETCH) then begin
     twave = x_calcfit(dindgen(npix), FITSTR=guessfit)
     wave = congrid(twave, npix+stretch, /interp)
     wave = wave[0:npix-1]
  end else begin
     wave =  x_calcfit(dindgen(npix), FITSTR=guessfit)
  endelse
  if keyword_set( LOG ) then begin
      ;; Check
      if median(wave) LT 20. then wave = 10^wave
  endif

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
  wave[npix-1-abs(shft):*] = 0.

  ;; Min and max
  wvmn = min([wave[abs(shft)+1], wave[npix-1-abs(shft)-1]],max=wvmx)
  gdwv = where(wave GT 0.,ngdwv)
  if ngdwv EQ 0 then stop
  mnpix = gdwv[0]
  mxpix = gdwv[ngdwv-1]
  ;; Find peaks JFH added TOLER=TOLER
  x_fndpeaks, spec, peak, NSIG = pksig, /silent, PKWDTH = pkwdth $
              , THIN = thin, MSK = msk, NORDB = fordr, TOLER = TOLER $
              , FWEIGHT = FWEIGHT, ICLSE = ICLSE, FGAUSS = fgausss, AUTOFIT = AUTOFIT
  if peak[0] EQ -1 then begin
      print, 'x_templarc:  No peaks found!  Returning..'
      flg = -1
      return
   endif

  ;; Mask
  rnd_pk = round(peak)
  peak = peak(where(msk[rnd_pk] EQ 1))
  npk = n_elements(peak)

  ;; Save all peaks
  if arg_present(ALL_PK) then all_pk = peak

  if keyword_set( SKIPLIN ) then return

  ;; Setup LINES
  gdlin = where(lines.wave GT wvmn AND $
                lines.wave LT wvmx AND lines.flg_qual GT 0 AND $
                lines.flg_qual LE MAXQUAL,ngdlin)

; ID lines
  ;; Loop on peaks
if (ngdlin ne 0L) then begin ; jm08jun10nyu
  for i=0L,npk-1 do begin
      ;; Require that a calibration line lies within +/- 3 pixels
      pk = round(peak[i])
      if pk GT mxpix or pk LT mnpix then continue
      ;;
      gdln = where( (lines[gdlin].wave - wave[(pk-mxoff)>mnpix])* $
                    (lines[gdlin].wave - wave[(pk+mxoff)<mxpix]) LE 0., $
                    nmatch)
      if nmatch NE 0 then begin
          gdln = gdlin[gdln]
          ;; Find the closest calibration line
          sep = abs(lines[gdln].wave - wave[pk])
          minsep = min(sep, imin)
          lines[gdln[imin]].flg_plt = 1 
          ;; Center up on the peak
          lines[gdln[imin]].pix = peak[i] 
;          IF i EQ npk-1L THEN stop
;          lines[gdln[imin]].wave_fit = x_calcfit(peak[i], FITSTR=guess_fit[qq])
      endif
   endfor
endif 
  if keyword_set(DEBUG) then stop


  ;; Restrict to the 'best' 10 lines
  if keyword_set( BEST ) then begin
      gdl = where(lines.flg_plt EQ 1, ngdl)
      msk = bytarr(ngdl)
      gd = where(spec GT 0., ngd)
      mn = min(gd, max=mx)
      ;; Split the spectrum in to best/2 sections
      nbst = (mx-mn) / (best/2)
      for jj=0L,(best/2)-1 do begin
          imn = (mn + jj*nbst) > (abs(shft)+1)
          imx = (mn + (jj+1)*nbst) < (npix-2-abs(shft))
          if wave[imn] GT wave[imx] then begin
              tmp = imn
              imn = imx
              imx = tmp
          endif
          ;; Find lines
          gdl_s = where(lines[gdl].wave GT wave[imn] AND $
                      lines[gdl].wave LT wave[imx], nng)
          print, imn, imx, wave[imn], wave[imx], nng
          case nng of 
              0: 
              1: msk[gdl_s] = 1B
              2: msk[gdl_s] = 1B
              else: begin
                  pkf = fltarr(nng)
                  for kk=0L,nng-1 do $
                    pkf[kk] = max(spec[lines[gdl_s[kk]].pix -1 + lindgen(3)])
                  srt = sort(pkf)
                  msk[gdl_s[srt[nng-2:nng-1]]] = 1B
              end
          endcase
      endfor
      bgd = where(msk EQ 0, nb)
      if nb NE 0 then lines[gdl[bgd]].flg_plt = 0
  endif
      
  return
end
