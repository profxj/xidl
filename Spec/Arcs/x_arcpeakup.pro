;+ 
; NAME: 
; x_arcpeakup   
;    Version 1.0
;
; PURPOSE:
;    ID's lines from an initial fit + a lines structure
;
; CALLING SEQUENCE:
;   
;   x_arcpeakup, spec, lines, guessfit, /FFT
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
;   x_arcpeakup, spec, lines, guessfit
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   26-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------

  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_arcpeakup, spec, lines, FFIT=ffit, WV=wv, NFRST=nfrst, NFIN=nfin, $
                 INTER=inter, SIG=sig

;  Error catching
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'x_arcpeakup, spec, lines, FFIT=, WV=, NFRST=, /INTER [v1.0]'
    	return
  endif 

; Optional Keywords
  npix = n_elements(spec)
  if not keyword_set( NFRST ) then nfrst = 3L
  if not keyword_set( NFIN ) then nfin = 5L
  if not keyword_set( MXOFF ) then mxoff = 3L
  if not keyword_set( SIG ) then sig = [3.,3.]

; FIRST fit
  frst_ffit = x_setfitstrct(NORD=nfrst, FLGREJ=1)
  gdlin = where(lines.flg_plt EQ 1)
  frst_fit = x_fitrej(lines[gdlin].pix,lines[gdlin].wave, FITSTR=frst_ffit)

; Calc for all lines
  x_fndpeaks, spec, peak, NSIG=3., /silent
  npk = n_elements(peak)
  wave = x_calcfit(dindgen(npix), FITSTR=frst_ffit)

; ID lines
  ;; Loop on peaks
  for i=0L,npk-1 do begin
      ;; Require that a calibration line lies within +/- 3 pixels
      pk = round(peak[i])
      gdln = where( (lines.wave - wave[(pk-mxoff)>0])* $
                    (lines.wave - wave[(pk+mxoff)<(npix-1)]) LE 0., $
                    nmatch)
      if nmatch NE 0 then begin
          ;; Find the closest calibration line
          sep = abs(lines[gdln].wave - wave[pk])
          minsep = min(sep, imin)
          lines[gdln[imin]].flg_plt = 1 
          ;; Center up on the peak
          lines[gdln[imin]].pix = peak[i] 
      endif
  endfor

; FINAL fit
  ffit = x_setfitstrct(NORD=nfin, FLGREJ=1)
  ffit.lsig = sig[0]
  ffit.hsig = sig[1]

  gdlin = where(lines.flg_plt EQ 1)
  fin_fit = x1dfit(lines[gdlin].pix,lines[gdlin].wave, FITSTR=ffit,$
                  INTER=inter)

  print, 'x_arcpeakup: RMS of solution is ', ffit.rms, ' Ang'

; RETURN
  if arg_present(WV) then wv = x_calcfit(dindgen(npix),FITSTR=ffit)

  return
end
