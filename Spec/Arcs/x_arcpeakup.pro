;+ 
; NAME: 
; x_arcpeakup   
;    Version 1.1
;
; PURPOSE:
;    Finds the fit iteratively to a set of lines given the Arc
;    spectrum and a line list.  This code has been superseded by
;    x_arctempl.
;
; CALLING SEQUENCE:
;   x_arcpeakup, spec, lines, FFIT=, WV=, NFRST=, NFIN=, /INTER, SIG=
;
; INPUTS:
;   spec       - Input arc spectrum
;   lines      - Arc line list structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  NFRST= -- Order of fit for the first try [default: 3L]
;  NFIN=  -- Order of fit for the final fit
;  MXOFF= -- Maximum pixel offset between a line and its expected
;            position
;  SIG=   -- Lower and upper sigma values for fitting
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
             'x_arcpeakup, spec, lines, FFIT=, WV=, NFRST=, /INTER [v1.1]'
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
