;+ 
; NAME:
; esi_lwdfitarc
;     Version 1.0
;
; PURPOSE:
;    Fits an Arc file to prepare for the wavelength map
;
; CALLING SEQUENCE:
;   
;  esi_lwdfitarc, esi, slit, /INTER, LINLIST=, REFROW=, CHK=
;
; INPUTS:
;   esi   -  ESI structure
;   indx  -  Indices of objects to process
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   FLAT  - Flat file
;   BIAS  - Bias frame
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_lwdfitarc, esi 
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro esi_lwdfitarc, esi, slit, LINLIST=linlist, INTER=inter, REFROW=refrow, CHK=chk

  COMMON ps_common, old_dname, old_pfont, file_name, opened, color, ledger
;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_lwdfitarc, esi, slit, LINLIST=, /INTER [v1.0]'
      return
  endif 

;  Optional Keywords
  
  if not keyword_set(LINLIST) then $
    linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/XeHgNe.lst'
  if not keyword_set(REFROW) then refrow = 820L

; Slit
  c_s = esi_slitnm(slit)


; Open Arc
  arc_fil = 'Arcs/ArcLWD_'+c_s+'.fits'
  img_arc = mrdfits(arc_fil, /silent)

  ;; Grab spectrum near row 820
  arc_spec = djs_median( img_arc[*,refrow-4:refrow+4], 2)
  npix = n_elements(arc_spec)
  xdat = findgen(npix)

  ;; Identify
  if keyword_set( INTER ) then begin
      x_identify, arc_spec, finfit, LINELIST=linlist
  endif else begin
      finfit = {fitstrct}
      finfit.func = 'POLY'
      finfit.nord = 6
      finfit.lsig = 2.5
      finfit.hsig = 2.5
      finfit.niter = 3
      finfit.minpt = 5
      finfit.maxrej = 5
      
      ;; Lines
      x_arclist, linlist, lines, /GDONLY

      ;; Peaks
      x_fndpeaks, arc_spec, peak, NSIG=4.
      npk = n_elements(peak)

      ;; Input wave soln
      if not keyword_set(ESILWDFIT) then $
        esilwdfit = getenv('XIDL_DIR')+ $
        '/ESI/pro/LWD/Arcs/AFIT_LWD10.fits'
      x_fitstrtofits, arcfit, esilwdfit, /reverse
      wave = x_calcfit(findgen(npix), FITSTR=arcfit)
      wave[0:475] = 0.

      ;; Add to lines
      for i=0L,npk-1 do begin
          ;; Require that a calibration line lies within +/- 3 pixels
          gdln = where( (lines.wave - wave[(peak[i]-3)>0])* $
                        (lines.wave - wave[(peak[i]+3)<(npix-1)]) LE 0., nmatch)
          if nmatch NE 0 then begin
              ;; Find the closest calibration line
              sep = abs(lines[gdln].wave - wave[peak[i]])
              minsep = min(sep, imin)
              lines[gdln[imin]].flg_plt = 1 
              ;; Center up on the peak
              center = peak[i]
              lines[gdln[imin]].pix = center 
              lines[gdln[imin]].wave_fit = x_calcfit(center, FITSTR=arcfit)
          endif
      endfor
      ;; Fit 
      gdfit = where(lines.flg_plt EQ 1, ngd)
      if ngd EQ 0 then begin
          print, 'esi_lwdfitarc: No good lines! Problem here'
          stop
      endif
      fit = x_fitrej(lines[gdfit].pix, lines[gdfit].wave, REJPT=rejpt, $
                     FITSTR=finfit)
      ;; PS File
      if keyword_set(opened) then begin
          if opened EQ 1 then ps_close, /noprint, /noid
      endif
      device, get_decomposed=svdecomp
      device, decompose=0
      filnm = 'Arcs/AFIT_LWD'+c_s+'.ps'
      ps_open, file=filnm, font=1, /color
      clr = getcolor(/load)
      plot, lines[gdfit].wave, fit-lines[gdfit].wave, psym=1, $
        charsize=1.8, $
        background=clr.white, color=clr.black, $
        xtitle='Wave', ytitle='Residual (Ang)', $
        xmargin=[12,2], ymargin=[6,2]
      if rejpt[0] NE -1 then $
        oplot, lines[gdfit[rejpt]].wave, $
        fit[rejpt]-lines[gdfit[rejpt]].wave, psym=2, color=clr.red
      xyouts, 0.2, 0.83, $
        'RMS = '+string(finfit.rms, FORMAT='(f5.3)'), /normal, $
        charsize=2.0
      ps_close, /noprint, /noid
      device, decomposed=svdecomp
  endelse

  ;; RMS
  print, 'esi_lwdfitarc: RMS = ', finfit.rms, ' Ang' 

  ;; CHK
  if keyword_set( CHK ) then begin
      wv = x_calcfit(findgen(2048), FITSTR=finfit)
      x_splot, wv, arc_spec, /block
  endif

; Output
  outfil = 'Arcs/AFIT_LWD'+c_s+'.fits'
  x_fitstrtofits, finfit, outfil

  print, 'esi_lwdfitarc: All done!'

  return
end

