;+ 
; NAME:
; x_arcimage   
;    Version 1.0
;
; PURPOSE:
;    Creates an arcimage given a tracestructure.  This code has been
;    superseded by x_mkaimg.
;
; CALLING SEQUENCE:
;   
;   arcimg = x_arcimage( trcstrct, arcfit, sz_img )
;
; INPUTS:
;   trcstrct - Trace structure (contains all of the info on the
;              position and value of various arc lines)
;   arcfit   - Template fit structure ; Includes fit instructions
;   sz_img   - Size of output image
;   [lines]  - Line list strurcture (Usually defined)
;
; RETURNS:
;   arcimg
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  YSTRT  - Row defining the arcfit [default = middle]
;  LINELIST  - Arc line list (in lieu of lines)
;  NSIG  - Sig of RMS from the fit that line must match
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   arcimg = x_arcimage( trcstrct, arcfit, imgsz)
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_arcimage, trcstr, arcfit, sz_img, lines, YSTRT=ystrt, $
                     LINELIST=linelist, NSIG=nsig


;  Error catching
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'img = x_arcimage( trcstr, arcfit, sz_img, [lines], LINELIST=, '
    print, '         NSIG=) [v1.1]'
    return, -1
  endif 

  if n_elements(sz_img) NE 2 then begin
      print, 'x_arcimage: need a 2-D image'
      return, -1
  endif

;  Optional Keywords

  if not keyword_set( NSIG ) then nsig = 3.
  if not keyword_set( YSTRT ) then ystrt = sz_img[1]/2L
  if not keyword_set( LINES ) then begin
      if not keyword_set( LINELIST ) then begin
          print, 'x_arcimage: Need to define lines or linelist!'
          return, -1
      endif
      x_arclist, linelist, lines
  endif

;  ID the key lines; those that are traced

  wvcen = x_calcfit(trcstr.xfit[ystrt,*], FITSTR=arcfit)
  ntrc = n_elements(wvcen)
  gdflg = lonarr(ntrc) - 1L
  for i=0L,ntrc-1 do begin
      mn = min(abs(wvcen[i]-lines.wave),imn)
      if mn LT nsig*sqrt(arcfit.rms) then begin
          lines[imn].flg_plt = i
          gdflg[i] = imn
      endif
  endfor

  gdtrc = where(gdflg NE -1L, ngd)
  pix = fltarr(ngd)
  wav = lines[gdflg[gdtrc]].wave
          
;  FITS

  xx = findgen(sz_img[0])
  tmpfit = arcfit
  arcimg = dblarr(sz_img[0], sz_img[1])
  for i=0L,sz_img[1]-1 do begin
      ; Set pixels
      pix = trcstr.xfit[i,gdtrc]
      ; Fit using ARCFIT options (Assuming rejection)
      fit = x_fitrej(pix, wav, FITSTR=tmpfit)
      ; Add to arcimg
      arcimg[*,i] = x_calcfit(xx, FITSTR=tmpfit)
  endfor

; Free memory
  delvarx, xx, pix, wav

  return, arcimg
end
