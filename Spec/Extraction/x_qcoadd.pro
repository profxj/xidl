;+ 
; NAME:
; x_qcoadd
;    Version 1.0
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   
;   x_qcoadd, wv1, fx1, wv2, fx2, wvt, wpt
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;   Updates slitstr for original positions
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_qskysub, slitstr, map
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   14-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_qcoadd, wv1, sp1, wv2, sp2, wvt, fxt


;  Error catching
  if  N_params() LT 5  then begin 
    print,'Syntax - ' + $
             'x_qcoadd, wv1, fx1, wv2, fx2, wvt, fxt [v1.0]'
    return
  endif 


;  Optional Keywords

;  Find offset for fx2

  tmp = (sp1 > 0.01) / (sp2 > 0.01)
  print, 'x_qcoadd: Fix this!'
  md = median(tmp[500:1500])  ; This will need to be fixed
  fx1 = sp1
  fx2 = sp2*md

;  Get all wv, fx

  npix1 = n_elements(fx1)
  npix2 = n_elements(fx2)
  allwv = reform( [wv1, wv2], npix1+npix2)
  allfx = reform( [fx1, fx2], npix1+npix2)

; BSPLIN

  fit = x1dfit(allwv, allfx, hsig=3., lsig=5., NORD=15, FUNC='BSPLIN', $
              FITSTR=fitstr)

; Start Looping
  fitstr.rms = sqrt(fitstr.rms)
  print, 'x_qcoadd: Should calculate the rms locally'
  msk1 = bytarr(npix1) + 1
  for i=0L,npix1-1 do begin
      fval = x_calcfit(wv1[i],FITSTR=fitstr)
      if abs(fval-fx1[i]) GT 5.*fitstr.rms then begin
          mn = min( abs(wv1[i]-wv2), imn)
          im = (imn-1) > 0
          ip = (imn+1) < (npix2-1)
          if fx1[i] GT fval then begin
              if fx2[im] GT 3*fitstr.rms+fval OR $
                fx2[imn] GT 3*fitstr.rms+fval OR $
                fx2[ip] GT 3*fitstr.rms+fval then continue else msk1[i] = 0
          endif else begin
              if fx2[im] LT fval - 3*fitstr.rms OR $
                fx2[imn] LT fval - 3*fitstr.rms OR $
                fx2[ip] LT fval - 3*fitstr.rms then continue else msk1[i] = 0
          endelse
      endif
  endfor
              
  ; Other spectrum
  msk2 = bytarr(npix1) + 1
  for i=0L,npix2-1 do begin
      fval = x_calcfit(wv2[i],FITSTR=fitstr)
      if abs(fval-fx2[i]) GT 5.*fitstr.rms then begin
          mn = min( abs(wv2[i]-wv1), imn)
          im = (imn-1) > 0
          ip = (imn+1) < (npix1-1)
          if fx2[i] GT fval then begin
              if fx1[im] GT 3*fitstr.rms+fval OR $
                fx1[imn] GT 3*fitstr.rms+fval OR $
                fx1[ip] GT 3*fitstr.rms+fval then continue else msk2[i] = 0
          endif else begin
              if fx1[im] LT fval - 3*fitstr.rms OR $
                fx1[imn] LT fval - 3*fitstr.rms OR $
                fx1[ip] LT fval - 3*fitstr.rms then continue else msk2[i] = 0
          endelse
      endif
  endfor

  ; Just the good pixels
  allmsk = reform( [msk1, msk2], npix1+npix2)
  gdpix = where(allmsk EQ 1, ngd)
  ; Sort
  srt = sort(allwv[gdpix])

  ; Final arrays
  wvt = allwv[gdpix[srt]]
  fxt = allfx[gdpix[srt]]

  return
end
