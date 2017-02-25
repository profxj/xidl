;+ 
; NAME:
; x_hca
;   Version 1.0
;
; PURPOSE:
;    Creates a plot which investigates hidden saturation
;
; CALLING SEQUENCE:
;   
; x_hca, datfil, z, vmnx, wvarr, YMNX=, PSFILE=, BW=
;
; INPUTS:
;   datfil - Data file
;   z      - redshift
;   vmnx   - velocity limits
;   wvarr  - Rest wavelengths of ions of interest
;
; RETURNS:
;   
; OUTPUTS:
;   Creates a Plot
;
; OPTIONAL OUTPUTS:
;  YMNX=  - Sets the y range of the plot
;  PSFILE= - Name of the plot file
;  /BW     - Make a B&W plot
;
; COMMENTS:
;
; EXAMPLES:
;   
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   3-Dec-2005 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro x_hca, datfil, z, vmnx, wvarr, YMNX=ymnx, PSFILE=psfile, BW=bw

  if  N_params() LT 4  then begin 
      print, 'Syntax - ' +$
        'x_hca, datfil, z, vmnx, wvarr, YMNX=, PSFILE=, /BW  (v1.1)'
      return
  endif 
;
  ;; Optional keywords
  cst = (10.d^14.5761)


  ;; Read data
  fx = x_readspec(datfil, inflg=inflg, wav=wav, sig=sig, /auto)

  ;; Velocity
  velo = x_allvelo(wav, z, wvarr, vmnx, all_pmnx=pmnx, NPIX=5000L)

  ;; Set YMAX
  if not keyword_set(YMNX) then begin
      ymnx = fltarr(2)
      getfnam, wvarr[0], fval
      pix = pmnx[0,0] + lindgen(pmnx[2,0])
      tau = alog(1. / (fx[pix] > sig[pix]) )
      Na = tau * cst / fval / wvarr[0]
      ymnx[1] = max(Na)*1.1
  endif

  norm = 10.^round(alog10(ymnx[1]))
  ymnx = ymnx / norm

  ;; PS FILE
  if keyword_set(PSFILE) then x_psopen, psfile, /maxs
  clr = getcolor(/load)
 

  plot, [0],  [0],  psym=psym, xrange=vmnx, $
    yrange=ymnx, color=clr.black, $
    xmarg=[10, 1], ymarg=[5,1], charsiz=1.8, $
    background=clr.white, ystyle=1, xstyle=1, $
    xtitle='Velocity (km/s)', ytitle='N!da!N (10!u'+ $
    string(round(alog10(norm)),format='(i2)')+'!N cm!u-2!N)', /nodata

  xcolors = x_setclrs()
  for qq=0L,n_elements(wvarr)-1 do begin
      getfnam, wvarr[qq], fval
      pix = pmnx[0,qq] + lindgen(pmnx[2,qq])
      tau = alog(1. / (fx[pix] > (0.5*sig[pix])) )
      Na = tau * cst / fval / wvarr[qq] / NORM
      if not keyword_set(BW) then $
        oplot, velo[*,qq], Na, color=xcolors[qq], psym=10 $
      else oplot, velo[*,qq], Na, psym=10, linestyle=qq
      ;; Label
      xyouts, vmnx[0] + (vmnx[1]-vmnx[0])*0.7, $
        ymnx[1]*(0.8 - qq*0.05), string(wvarr[qq],format='(f8.3)'), $
        color=xcolors[qq], charsize=1.5
  endfor

  ;;
  if keyword_set(PSFILE) then x_psclose

return
end

