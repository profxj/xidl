;+ 
; NAME:
; fig_1darcfit
;     Version 1.1
;
; PURPOSE:
;-
;------------------------------------------------------------------------------

pro fig_1darcfit, WV=wv, FIT=fit, REJPT=rejpt,$
                  ORDR=ordr, RMS=rms, FORDR=fordr, DWV=dwv

  if not keyword_set(PSFILE) then psfile='fig_1darcfit.ps'
  if not keyword_set(ARCFIL) then arcfil = $
    getenv('MIKE_PAP')+'Arcs/Fits/mb0005_fit.idl'
  if not keyword_set(OSTRFIL) then ostrfil = $
    getenv('MIKE_PAP')+'Flats/OStr_B_01.fits'
  if not keyword_set(CSIZE) then csize = 2.4
  if not keyword_set(LSZ) then lsz = 1.2

  ymnx= [-0.5,0.5]

  restore, arcfil
  ostr = xmrdfits(ostrfil, 1, /silent)

  x_psopen, psfile, /maxs
  !p.multi=[0,3,2]
  clr = getcolor(/load)

  sz_arc = [1024L,2048]

  ;; Loop
  nordr = 6L
  for ii=0L,nordr-1 do begin
      ;; Order
      ordr = ostr[ii].order

      if rejstr[ii].ngdf EQ 0L then continue
      gdfit = rejstr[ii].gdfpt[0:rejstr[ii].ngdf-1]
      ;; Rejpt
      if rejstr[ii].nrej NE 0 then rejpt = rejstr[ii].rejpt[0:rejstr[ii].nrej-1] $
      else rejpt = -1L
      ;; Dwv
      wv = 10^x_calcfit(dindgen(sz_arc[1]),fitstr=all_arcfit[ii])
      dwv = abs(median(wv - shift(wv,1)))
      ;; Subroutine 
      fit=10^x_calcfit(double(rejstr[ii].gdfpx[0:rejstr[ii].ngdf-1]), $
                       FITSTR=all_arcfit[ii])

      ;; Set the wavelengths
      wv = sv_lines[ii].wv[0:sv_lines[ii].nlin-1]
      if rejstr[ii].nrej NE 0 then begin
          wv = [wv, rejstr[ii].rejwv[0:rejstr[ii].nrej-1]]
          wv = wv[sort(wv)]
      endif
          

      ;; All points
      if n_elements(wv) LE 1 then NODATA=1
      mn = min(wv, max=mx)
      xlbl = mn + (mx-mn)*0.05

      plot, [wv], [(fit-wv)/dwv], psym=1, $
            background=clr.white, color=clr.black, $
            xtitle='Wavelength (Ang)', ytitle='Residual (Pix)', $
            xmargin=[8,1], $
            ymargin=[4,0], yrange=[-0.5, 0.5], xrange=[mn-1,mx+1.], $
            xstyle=1, ystyle=1, xtickinterval=(fix(wv[0]/1000.)-1)*10, $
            NODATA=nodata, charsize=csize
      ;; MASK 
      if rejpt[0] NE -1 then $
        oplot, [wv[rejpt]], [(fit[rejpt]-wv[rejpt])/dwv], psym=2, color=clr.red
      ;; Order

      ;; 
      oplot, [-9e9, 9e9], [0., 0.], linestyle=1, color=clr.blue
      
      xyouts, xlbl, ymnx[1]-0.1, 'Order = '+string(ordr, FORMAT='(i3)')+$
              '   !9Dl!X = '+string(dwv, FORMAT='(f6.4)'), charsize=lsz
      RMS=all_arcfit[ii].rms*sv_lines[ii].wv[0]*alog(10.) / dwv  ; pix
      xyouts, xlbl, ymnx[0]+0.1, 'RMS(pix)='+string(rms, FORMAT='(f6.4)'),$
              charsize=lsz
;              '   Fit Ord='+string(all_arcfit[ii].nord, $
;                                   FORMAT='(i2)'), charsize=lsz
      
  endfor
  x_psclose
  !p.multi=[0,1,1]

  return
end

