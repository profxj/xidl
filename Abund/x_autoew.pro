;+ 
; NAME:
; x_autoew   
;   Version 1.0
;
; PURPOSE:
;    Measures the EW of a set of lines passed in a structure
;
; CALLING SEQUENCE:
;   
;   x_autoew, fil, strct, /PLOT, LSNR=, RESSIG=, INFLG=, CUTWV=
;
; INPUTS:
;  filnm -- Filename of the spectrum
;
; RETURNS:
;
; OUTPUTS:
;  strct -- Structure containing info on all the lines
;
; OPTIONAL KEYWORDS:
;  LSNR   --  SNR for auto-id of lines (default: 5)
;  /PLOT  -- Show the fits
;  SIGMA  -- Gaussian sigma of the lines.  Should match the resolution
;  INFLG  -- Flag to pass to x_readspec for reading the data file
;  CUTWV= -- Cut out wavelength regions at end of the spectrum
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
; REVISION HISTORY:
;    Written by JXP
;-
;------------------------------------------------------------------------------

pro x_autoew, filnm, strct, PLOT=plot, LSNR=lsnr, RESSIG=RESSIG, INFLG=inflg, $
              CUTWV=cutwv

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'x_autoew, filnme, strct, /PLOT, LSNR=, RESSIG=, INFLG= [v1.1]'
    return
  endif 

  if not keyword_set( LSNR ) then lsnr = 5.
  if not keyword_set(RESSIG) then RESSIG = 3.

  fx = x_readspec(filnm, inflg=INFLG, wav=wv, sig=sig)
  if keyword_set(PLOT) then clr = getcolor(/load)

  ;; Find lines
  x_findgauss, 1.-fx, 1./sig^2, xpeak=xpeak, ypeak=ypeak, $
               xsn=xsn, sn=sn, nfind=300L, SIGMA=RESSIG
  good = where(sn GT lsnr,ngood)
  if ngood NE 0 then begin
      srt = sort(sn[good])
      gdd = good[srt]
      ;; Cut out endpoints of spectrum
      if keyword_set(CUTWV) then begin
          fgood = where(wv[round(xpeak[gdd])] LT CUTWV[1] and $
                        wv[round(xpeak[gdd])] GT CUTWV[0])
          gdd = gdd[fgood]
      endif
      printcol, wv[round(xpeak[gdd])], sn[gdd]
      ngd = n_elements(gdd)
      ;; Recenter to closest pixel
      for jj=0L,ngd-1 do begin
          pix = round(xpeak[gdd[jj]]) + lindgen(5) - 2
          mn = min(fx[pix],imn)
          xpeak[gdd[jj]] = pix[imn]
      endfor
  endif 

  tmp = { $
        dfil: '', $
        flg: 0, $
        obswv: 0.d, $
        nam: '', $
        wrest: 0.d, $
        zabs: 0.d, $
        ew: 0.,  $
        sigew: 0. $
        }
  strct = replicate(tmp, ngd) 
  strct.dfil = filnm
  

  ;; Fit with a Gaussian
  for qq=0L,ngd-1 do begin
      cpix = round(xpeak[gdd[qq]])
      obswv = wv[cpix]

      ;; Identify pixels for fitting
      poss_pix = cpix + lindgen(51) - 21
      cti = where(fx[poss_pix] GT 0.95, ncti)
      if ncti EQ 0 then stop
      dff = min(abs(cpix - poss_pix[cti])) > round(RESSIG*2)
      dff = dff + 1
;      print, qq, dff
      px = lindgen(dff*2 + 1) + cpix - dff

      ;; Profile
      gprof = -1.*(fx[px] - 1.)
      wcen = obswv
      dwv = abs(wv[cpix]-wv[cpix+1])
      gsssig = 4.

      ;; FIT
      yfit = gaussfit(wv[px], gprof, acoeff, $
                      estimates=[max(gprof), wcen, gsssig], $
                      sigma=sigma, nterms=3, $
                      measure_errors=sig[px])

      ;; Calc EW
      ewval = acoeff[0] * acoeff[2] * sqrt(!pi*2.) ; A
      sigew1 = sqrt(total( ((sig[px])*dwv)^2)) ; A
      sigew2 = sqrt(!pi*2.) * sqrt( (acoeff[0]*sigma[2])^2 + $
                                    (acoeff[2]*sigma[0])^2 )
;      print, wv[cpix], sigew1, sigew2

      ;; Write to Struct
      strct[qq].obswv = acoeff[1]
      strct[qq].ew = ewval
      strct[qq].sigew = sigew1 > sigew2

      ;; Quick plot
      if keyword_set(PLOT) then begin
          ppx = lindgen(51) + cpix - 25
          plot, wv[ppx], fx[ppx], color=clr.black, background=clr.white, psym=10
;          oplot, wv[ppx], conti[ppx], color=clr.red
          oplot, wv[px], 1.-yfit, color=clr.blue

          oplot, replicate(wcen,2), [-9e9,9e9], color=clr.black, linestyle=2, thick=3

          ;; Value
          xyouts, 0.1, 0.3, 'EW = '+string(strct[qq].ew, format='(f5.3)'), $
                  /normal, charsize=lsz, color=clr.black
          xyouts, 0.1, 0.2, '!9s!X(EW) = '+string(strct[qq].sigew, format='(f5.3)'), $
                  /normal, charsize=lsz, color=clr.black
          wait, 0.5
      endif
  endfor

  return
end
