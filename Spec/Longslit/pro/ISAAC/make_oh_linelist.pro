;+
; NAME:
;   tspec_mk_oh_linelist
;
; PURPOSE:
;   Generate an OH line list from Rouss for TSpec
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Jan-2013  Written by JXP
;-
;------------------------------------------------------------------------------
pro make_oh_linelist, instrument, modelfile, linefile
  
  if  N_params() LT 2  then begin 
     print, 'Syntax - ' + $
            'tspec_mk_oh_linelist, modelfile, linefile'
     return
  endif   

  CASE INSTRUMENT OF
     'ISAAC': BEGIN
        slit_width = 0.6d ;; hard wiring everything for 0.6" slit for now
        plate_scale = 0.147d
        fnslit =slit_width/plate_scale
        pkwdth = 1.3d*nslit
        toler = nslit/2.0d > 2.0d
        THIN = 1
        FWEIGHT = 0
        NSIG = 5.0d
     END
     'TRIPLESPEC': BEGIN
        ;; Fill in analgous Tspec parameters
     END
     ELSE: message, 'Unknown instrument?'
  ENDCASE
  
  ;; Read in sky model of Rousselot linelist 
  model = xmrdfits(modelfile, 0, head)
  npix = n_elements(model)
  crval1 = sxpar(head, 'CRVAL1')
  cdelt1 = sxpar(head, 'CDELT1')
  dcflag = sxpar(head, 'DC-FLAG')
  IF dcflag EQ 1 THEN wave = 10.d^(crval1 + cdelt1*dindgen(npix)) $
  ELSE wave = crval1 + cdelt1*dindgen(npix) 
  
  ;; Define clean lineless regions for fitting continuum
  fit_wv = [ [0.8137525, 0.8244728],  $
             [0.8621566, 0.8739599], $
             [0.9137012, 0.9278867],$
             [1.060408, 1.063224], $  
             [1.079467, 1.082174], $
             [1.188598, 1.192172], $
             [1.207332, 1.210255], $
             [1.495439, 1.498038], $
             [1.535289, 1.537021], $
             [1.756367, 1.761998], $
             [1.862055, 1.867686], $
             [2.076419, 2.082700],$
             [2.136626, 2.146156], $ 
             [2.235341, 2.242704], $   
             [2.256565, 2.263712], $ 
             [2.279305, 2.288401], $    
             [2.332366, 2.343627], $    
             [2.369833, 2.376170]]  
  szf = size(fit_wv,/dimen)
  gdi = [0]
  for ii=0L,szf[1]-1 do begin
     morei = where( wave GT fit_wv[0,ii] AND wave LT fit_wv[1,ii], nmore)
     if nmore GT 0 then gdi = [gdi, morei]
  endfor
  
  
  afunc = 'LEGEND'
  ngd = n_elements(gdi)
  tmpfit = x_fitrej(wave[gdi], model[gdi], afunc, 8, $
                    ivar=replicate(1., ngd), NITER=1, $
                    hsigma=3., lsigma=3., rms=rms, FITSTR=fitstr, $
                    /SILENT)
  autofit = x_calcfit(wave, FITSTR = fitstr)

  x_fndpeaks, model, peak, NSIG = NSIG, /silent, PKWDTH = pkwdth, THIN = THIN $
              , FWEIGHT = FWEIGHT, TOLER = TOLER, AUTOFIT = autofit, RMS = rms
  ;print, peak[0:10]
  ;x_splot, model, xtwo=peak, ytwo=model[round(peak)], /bloc, psym2=1
  ;stop

  pkwave = interpol(wave, dindgen(n_elements(wave)), peak)
  m_pkwave = interpol(model, dindgen(n_elements(wave)), peak)
  x_splot, wave, model, xtwo=pkwave, ytwo=m_pkwave, /bloc, psym2=1

  label= replicate(' OH', n_elements(peak))
  writecol, linefile, pkwave, m_pkwave, label

  return
end

;------------------------------------------------------------------------------
