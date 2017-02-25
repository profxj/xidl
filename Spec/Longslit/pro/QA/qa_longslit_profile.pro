;+
; NAME:
;   qa_longslit_profile
;
; PURPOSE:
;   Generate QA for the object profile fitting
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OUTPUTS:
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
;   20-Apr-2005  Written by J. Hennawi Berkeley
;-
;------------------------------------------------------------------------------

pro qa_longslit_profile, x_tot, y_tot, model_tot, l_limit, r_limit, ind = ind1 $
                         , title = title, XTRUNC = XTRUNC1, XRANGE = XRANGE1 $
                         , ORIG_MODEL = ORIG_MODEL_TOT, YRANGE = YRANGE1

  clr = getcolor(/load)
  green = clr.green
  red = clr.red
  blue = clr.blue
  white = clr.white
  cyan = clr.cyan
  lightcyan = clr.TURQUOISE
  yellow = clr.yellow

  IF KEYWORD_SET(XRANGE1) THEN XRANGE = XRANGE1 
  IF KEYWORD_SET(XTRUNC1) THEN XTRUNC = XTRUNC1 ELSE XTRUNC = 1d6
  IF n_elements(ind1) GT 0 THEN ind = ind1 $
  ELSE ind = lindgen(n_elements(x_tot))

   x = x_tot[ind]
   y = y_tot[ind]
   model = model_tot[ind]
   max_model = max(model) <  2.0D
   IF KEYWORD_SET(ORIG_MODEL_TOT) THEN ORIG_MODEL = orig_model_tot[ind]
   IF KEYWORD_SET(XRANGE) THEN BEGIN
      minx = -XRANGE
      maxx = XRANGE
   ENDIF ELSE BEGIN
      goodpix = WHERE(model GT 0.001*max_model, ngood)
      IF ngood GT 0 THEN BEGIN
         minx = 1.3*min(x[goodpix]) > (-XTRUNC)
         maxx = 1.3*max(x[goodpix]) < XTRUNC
      ENDIF ELSE BEGIN
         minx = -5.0D
         maxx = 5.0D
      ENDELSE
   ENDELSE


   nsamp = 150L
   half_bin = (maxx - minx)/nsamp/2.

   IF KEYWORD_SET(YRANGE1) THEN YRANGE = YRANGE1 $
   ELSE BEGIN
      ymax = (max(model)*1.5 > 0.3)
      ymin = -0.1*ymax < (-0.05)
      YRANGE = [ymin, ymax]
   ENDELSE
   plot, [0], [0], yr = yrange, xr = [minx, maxx], /nodata, title = title, $
       /xs, /ys
   plot_mid = (findgen(nsamp)+0.5)/nsamp * (maxx - minx) + minx

   y20 = fltarr(nsamp)
   y80 = fltarr(nsamp)
   y50 = fltarr(nsamp)
   model_samp = fltarr(nsamp)
   nbin = fltarr(nsamp)
   
   for i = 0L, nsamp - 1L do begin
     dist = abs(x-plot_mid[i])
     close = where(dist LT half_bin, nclose)
     nbin[i] = nclose
     if close[0] NE -1 then begin
         closest = min(dist[close], pl)
         model_samp[i] =model[close[pl]] 
     endif
     if nclose GT 3 then begin
         s = sort(y[close])
         y50[i] = y[close[s[(nclose-1)*0.5]]]
         y80[i] = y[close[s[(nclose-1)*0.8]]]
         y20[i] = y[close[s[(nclose-1)*0.2]]]
         oplot, [plot_mid[i], plot_mid[i]], [y20[i], y80[i]], thick = 2.0 $
                , color = lightcyan
     endif
  endfor
   icl = WHERE(nbin GT 3, ngd)
   IF ngd GT 0 THEN $
      oplot, plot_mid[icl], y50[icl], ps = 2, thick = 3, color = green $
   ELSE oplot, plot_mid, y50, ps = 2, thick = 3, color = green
   IF KEYWORD_SET(ORIG_MODEL_TOT) THEN $
      oplot, x[sort(x)], orig_model[sort(x)], thick = 1, color = blue
   oplot, x[sort(x)], model[sort(x)], thick = 3, color = red 
   
   ;;oplot, x_tot[sort(x_tot)], orig_model[sort(x_tot)], thick = 1 $
   ;;          , color = blue
   ;;oplot, x_tot[sort(x_tot)], model_tot[sort(x_tot)], thick = 3,
   ;;color = red 
   

   oplot, x, y, color = white, psym = 3

      
   if keyword_set(l_limit) then $
      oplot, [l_limit, l_limit], yrange, color = green
   if keyword_set(r_limit) then $
      oplot, [r_limit, r_limit], yrange, color = green

   return
end
