;+
; NAME:
;   long_slits
;
; PURPOSE:
;   Simple code to find slit edges.  Obsolete?
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
pro lris_slits, flat, lhedg, rhedg


   ncol = (size(flat))[1]
   nrow = (size(flat))[2]
   saw = flat[0:ncol-2,*] - flat[1:*,*]

   saw[ncol/2-2:ncol/2+2,*] = 0

   saw_cross = median(djs_avsigclip(saw,2),3)

   flat_med = djs_avsigclip(smooth(flat,5),2)
   flat_med = 0.5*(flat_med[0:ncol-2] + flat_med[1:*])

   saw_ratio = saw_cross/(flat_med + (flat_med EQ 0)) * (flat_med GT 100)
   saw_ratio[0:3] = 0
   saw_ratio[ncol-5:ncol-2]= 0

  temp_r = trace_crude(saw_ratio#replicate(1, 2), nave = 1, thresh = 0.02)
  startr = transpose(temp_r[0, *])
  r_edges = n_elements(startr)
 
  temp_l = trace_crude(-1.0*saw_ratio#replicate(1, 2), nave = 1, thresh = 0.02)
  startl = transpose(temp_l[0, *])
  l_edges = n_elements(startl)

  if r_edges NE l_edges then begin
     print, 'Fix the arrays in startr, startl'
     plot, saw_ratio
     oplot, startr, saw_ratio[startr],ps=1
     oplot, startl, saw_ratio[startl],ps=1
     stop
  endif

  xr = trace_crude(saw, xstart=startr, yset=yr, xerr=rerr)
  xl = trace_crude(-1.0*saw, xstart=startl, yset=yl, xerr=lerr)

  mean_diff = djs_median(xr-xl,1)

  x = 0.5*(xr+xl)
  xy2traceset, yr, x, tset, ncoeff=3, yfit=xfit, $
         invvar=(rerr LT 900) * (lerr LT 900)

  lhedg = xfit - 0.5*mean_diff ## replicate(1,nrow) - 0.5
  rhedg = xfit + 0.5*mean_diff ## replicate(1,nrow) - 0.5

  return
end
