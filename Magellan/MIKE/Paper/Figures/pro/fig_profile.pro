pro fig_profile

if NOT keyword_set(profile_temp) then restore, 'data/profile_ordr47_iter1.dat'

good = where(profile_ivar GT 0)

print, min(slit_dist[good]), max(slit_dist[good])

slit_bin = findgen(251)*0.01 - 1.25
profile_med = fltarr(251)
profile_75 = fltarr(251)
profile_25 = fltarr(251)

for i=0,250 do begin

   in = where(slit_dist[good] GE slit_bin[i] - 0.05 AND $
              slit_dist[good] LT slit_bin[i] + 0.05, nin)
  
   if nin LT 10 then continue

   s = sort(profile_temp[good[in]])
   profile_med[i] = profile_temp[good[in[s[nin*0.5]]]] 
   profile_75[i] = profile_temp[good[in[s[nin*0.75]]]] 
   profile_25[i] = profile_temp[good[in[s[nin*0.25]]]] 
endfor

mean_profile = bspline_valu(slit_bin, profile_set, $
         x2=slit_bin*0 + 0.5*(profile_set.xmin + profile_set.xmax))
top_profile = bspline_valu(slit_bin, profile_set, $
         x2=slit_bin*0 + profile_set.xmax)
bot_profile = bspline_valu(slit_bin, profile_set, $
         x2=slit_bin*0 + 0.5*profile_set.xmin)

x_psopen, 'fig_profile.ps', /maxs, /portrait

plot, slit_bin, profile_med, xr=[-1,1], yr=[-0.5, 3.5],/ys,ps=10, $
   xtitle='Slit position relative to Object Centroid (HSL)', $
;$   xtitle='Slit Position relative to Object Trace (pixels)', $
   ytitle='Object Profile Cross-section', $
   chars=1.5

oplot, slit_bin, profile_med, ps=10
errplot, slit_bin, profile_25, profile_75, width=0.003
djs_oplot, slit_bin, mean_profile, color='green'
djs_oplot, slit_bin, top_profile, color='red'
djs_oplot, slit_bin, bot_profile, color='red'

djs_oplot, -1.0*[cutoff_aper[0], cutoff_aper[0]] + profile_cen, $
                       [0.0, 1.0], thick=8, color='blue', linesty=2
djs_oplot, [cutoff_aper[1], cutoff_aper[1]] + profile_cen, $
                       [0.0, 1.0], thick=8, color='blue', linesty=2

djs_oplot, [-1.0*base_aper[0], -1.0*base_aper[0], $
                   base_aper[1], base_aper[1]] + profile_cen, $
                       [2.0, 0.0, 0.0, 2.0], thick=8, color='burlywood', $
           linestyle=1

djs_oplot, [profile_lwhm, profile_lwhm, profile_rwhm, profile_rwhm], $
       [profile_max/2, profile_max/2.2, profile_max/2.2, profile_max/2], $
                  thick = 6

xyouts, [0.0], [profile_max/2], $
                [string(spatial_fwhm, '"', format='(f5.2,a)')], $
                 alignment=0.5, charsize=1.5

x_psclose

end
