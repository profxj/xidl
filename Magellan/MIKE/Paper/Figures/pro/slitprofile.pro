pro slitprofile

;restore, 'rtwiflat_order52.dat'
;psfile = 'twiflat_order52.ps'
;restore, 'rtwiflat_order67.dat'
;psfile = 'twiflat_order67.ps'
;restore, 'rslitflat_order52.dat'
;psfile = 'slitflat_order52.ps'
restore, 'data/rslitflat_order67.dat'
psfile = 'slitflat_order67.ps'
;restore, 'bslitflat_order82.dat'
;psfile = 'slitflat_order82.ps'
;restore, 'bslitflat_order95.dat'
;psfile = 'slitflat_order95.ps'

lower = where(ywave LT 600.)
tlower = x_mkslitqa(slit_frac[lower], profile[lower], mask, $
                 ordr_str[q], per05=0.25, per95=0.75)

mid = where(ywave GT 600. AND ywave LT 1448.)
tmid = x_mkslitqa(slit_frac[mid], profile[mid], mask, $
                 ordr_str[q], per05=0.25, per95=0.75)

upper = where(ywave GT 1448.)
tupper = x_mkslitqa(slit_frac[upper], profile[upper], mask, $
                 ordr_str[q], per05=0.25, per95=0.75)

ss = poly(profile_x, ab)

x_psopen, psfile, /color, /square, /portrait
djs_plot, tlower.cen, tlower.median, yr=[0.95,1.04],ps=10, xr=[-1,1], $
       xtitle='Fractional slit position ', /nodata, $
       ytitle='Normalized flat-field profile', chars=1.5
;djs_oplot, tmid.cen, tmid.edgel * 0.7 + tmid.mid*0.3 , color='blue'
;djs_oplot, !x.crange, [0.98,0.98], lines=2
djs_oplot, tmid.cen, tmid.median, ps=10
djs_oplot, tmid.cen, tmid.mid
djs_oplot, !x.crange, [1.00,1.00], lines=2
;djs_oplot, tupper.cen, tupper.median+0.02, ps=10
;djs_oplot, tupper.cen, tmid.edger * 0.7 + tmid.mid*0.3, color='red'
djs_oplot, tupper.cen, tmid.edger, color='red'
djs_oplot, tmid.cen, tmid.edgel , color='blue'
;djs_oplot, !x.crange, [1.02,1.02], lines=2
djs_oplot, profile_x, ss, lines=3

x_psclose
end
