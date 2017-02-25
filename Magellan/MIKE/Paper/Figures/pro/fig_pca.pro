pro fig_pca

restore, 'data/red_edgeflat.dat'
ordr = mrdfits('data/OStr_R_02sep04.fits',1)
;save, all_ordr, aset, gd, est_center, pca_out, eigenvec, filename='data/r.dat'

ycol = findgen(2048)
x0_all = poly(all_ordr, aset)
nordr = n_elements(all_ordr)
ngd = n_elements(gd)

x_psopen, 'fig_pca.ps', /maxs, /portrait
plot, [0,1],[0,1], /nodata, ys=5, xs=5
xyouts, -0.09, 0.55, 'Eigenfunction Response', $
       orient=90., align=0.5, /noclip, chars=1.5
xyouts, 0.97, 0.55, 'PCA coefficients (pixels)', $
       orient=-90., align=0.5, /noclip, chars=1.5

x1 = 0.3
x2 = 0.9
y1 = 0.80 & y2 = y1 + 0.14
;plot, ordr.order, x0_all, xr=[73, 38], /xs, yr=[-50, 1050], $
;       position=[x1,y1,x2,y2], ys=9, ychars=0.01, xchars=0.01, /noerase
;axis, yaxis=1, /noerase, yr=[-50, 1050], /ys, ychars=1.4
;oplot, ordr[gd].order, est_center, ps=1

y1 = 0.20 & y2 = 0.9
plot, ordr.order, poly(x0_all, pca_out.coeff0), position=[x1,y1,x2,y2], $
            xr=[73, 38], /xs, yr=[-28, 28],ys=9, /noerase, $
            xchars=1.4,  ychars=0.01, xtitle='Order Number'
axis, yaxis=1, /noerase, yr=!y.crange, /ys, ychars=1.4
oplot, ordr[gd].order, pca_out.hidden[0,*], ps=1
oplot, ordr[gd].order, pca_out.hidden[1,*]
oplot, ordr[gd].order, pca_out.hidden[2,*]
oplot, ordr[gd].order, pca_out.hidden[3,*]
oplot, ordr[gd].order, pca_out.hidden[4,*]
xyouts, replicate(71.5,5), pca_out.hidden[*,0] - 0.7, $
        ['1','2','3','4','5'], chars=1.4

y1 = 0.76 & y2 = y1 + 0.14
plot, ycol, eigenvec[*,0], position=[0.1,y1,x1,y2], $
            /noerase, /xs, xchars=0.01, yr=[-1.2,1.2], /ys, ychars=1.4
xyouts, [1800], [0.7], ['1'], chars=1.4

y1 = 0.62 & y2 = y1 + 0.14
;plot, ordr.order, poly(x0_all, pca_out.coeff1), position=[x1,y1,x2,y2], $
;            xr=[73, 38], /xs, yr=[-23.15, -22.95],ys=9, /noerase, $
;            xchars=0.01,  ychars=0.01
;axis, yaxis=1, /noerase, yr=!y.crange, /ys, ychars=1.4
;oplot, ordr[gd].order, pca_out.hidden[1,*], ps=1
plot, ycol, eigenvec[*,1], position=[0.1,y1,x1,y2], $
            /noerase, /xs, xchars=0.01, yr=[-1.2,1.2], /ys, ychars=1.4
xyouts, [1800], [0.7], ['2'], chars=1.4

y1 = 0.48 & y2 = y1 + 0.14
;plot, ordr.order, replicate(pca_out.high_fit[0,2],nordr), $
;            position=[x1,y1,x2,y2], $
;            xr=[73, 38], /xs, yr=[-3.22, -3.18],ys=9 , /noerase, $
;            xchars=0.01,  ychars=0.01
;axis, yaxis=1, /noerase, yr=!y.crange, /ys, ychars=1.4
;oplot, ordr[gd].order, pca_out.hidden[2,*], ps=1
plot, ycol, eigenvec[*,2], position=[0.1,y1,x1,y2], $
            /noerase, /xs, xchars=0.01, yr=[-1.2,1.2], /ys, ychars=1.4
xyouts, [1800], [0.7], ['3'], chars=1.4
y1 = 0.34 & y2 = y1 + 0.14
;plot, ordr.order, replicate(pca_out.high_fit[0,3],nordr), $
;            position=[x1,y1,x2,y2], $
;            xr=[73, 38], /xs, yr=[1.786, 1.804],ys=9 , /noerase, $
;            xchars=0.01,  ychars=0.01
;axis, yaxis=1, /noerase, yr=!y.crange, /ys, ychars=1.4
;oplot, ordr[gd].order, pca_out.hidden[3,*], ps=1
plot, ycol, eigenvec[*,3], position=[0.1,y1,x1,y2], $
            /noerase, /xs, xchars=0.01, yr=[-1.2,1.2], /ys, ychars=1.4
xyouts, [1800], [0.7], ['4'], chars=1.4
y1 = 0.20 & y2 = y1 + 0.14
;plot, ordr.order, replicate(pca_out.high_fit[0,4],nordr), $
;            position=[x1,y1,x2,y2], $
;            xr=[73, 38], /xs, yr=[-0.559, -0.536],ys=9 , /noerase, $
;            ychars=0.01, xtitle='Order Number', xchars=1.4
;axis, yaxis=1, /noerase, yr=!y.crange, /ys, ychars=1.4
;oplot, ordr[gd].order, pca_out.hidden[4,*], ps=1
plot, ycol, eigenvec[*,4], position=[0.1,y1,x1,y2], $
            /noerase, /xs, yr=[-1.2,1.2], /ys, $
            xtitle='Row Number (pixels)', xchars=1.4, ychars=1.4
xyouts, [1800], [0.7], ['5'], chars=1.4

x_psclose
end
