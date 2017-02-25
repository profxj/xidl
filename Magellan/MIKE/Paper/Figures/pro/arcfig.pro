
   colors = GetColor(/Load, Start=1)
   !p.multi=[0,1,1]
   !p.background = colors.white
   !p.charsize = 1.5
   x_psopen, 'arcimg.ps'
   loadct, 0
   arc = readfits ('data/r0033.fits') 
   arc[821,*]=0
   tvscl, 2000. - ((arc[416:966,838:1388]<2000.) > 1100)
   x_psclose

   arc = readfits ('data/r0033.fits')
   !p.multi=[0,1,1]
   !p.background = colors.white
   x_psopen, 'arcplot.ps'
   y = arc[821,838:1388]
   plot, findgen(1388-838+1)+838, y, color=colors.black $
         , xtitle ='Row number (y)', ytitle = 'DN/pixel' $
         , ystyle=1, yrange=[500,5000] $
         , xstyle=1, xrange=[838, 1388]
   x_psclose

end
