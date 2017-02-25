; from mike_flux 

     window, 0 
     set_plot, 'x'
     device, decomposed=0

;     hcopy = 0
     if (hcopy eq 1) then x_psopen, 'r_nrmlzd_blzfunc.ps', /maxs, /landscape
     yrng = [0.,1.1]
     yttl= 'Normalized Sensitivity'

     colors = GetColor(/Load, Start=1)
     !p.background = colors.white
     !p.charsize = 1.5

     mike = mike_ar('$MIKE_PAP/strct.fits')
                                ; mike_calibstd, mike, 2, ESOFIL='feg274.dat'   
                                ; mike_calibstd, mike, 3, ESOFIL='ffeige110.dat'
                                ; mike_calibstd, mike, 10, HSTFIL='fhr1996.dat'
     
;; red
     restore, '$MIKE_PAP/Extract/Sens_mr0020.idl'
     objstr = xmrdfits('$MIKE_PAP/Extract/Obj_mr0006.fits' $
                       , 1, STRUCTYP='mikeobjstrct', /silent)
     ordrstr = xmrdfits('$MIKE_PAP/Flats/OSTr_R_01.fits',1)
     idx = where(mike.frame eq '0006' and mike.side eq 2)
     
     full_rtio_arr = fltarr(n_elements(objstr),2050) ; objstr[*].npix)
     gpx_arr = intarr(n_elements(objstr),2050) 
     peak  = findgen(n_elements(objstr))*0.  ; for peak value and pixel #
     peakn  = bindgen(n_elements(objstr))*0  ; for peak value and pixel #
     
     ;; Legendre polynomials
     for jj=0L,n_elements(objstr)-1 do begin
         ;; Match up
         mtch = where(objstr[jj].order EQ ordr_fit,nmt)
         if nmt EQ 0 then begin
             print, 'mike_flux: No match to order = ', objstr[jj].order
             continue
         endif
         mtch = mtch[0]
         gpx   = lindgen(objstr[jj].npix)
         
         full_rtio = x_calcfit(objstr[jj].wave[gpx], fitstr=tot_fit[mtch])
         objstr[jj].flux[gpx] = objstr[jj].fx[gpx]/full_rtio/mike[idx].exp
         b = where(objstr[jj].var[gpx] GT 0.)
         objstr[jj].sig[gpx[b]] = sqrt(objstr[jj].var[gpx[b]]) $
                                 / full_rtio[b] / mike[idx].exp
         objstr[jj].flg_flux = 1
         full_rtio_arr[jj,0:objstr[jj].npix-1] =full_rtio
         
         gpx_arr[jj,0:n_elements(gpx)-1] = gpx
         peak[jj] = max(full_rtio, n)
         peakn[jj] = n
     endfor
     
     jj = 10
     n = peakn[jj]
     nlo = peakn[jj+1]
     nhi = peakn[jj-1]
     centw =     objstr[jj].wave[gpx_arr[jj,n]]
     low = abs( objstr[jj+1].wave[gpx_arr[jj+1,nlo]] + centw ) * 0.5
     hiw = abs( objstr[jj-1].wave[gpx_arr[jj-1,nhi]] + centw ) * 0.5
     rel_wav = (objstr[jj].wave[gpx_arr[jj,0:objstr[jj].npix-1]] - centw ) $
               / ((hiw-low) )      ;   objstr[jj].order/51. 
     plot, rel_wav $
           , full_rtio_arr[jj,gpx_arr[jj,0:objstr[jj].npix-1]]  /peak[jj] $
           , color= colors.black, /noerase $
           , xrange = [-1.2, 1.2], xstyle=1 $
           , yrange = yrng,  ystyle=1 $
           , xtitle = 'Relative position in order in units of FSR', ytitle = yttl
     
     norm = objstr[jj].wave[gpx_arr[jj,peakn[jj]]] * objstr[jj].order 

     i = 10
     for jj=2L,n_elements(objstr)-2 do begin
         i = i+2
         n = peakn[jj]
         nlo = peakn[jj+1]
         nhi = peakn[jj-1]
         centw =     objstr[jj].wave[gpx_arr[jj,n]]
         low = abs( objstr[jj+1].wave[gpx_arr[jj+1,nlo]] + centw ) * 0.5
         hiw = abs( objstr[jj-1].wave[gpx_arr[jj-1,nhi]] + centw ) * 0.5
         rel_wav = (objstr[jj].wave[gpx_arr[jj,0:objstr[jj].npix-1]] - centw ) $
                   / ((hiw-low))  ;   objstr[jj].order/51. 
         oplot, rel_wav $
               , full_rtio_arr[jj,gpx_arr[jj,0:objstr[jj].npix-1]]  /peak[jj] $
               , color= i+10
     endfor
     
;  for jj= 0L, n_elements(objstr)-1 do  $
;   print, jj, peak[jj] $
;          , objstr[jj].order $
;          , objstr[jj].wave[gpx_arr[jj,n]] $
;          , objstr[jj].wave[gpx_arr[jj,n]] * objstr[jj].order $
;          , objstr[jj].wave[gpx_arr[jj,n]] * objstr[jj].order /norm $
;          , full_rtio_arr[jj,gpx_arr[jj,n]]/peak[jj]

     jj = 10
     n = peakn[jj]
     nlo = peakn[jj+1]
     nhi = peakn[jj-1]
     centw =     objstr[jj].wave[gpx_arr[jj,n]]
     low = abs( objstr[jj+1].wave[gpx_arr[jj+1,nlo]] + centw ) * 0.5
     hiw = abs( objstr[jj-1].wave[gpx_arr[jj-1,nhi]] + centw ) * 0.5
     rel_wav = (objstr[jj].wave[gpx_arr[jj,0:objstr[jj].npix-1]] - centw ) $
               / ((hiw-low))      ;   objstr[jj].order/51. 
     plot, rel_wav $
           , full_rtio_arr[jj,gpx_arr[jj,0:objstr[jj].npix-1]]  /peak[jj] $
           , color= colors.black, /noerase $
           , xrange = [-1.2, 1.2], xstyle=1 $
           , yrange = yrng,  ystyle=1, thick=2
     xyouts, [-1.05], [1.0], 'Red side', color = colors.black

     if (hcopy eq 1) then x_psclose

     
end



