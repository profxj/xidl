pro figflux_1
;     sd = 'r'
     hcopy = 0

     set_plot, 'x'
     !p.font = 0
     device, decomposed=0

     A = FINDGEN(17) * (!PI*2/16.)
     usersym, COS(A), SIN(A), thick = 2.0

     !x.margin=[10,10]
     !y.margin=[4,2]

     for m = 0L,4 do begin 
         if (m gt 0) then stop 

         if (m eq 0) then name = 'GD71'          ;  terrible for b
         if (m eq 1) then name = 'EG274' ; not great in overlap. good in balmer lines, bad b
         if (m eq 2) then name = 'Feige110' ; okay for b
         if (m eq 3) then name = 'HR1996'   ; 5e5 flux!  looks good.    okay for b
         if (m eq 4) then name = 'EG21'   ; good everywhere ,good b  and big balmer lines!
         ;if (m eq 5) then name = 'HE0310m6060'   ;not enough flux for b

     for q = 1L, 2 do begin 
         if q then begin 
             sd = 'r'
             xrng = [4700,9000]
             nsig = 1.7  
             nsigarr = findgen(130) * 0. + nsig
             badord = [3,6,9,11]
             nsigarr[badord]= 10.

             if (hcopy gt 0.5) then begin ; this loop outside q test if one plot per side
                 !p.charsize = 1.5
                 !p.thick = 6
                 !x.thick = 6
                 !y.thick = 6
              ;  x_psopen , name+sd+'.ps', /maxs  ; one per side
                 x_psopen , name+'.ps', /maxs 
                 colors = GetColor(/Load, Start=1)
             endif else begin
                 colors = GetColor(/Load, Start=1)
                 !p.background = colors.white
                 !p.color = colors.black
                 !p.charthick = 3
                 !p.thick = 3
                 !x.thick = 3
                 !y.thick = 3
             endelse

         endif  else begin 
             sd = 'b'
             xrng = [3200,4800]
             nsig = 1.7
             nsigarr = findgen(130) * 0. + nsig
         endelse

         mike = mike_ar('$MIKE_PAP/strct.fits')

         star = xmrdfits('$MIKE_PAP/FSpec/'+name+'a_'+sd+'.fits.gz',1, /silent)

         ord = where (star.wave[10,*] gt 0 )
         gd = lindgen(130)
         for i = 0, n_elements(ord)-1 do begin 
             n = where(star.wave[*,ord[i]] gt 0.  ,numn)
             gd[i]= numn-5
             if (n[0] gt 0) then print, 'Pixel array is not starting at 0.'
         endfor

         HDR = 'sdf'
         err = mrdfits('$MIKE_PAP/FSpec/'+name+'a_'+sd+'_E.fits',0,HDR, /silent)
         flux = mrdfits('$MIKE_PAP/FSpec/'+name+'a_'+sd+'_F.fits',0,HDR, /silent)
         print, name, ':  ARIMASS,  HA =  ', sxpar(HDR,'AIRMASS'), '  ', sxpar(HDR,'HA')
         cdelt = sxpar(HDR,'CDELT1')
         crpix = sxpar(HDR,'CRPIX')
         crval = sxpar(HDR,'CRVAL1')
         wave = 10.^(crval + findgen(n_elements(flux)-1)*cdelt)

         if (q eq 1 and hcopy lt 1.5 and sd eq 'r')  then !p.multi=[0,1,2] & $
           plot, [0],xstyle=4, ystyle = 4
         if (q eq 1 and hcopy lt 1.5)  then !p.multi = [1,1,2]
         if (q eq 2 and hcopy lt 1.5)  then !p.multi = [2,1,2]
         plot, [0,1],[0,1] , color= colors.black $
               , yrange=[min(err),max(err)*5 < avg(err)*10], ystyle =5 $
               , xrange= xrng, xstyle = 1   $
               , title= name+':  ARIMASS= '+string(sxpar(HDR,'AIRMASS'))+ $
               ' HA: '+string(sxpar(HDR,'HA'))
         axis, yaxis = 1, yrange=[min(err),max(err)*5 < avg(err)*10] $
               , ytitle = 'error array', color = colors.black, /noerase
         oplot, wave, err, color= colors.darkgray

         ii=0L
         for i = ii, n_elements(ord)-1 do begin ;  1L
             temp = smooth (sqrt(star.var[*,ord[i] ]), 15, /edge_trunc)
             bttr = where (sqrt(star.var [0:gd[i],ord[i]]) lt $
                           min(sqrt(star.var [0:gd[i],ord[i]]) * nsigarr[i]))
             if i then clr = colors.red 
             oplot, star.wave[ bttr,ord[i] ], temp[bttr], color=clr
             clr = colors.black
         endfor

         if (q eq 1 and hcopy lt 1.5)  then !p.multi = [1,1,2]
         if (q eq 2 and hcopy lt 1.5)  then !p.multi = [2,1,2]
         plot, [0,1], [0,1], color= colors.black, /noerase $
               , yrange=[min(flux)>0,max(flux)< avg(flux)*3], ystyle = 9 $
               , xrange= xrng, xstyle = 1 $  
               , ytitle = 'flux'
         oplot, wave, flux, color= colors.gray

         ii=0L
         clr = colors.black
         for i = ii, n_elements(ord)-1 do begin ;  1L
             temp = smooth(star.fx[*,  ord[i] ], 15, /edge_trunc)
             if i then clr = colors.red 
             bttr = where (sqrt(star.var [0:gd[i],ord[i]]) lt $
                           min(sqrt(star.var [0:gd[i],ord[i]]) * nsigarr[i])) 
             oplot, star.wave[ bttr , ord[i] ] ,temp[bttr] , color= clr
             clr = colors.black
         endfor

         ;  if (hcopy gt 0.5) then x_psclose  ; if one plot per side

     endfor

     if (hcopy gt 0.5) then x_psclose

 endfor

 end
