pro fig_eg21, FEIGE=feige, YMNX=ymnx, REDUCE=reduce
;     sd = 'r'

 
  if not keyword_set(FEIGE) then begin
      psfil='fig_eg21.ps' 
      if keyword_set(REDUCE) then psfil='fig_eg21_rdx.ps'
      name = 'EG21'             ; good everywhere ,good b  and big balmer lines!
      readcol, 'data/feg21.dat', eg21_wav, eg21_flx
  endif else begin
      psfil='fig_feige110.ps'
      if keyword_set(REDUCE) then psfil='fig_feige110_rdx.ps'
      name = 'Feige110'             ; good everywhere ,good b  and big balmer lines!
      readcol, getenv('XIDL_DIR')+'Spec/Flux/HSTFIL/ffeige110.dat', eg21_wav, eg21_flx
  endelse
  

  x_psopen, psfil, /maxs
  clr = getcolor(/load)
  !p.multi=[0,1,2] 

  xmrg = [8,4]
  ymrg = [4,1]
  csz = 1.5
  lsz = 1.5

  for q = 1L, 2 do begin 
      if q EQ 2 then begin 
          sd = 'r'
          xrng = [4700,9000]
          nsig = 1.7  
          nsigarr = findgen(130) * 0. + nsig
          badord = [3,6,9,11]
          nsigarr[badord]= 10.
          sclwv = 6000.
          egoff = 10.
          yoff = 0.
      endif  else begin 
          sd = 'b'
          xrng = [3200,4800]
          nsig = 1.7
          nsigarr = findgen(130) * 0. + nsig
          sclwv = 4600.
          egoff = 20.
          yoff = 20.
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
      flux = x_readspec('$MIKE_PAP/FSpec/'+name+'a_'+sd+'_F.fits',wav=wave,sig=err,$
                        /AUTO, HEAD=hdr)
      ;; Scale
      flux = flux / 10.
      err = err / 10.
      star.fx = star.fx / 10.
      star.var = star.var / 100.

      print, name, ':  ARIMASS,  HA =  ', sxpar(HDR,'AIRMASS'), '  ', sxpar(HDR,'HA')
;      cdelt = sxpar(HDR,'CDELT1')
;      crpix = sxpar(HDR,'CRPIX')
;      crval = sxpar(HDR,'CRVAL1')
;      wave = 10.^(crval + findgen(n_elements(flux)-1)*cdelt)

;      if (q eq 1 and hcopy lt 1.5 and sd eq 'r')  then !p.multi=[0,1,2] & $
;        plot, [0],xstyle=4, ystyle = 4
;      if (q eq 1 and hcopy lt 1.5)  then !p.multi = [1,1,2]
;      if (q eq 2 and hcopy lt 1.5)  then !p.multi = [2,1,2]

      ;; Star
      plot, [0,1], [0,1], color= clr.black $
            , yrange=[min(flux)>0,max(flux)< avg(flux)*3 + egoff + yoff], $
            ystyle = 9, xrange= xrng, xstyle = 1, xmarg=xmrg $  
            , ytitle = 'Relative flux (erg s!u-1!N A!u-1!N cm!u-2!N)', $
            xtitle='Wavelength (Ang)', charsiz=csz, ymargi=ymrg
      print, 'xrng', xrng
      npix = n_elements(wave)
      if keyword_set(REDUCE) then $
        oplot, wave[0:npix-1:4], flux[0:npix-1:4], color= clr.gray $
      else $
        oplot, wave, flux, color= clr.gray
      
      ii=0L
      iclr = clr.black
      for i = ii, n_elements(ord)-1 do begin ;  1L
          temp = smooth(star.fx[*,  ord[i] ], 15, /edge_trunc)
          if i then iclr = clr.red 
          bttr = where (sqrt(star.var [0:gd[i],ord[i]]) lt $
                        min(sqrt(star.var [0:gd[i],ord[i]]) * nsigarr[i])) 
          oplot, star.wave[ bttr , ord[i] ] ,temp[bttr] , color= iclr
          iclr = clr.black
      endfor

      ;; 'Right answer'
      mn = min(abs(eg21_wav-sclwv),imn)
      px = where(abs(wave-sclwv) LT 20.)
      scl = median(flux[px]) / eg21_flx[imn]

      oplot, eg21_wav, eg21_flx*scl + egoff, color=clr.blue, linesty=1

      ;; Error
      plot, [0,1],[0,1] , color= clr.black, ymargi=ymrg $
            , yrange=[min(err),max(err)*5 < avg(err)*10], ystyle =5 $
            , xrange= xrng, xstyle = 5, xmarg=xmrg, /nodat, /noerase, charsiz=csz
;            , title= name+':  AIRMASS= '+string(sxpar(HDR,'AIRMASS'))+ $
;            ' HA: '+string(sxpar(HDR,'HA')), 
      axis, yaxis = 1, yrange=[min(err),max(err)*5 < avg(err)*10] $
            , ytitle = 'Relative Error', color = clr.black, /noerase, charsiz=csz

      if not keyword_set(REDUCE) then begin
          oplot, wave, err, color= clr.darkgray
      
          ii=0L
          for i = ii, n_elements(ord)-1 do begin ;  1L
              temp = smooth (sqrt(star.var[*,ord[i] ]), 15, /edge_trunc)
              bttr = where (sqrt(star.var [0:gd[i],ord[i]]) lt $
                            min(sqrt(star.var [0:gd[i],ord[i]]) * nsigarr[i]))
              if i then iclr = clr.red 
              oplot, star.wave[ bttr,ord[i] ], temp[bttr], color=iclr
              iclr = clr.black
          endfor
      endif

      ;; Label
      if q EQ 1 then begin
          if not keyword_set(YMNX) then ymnx = [min(err),max(err)*5 < avg(err)*10]
          ;; Name
          x1 = xrng[0] + 0.8*(xrng[1]-xrng[0])
          x2 = xrng[0] + 0.83*(xrng[1]-xrng[0])
          xyouts, x1, 0.85*ymnx[1], name, color=clr.black, charsiz=lsz
          ;; MIKE
          oplot, [x1,x2], replicate(0.80*ymnx[1],2), color=clr.black
          oplot, [x1,x2], replicate(0.79*ymnx[1],2), color=clr.red
          oplot, [x1,x2], replicate(0.78*ymnx[1],2), color=clr.darkgray
          xyouts, x2+10., 0.77*ymnx[1], 'MIKE', color=clr.black, charsiz=lsz
;          xyouts, x2, ymnx[1]*0.85, 'MIKE', color=clr.black, charsiz=lsz

          ;; Calib
          oplot, [x1,x2], replicate(0.72*ymnx[1],2), color=clr.blue, linesty=1
          xyouts, x2+10., 0.70*ymnx[1], 'Calib', color=clr.black, charsiz=lsz
      endif
          
      
  endfor
  x_psclose
  !p.multi=[0,1,1]

 end
