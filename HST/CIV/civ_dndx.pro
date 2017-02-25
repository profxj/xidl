;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; civ_dndx.pro               
; Author: Kathy Cooksey                      Date: 28 Jan 2009
; Project: 
; Description: 
; Input: 
; Optional:
; Output: 
; Example:
; History:
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@civ_sensitivity
@civ_calcewn

pro civ_dndx_plot,strct_fil, sens_fil, psfil, siiv=siiv, ncolm=ncolm, $
                  csize=csize,psize=psize,lthick=lthick,$
                  xrng=xrng,yrng=yrng, fit_fil=fit_fil, label=label
  
  if size(strct_fil,/type) eq 8 then strct = strct_fil $
  else strct = xmrdfits(strct_fil,1,/silent)
  nciv = n_elements(strct)


  ;; Params
  ge_sign = '>'
  if not keyword_set(csize) then csize = 2.
  if not keyword_set(psize) then psize = 2.
  if not keyword_set(lthick) then lthick = 5
  if not keyword_set(yrng) then yrng = [0.,15]

  if keyword_set(siiv) then begin
     if keyword_set(ncolm) then begin
        ytitle = 'd!8N!X/d!8X!X (!8N!X(Si!E+3!N) '+ge_sign+' !8N!X!Dlim!N)'
        xtitle = 'log !8N!X!Dlim!N'
     endif else begin
        ytitle = 'd!8N!X/d!8X!X (!8W!X!Dr,1393!N '+ge_sign+' !8W!X!Dlim!N)'
        xtitle = '!8W!X!Dlim!N (m'+STRING("305B) +')' 
     endelse 
  endif else begin
     if keyword_set(ncolm) then begin
        ytitle = 'd!8N!X/d!8X!X (!8N!X(C!E+3!N) '+ge_sign+' !8N!X!Dlim!N)'
        xtitle = 'log !8N!X!Dlim!N'
     endif else begin
        ytitle = 'd!8N!X/d!8X!X (!8W!X!Dr,1548!N '+ge_sign+' !8W!X!Dlim!N)'
        xtitle = '!8W!X!Dlim!N (m'+STRING("305B) +')' 
     endelse 
  endelse 

  if keyword_set(ncolm) then begin
     civ_dndx, rslt, sens_fil, strct, siiv=siiv, /ncolm, data=data
     if not keyword_set(xrng) then xrng = [13.,15.5]
     xloc = 0.75
  endif else begin
     civ_dndx, rslt, sens_fil, strct, siiv=siiv, /ew, data=data
     if not keyword_set(xrng) then xrng = [10.,100.]
     xloc = 0.8
  endelse 

  ;; Set correct histogram array
  xplt = fltarr(2*nciv+1)
  yplt = fltarr(2*nciv+1)
  for jj=0,nciv-2 do begin
     xplt[2*jj] = data[jj]
     xplt[2*jj+1] = data[jj+1]
     yplt[2*jj] = rslt[jj,0]
     yplt[2*jj+1] = rslt[jj,0]
  endfor 
  xplt[2*nciv] = data[nciv-1]
  yplt[2*nciv] = 0

  x_psopen,psfil,/maxs
  !p.multi = [1,1,1]
  !x.margin = [8,3]
  !y.margin = [4,2]

  clr = getcolor(/load)

  plot,xrng,yrng,/nodata,/ystyle,/xstyle,background=clr.white,color=clr.black,$
       ytitle=ytitle,xtitle=xtitle,charsize=csize

  
  ;; Overplot errors 
  oploterror,data,rslt[*,0],rslt[*,1],/lobar,errcolor=pclr,thick=lthick,psym=3
  oploterror,data,rslt[*,0],rslt[*,2],/hibar,errcolor=pclr,thick=lthick,psym=3

  ;; dN/dX
  oplot,xplt,yplt,color=clr.black,psym=10.,thick=lthick

  ;; Fit
  if keyword_set(fit_fil) then begin
     if size(fit_fil,/type) eq 7 then $
        fitstrct = xmrdfits(fit_fil,1,/silent) $
     else fitstrct = fit_fil

     ;; Linear space first 
     yfit = dblarr(1000L)
     if keyword_set(ncolm) then begin
        xfit = 10^xrng[0] + dindgen(1000L) * (10^xrng[1]-10^xrng[0])/1000L
        civ_dndx, yfit, sens_fil, siiv=siiv, ncolm=xfit, fit_fil=fitstrct
        oplot,alog10(xfit),yfit[*,0], linestyle=0,color=clr.red,thick=lthick 
     endif else begin 
        xfit = xrng[0] + dindgen(1000L) * (xrng[1]-xrng[0])/1000L
        civ_dndx, yfit, sens_fil, siiv=siiv, ew=xfit, fit_fil=fitstrct
        oplot,xfit,yfit[*,0], linestyle=0,color=clr.red,thick=lthick 
     endelse 

     ;; Legend (upper right)
     if keyword_set(label) then begin
        dx = xrng[1] - xrng[0]
        dy = yrng[1] - yrng[0]
        xyouts, xloc*dx + xrng[0], 0.93*dy + yrng[0], $
                '!9a!X = '+strtrim(string(fitstrct.alpha,format='(f5.2)'),2),$
                charsize=csize, color=clr.black
        xyouts, xloc*dx + xrng[0], 0.88*dy + yrng[0], $
                '!8k!X = '+strtrim(string(fitstrct.coeff,format='(e9.2)'),2),$
                charsize=csize, color=clr.black 
        xyouts, xloc*dx + xrng[0], 0.83*dy + yrng[0], $
                'P!DKS!N = '+strtrim(string(fitstrct.prob_ks,format='(f5.3)'),2),$
                charsize=csize, color=clr.black

     endif                      ; /label
  endif                         ; fit_fil=

  x_psclose
  print,'civ_dndx: created ',psfil

end                             ; civ_dndx_plot


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro civ_dndx, rslt, sens_fil, strct_fil, siiv=siiv, data=data, xdata=xdata,$
              sigdata=sigdata, ew=ew, ncolm=ncolm, fit_fil=fit_fil,dndz=dndz, $
              sigcoeff_dndx=sigcoeff_dndx, sigalpha_dndx=sigalpha_dndx,$
              intlim_max=intlim_max,_extra=extra

  if keyword_set(siiv) then dblt_name = 'SiIV' $
  else dblt_name = 'CIV'

  if keyword_set(ncolm) then begin
     nbin = n_elements(ncolm) 
     if nbin eq 1 and ncolm[0] eq 1 then ret_obs = 1
  endif else begin
     if not keyword_set(ew) then begin
        ew = 1 
        nbin = 1
     endif else nbin = n_elements(ew) 
     if nbin eq 1 and ew[0] eq 1 then ret_obs = 1
  endelse 

  if keyword_set(fit_fil) then begin
     ;;;;;;;;;
     ;; If fit_fil provided, return fit value
     ;;;;;;;;;
     if size(fit_fil,/type) eq 7 then begin
        fitstrct = xmrdfits(fit_fil,1,/silent) 
        errstrct = xmrdfits(fit_fil,3,/silent) 
     endif else begin
;        stop,'civ_dndx: if want fit dN/dX errors, input file' 
        fitstrct = fit_fil
     endelse 

     if keyword_set(ncolm) then data = ncolm $
     else data = ew

     rslt = dblarr(nbin,3)
     for ii = 0,nbin-1 do begin
        alpha1 = 1. + fitstrct.alpha
        
        C =  fitstrct.coeff/(alpha1 * fitstrct.datanorm^fitstrct.alpha)
        if keyword_set(intlim_max) then begin
           if intlim_max eq 0. then $
              rslt[ii,0] =  -C * data[ii]^alpha1 $ ; intlim_max = inf
           else begin
              if intlim_max lt 0 then $
                 rslt[ii,0] =  C * (fitstrct.intlim[2]^alpha1 - data[ii]^alpha1) $
              else $
                 rslt[ii,0] =  C * (intlim_max^alpha1 - data[ii]^alpha1) 
           endelse 
        endif else rslt[ii,0] =  -C * data[ii]^alpha1 ; intlim_max = inf
        
        if keyword_set(errstrct) then begin
           ;; This is right
           dndx = dblarr(n_elements(errstrct.surf)) ; minimize computation
           alpha_grid = dndx
           coeff_grid = dndx
           tmp = fitstrct       ; be safe
           logL_grid = reform(errstrct.logL,n_elements(errstrct.coeff_grid),$
                              n_elements(errstrct.alpha_grid)) ; for reference
           
           ;; Be consistent in error analysis
           if keyword_set(intlim_max) then intmax = intlim_max $
           else begin
              ;; Check quality of error analysis
              if fitstrct.alpha+fitstrct.sigalpha[1] ge -1. then begin
                 ;; This will be slightly inconsistent with the
                 ;; calculated value, which has assumed infinite upper limit
                 intmax = fitstrct.intlim[2] 
                 print,'civ_dndx: shallow slope; use upper integration limit',$
                       intmax
              endif else intmax = 0
           endelse 

           for nn=0L,n_elements(errstrct.surf)-1 do begin 
              mm = array_indices(logL_grid,errstrct.surf[nn])
              alpha_grid[nn] = errstrct.alpha_grid[mm[1]] 
              coeff_grid[nn] = errstrct.coeff_grid[mm[0]] 
              tmp.coeff = coeff_grid[nn]
              tmp.alpha = alpha_grid[nn]

              if keyword_set(ncolm) then $
                 civ_dndx, tmprslt, siiv=siiv, ncolm=data[ii], fit_fil=tmp, $
                           intlim_max=intmax $    ; won't call this section
              else civ_dndx, tmprslt, siiv=siiv, ew=data[ii], fit_fil=tmp, $
                             intlim_max=intmax

              dndx[nn] = tmprslt[0] 
           endfor

           ;; Find extrema of error ellipse values
           dndx_mn = min(dndx,imn,max=dndx_mx,subscript_max=imx) 

           rslt[ii,1] = rslt[ii,0] - dndx_mn
           rslt[ii,2] = dndx_mx - rslt[ii,0]

           ;; Return locations of these limits
           sigcoeff_dndx = dblarr(2)
           sigalpha_dndx = dblarr(2)
           mm = array_indices(logL_grid,errstrct.surf[imn])
           sigcoeff_dndx[0] = fitstrct.coeff - errstrct.coeff_grid[mm[0]]
           sigalpha_dndx[0] = fitstrct.alpha - errstrct.alpha_grid[mm[1]]
           mm = array_indices(logL_grid,errstrct.surf[imx])
           sigcoeff_dndx[1] = errstrct.coeff_grid[mm[0]] - fitstrct.coeff
           sigalpha_dndx[1] = errstrct.alpha_grid[mm[1]] - fitstrct.alpha
        endif else begin
           ;; This is wrong
;           stop,'civ_dndx: do not compute fit dN/dX error this way'
           ;; Lower error
           alpha1 = 1. + fitstrct.alpha - fitstrct.sigalpha[0]
           rslt[ii,1] =  rslt[ii,0] + (fitstrct.coeff + $
                                       fitstrct.sigcoeff[1])/alpha1 * $
                         data[ii]^alpha1 / fitstrct.datanorm^fitstrct.alpha
           ;; Upper error
           alpha1 = 1. + fitstrct.alpha + fitstrct.sigalpha[1]
           rslt[ii,2] =  -(fitstrct.coeff - fitstrct.sigcoeff[0])/alpha1 * $
                         data[ii]^alpha1 / fitstrct.datanorm^fitstrct.alpha - $
                         rslt[ii,0]
        endelse 
     endfor                     ; loop nbin
  endif else begin              ; fit_fil=
     ;;;;;;;;;
     ;; Observed results
     ;;;;;;;;;
     ;; Read in
     if not keyword_set(xdata) then begin 
        ;; Don't overwrite if input
        if size(sens_fil,/type) eq 7 then sens = xmrdfits(sens_fil,1,/silent) $
        else sens = sens_fil          
     endif 
     if not keyword_set(data) then begin
        civ_group, strct, strct_fil, _extra=extra ; sample selection
        ion_indx = where(strtrim(strct[0].ion,2) ne '',dblt_flag)
        ion_indx = ion_indx[0]

        ;; Set up observed data
        if keyword_set(ncolm) then begin
           data = civ_calcewn_ndblt(strct,dblt_name,/log,signcolm=sigdata) 
           if keyword_set(sens) then $
              xdata = civ_sensitivity_x(sens,ncolm=data,z=zdata,$
                                        signcolm=sigdata,sigx=sigx,sigz=sigz)
        endif else begin
           data = strct.ew[ion_indx] 
           sigdata = strct.sigew[ion_indx] 
           if keyword_set(sens) then $
              xdata = civ_sensitivity_x(sens,ew=data,z=zdata,$
                                        sigew=sigdata,sigx=sigx,sigz=sigz)
        endelse 
        
        srt = sort(data)        ; order
        data = data[srt]
        sigdata = sigdata[srt]
        strct = strct[srt]
        xdata = xdata[srt]
        if keyword_set(dndz) then begin
           xdata = zdata[srt]
           sigx = sigz
        endif 
        sigx = sigx[srt,*]
     endif 
     if not keyword_set(xdata) then begin 
        ;; Create if does not exist
        if keyword_set(ncolm) then begin
           xdata = civ_sensitivity_x(sens,ncolm=data,z=zdata,$
                                     signcolm=sigdata,sigx=sigx,sigz=sigz)
        endif else begin
           xdata = civ_sensitivity_x(sens,ew=data,z=zdata,$
                                     sigew=sigdata,sigx=sigx,sigz=sigz)
        endelse      
        srt = sort(data)        ; order
        data = data[srt]
        sigdata = sigdata[srt]
        xdata = xdata[srt]
        if keyword_set(dndz) then begin
           xdata = zdata[srt]
           sigx = sigz
        endif 
        sigx = sigx[srt,*]
     endif 

     nciv = n_elements(data)
     rslt = dblarr(nciv,3)
     sigdndx = dblarr(2)
     sigpoiss = dblarr(2)

     ;; Observed dN/dX
     for ii=0,nciv-1 do begin
        rslt[ii,0] = total(1./xdata[ii:nciv-1])
        
        ;; Errors: component depending on X pathlength (which depends
        ;; on EW or ncolm)
        sigdndx[0] = sqrt(total((sigx[ii:nciv-1,0]/xdata[ii:nciv-1]^2)^2))
        sigdndx[1] = sqrt(total((sigx[ii:nciv-1,1]/xdata[ii:nciv-1]^2)^2))

        ;; Component depending on Poisson statistics
        p = x_poisscl(nciv-ii,0.683,/silent)
        sigpoiss[0] = ((nciv-ii)-p[1])
        sigpoiss[1] = (p[0]-(nciv-ii))
        
        ;; Add errors in quadrature: f(x,y) = x/y
        ;; (sig_f/f)^2 = (sig_x/x)^2 + (sig_y/y)^2
        ;; But know that errors anti-correlated
        rslt[ii,1] = rslt[ii,0]*sqrt((sigdndx[0]/rslt[ii,0])^2 + $
                                     (sigpoiss[0]/(nciv-ii))^2)
        rslt[ii,2] = rslt[ii,0]*sqrt((sigdndx[1]/rslt[ii,0])^2 + $
                                     (sigpoiss[1]/(nciv-ii))^2)
     endfor                     ; loop rslt 
     
     ;; Interpolate
     if not keyword_set(ret_obs) then begin
        obsrslt = rslt
        obsdata = data
        if keyword_set(ncolm) then data = ncolm[sort(ncolm)] $
        else data = ew[sort(ew)]
        rslt = dblarr(nbin,3)

        for ii=0,nbin-1 do begin
           if data[ii] lt obsdata[1] then begin
              rslt[ii,*] = obsrslt[0,*]
           endif else begin
              if data[ii] gt obsdata[nciv-1] then begin
                 rslt[ii,*] = obsrslt[nciv-1,*]
              endif else begin
                 rslt[ii,0] = interpol(obsrslt[*,0], obsdata, data[ii]) 
                 rslt[ii,1] = interpol(obsrslt[*,0]-obsrslt[*,1], obsdata, $
                                       data[ii]) 
                 rslt[ii,2] = interpol(obsrslt[*,0]+obsrslt[*,2], obsdata, $
                                       data[ii]) 
                 rslt[ii,1] = rslt[ii,0] - rslt[ii,1]
                 rslt[ii,2] = rslt[ii,2] - rslt[ii,0]
              endelse 
           endelse 
        endfor
     endif                      ; interpolate

  endelse                       ; observed

end
