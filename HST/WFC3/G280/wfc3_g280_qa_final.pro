;+ 
; NAME:
;  wfc3_g280_qa_final
;
; PURPOSE:
;   Creates a QA file of the reduction procedure
;
; CALLING SEQUENCE:
;   
;   wfc3_g280_qa_final, wfc3_g280_strct, ii, specim, IMGSZ=, $
;                       QADIR=, WAVERG=, STRETCH=, $
;                       SIGMA=
;
; INPUTS:
;
;   wfc3_g280_strct -- the wfc3_g280 structure
;   ii -- the index of the object in the structure
; 
; RETURNS:
;
; OUTPUTS:
;   Creates a QA file for the requested object
;
; OPTIONAL KEYWORDS:
;  IMGSZ= -- Size of the direct image
;  QADIR= -- Directory for the QA file
;  WAVERG= -- Wavelength Range to plot
;  STRETCH= -- Stretch for plotting the 2D spectra
;  SIMGA= -- Sigma used in extraction
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
; wfc3_g280_qa_final, wfc3_g280_strct, ii, SPECIM=specim, IMGSZ=imgsz, $
;                     QADIR=qadir, WAVERG=waverg, STRETCH=stretch, $
;                     SIGMA=tsigma
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   10-Jun-2016 Written by MN
;------------------------------------------------------------------------------

pro wfc3_g280_qa_final, wfc3_g280_strct, ii, SPECIM=specim, IMGSZ=imgsz, $
                        QADIR=qadir, WAVERG=waverg, STRETCH=stretch, $
                        SIGMA=sigma, BOXCAR=boxcar, BOX_SIZE=box_size
  
  if keyword_set(waverg) then wvrg=waverg else $
     wvrg=[1800d,6500d]
  if not keyword_set(stretch) then stretch=10
  if not keyword_set(imgsz) then imgsz=7L
  
  thick=5
  psfile=qadir+wfc3_g280_strct(ii).name+'.eps'
  ps_start, psfile, /encapsulated, /nomatch, xsize=9, ysize=7

  ext=7-3*wfc3_g280_strct(ii).chip
  specim=xmrdfits(wfc3_g280_strct(ii).spec_fil,ext)
  img=xmrdfits(wfc3_g280_strct(ii).img_fil,ext)
  ;; trim direct image
  
  x0=round(wfc3_g280_strct(ii).x0)
  y0=round(wfc3_g280_strct(ii).y0)
  img_trim=img[x0-imgsz:x0+imgsz,y0-imgsz:y0+imgsz]
  
  ;; loop over the beams
  for jj=0L, 1 do begin
  
     case jj of
        0: begin
           cnt=wfc3_g280_strct(ii).cnta
           trace_x=wfc3_g280_strct(ii).trace_xa
           trace_y_orig=wfc3_g280_strct(ii).trace_ya_orig
           trace_yfit=wfc3_g280_strct(ii).trace_ya_fit
           trace_sigma_fit=wfc3_g280_strct(ii).trace_sigma_fita
           wave=wfc3_g280_strct(ii).wavea
           badpix=wfc3_g280_strct(ii).badpixa
           bm='Beam A'
           clr='red3'
        end
        1: begin
           cnt=wfc3_g280_strct(ii).cntc
           trace_x=wfc3_g280_strct(ii).trace_xc
           trace_y_orig=wfc3_g280_strct(ii).trace_yc_orig
           trace_yfit=wfc3_g280_strct(ii).trace_yc_fit
           trace_sigma_fit=wfc3_g280_strct(ii).trace_sigma_fitc
           wave=wfc3_g280_strct(ii).wavec
           badpix=wfc3_g280_strct(ii).badpixc
           bm='Beam C'
           clr='blu3'
        end
        else: stop
     endcase
     
     ;; sigma
     tsigma=median(trace_sigma_fit(0:cnt-1))*1.<2.
     if keyword_set(SIGMA) then sigma=sigma else sigma=tsigma
     if keyword_set(boxcar) then sigma=box_size
     
     ;; Image
     dim=size(specim,/dim)
     tmpgd=where(trace_x gt 0 and trace_x lt dim(0) and $
                 wave gt wvrg[0] and wave lt wvrg[1])
     gd=min(tmpgd) + lindgen(max(tmpgd)-min(tmpgd))
     xmn = min(trace_x(gd), max=xmx)
     xmnx = round([(dim(0)-1)<xmn>0L,(dim(0)-1)<xmx>0L])
     ymn = min(trace_y_orig(gd), max=ymx)
     YBUFF = 30L
     ymn = (dim(1)-1)<(ymn-YBUFF)>0L 
     ymx = (dim(1)-1)<(ymx+YBUFF)>0L
     ymnx = round([ymn,ymx])
     img_cut = specim[xmnx[0]:xmnx[1], ymnx[0]:ymnx[1]]
     szc = size(img_cut,/dimen)

     ;; if beam A then flip the image
     if jj eq 0 then begin
        img_cut=reverse(img_cut)
        xmnx=reverse(xmnx)
     endif else xgd=gd
     
     ;; position of plot
     pos=[0.10,0.75-jj*0.25,0.95,0.95-jj*0.25]

     ;; plot the image
     if jj eq 1 then xtit='Column' else xtit=''
     cgimage, img_cut, pos=pos, stretch=stretch, /noerase, $
              cti=0, /rev
     cgplot, [0], [0], xs=1, ys=1, xr=xmnx, yr=ymnx, $
             pos=pos, /nodata, /noerase, xtitle=xtit, $
             yminor=2, yticks=4, xminor=2, xticks=4, $
             xthick=thick, ythick=thick, ytickint=20L,$
             xtickinterval=50, xcharsize=0.8, ycharsize=0.8

     ;; overlay the traces
     cgplot, trace_x(gd), trace_y_orig(gd), color='red', $
             thick=1, linest=1, /ov
     cgplot, trace_x(gd), trace_yfit(gd), color='cyan', $
             thick=1.2, linest=0, /ov
     badpix=badpix(where(badpix ne 0, nbad))
     if nbad ne 0 then $
        cgplot, trace_x(badpix), trace_yfit(badpix), color='black', $
                psym=16, symsize=0.2, /ov
     cgplot, trace_x(gd), trace_yfit(gd)+sigma, $
             color='orange', thick=1, linest=2, /ov
     cgplot, trace_x(gd), trace_yfit(gd)-sigma, $
             color='orange', thick=1, linest=2, /ov
     cgplot, trace_x(gd), trace_yfit(gd)+6*trace_sigma_fit(gd), $
             color='yellow', thick=1, linest=2, /ov
     cgplot, trace_x(gd), trace_yfit(gd)+6*trace_sigma_fit(gd)+19, $
             color='yellow', thick=1, linest=2, /ov
     cgplot, trace_x(gd), trace_yfit(gd)-6*trace_sigma_fit(gd), $
             color='yellow', thick=1, linest=2, /ov
     cgplot, trace_x(gd), trace_yfit(gd)-6*trace_sigma_fit(gd)-19, $
             color='yellow', thick=1, linest=2, /ov
     
     al_legend, bm, psym=15, color=clr, textcolor=clr, /right, /top, $
                bthick=thick, outline_color='black', $
                background='white', charsize=0.8     
  endfor

  ;; plot the position of the fit
  pos=[0.12,0.51,0.19,0.60]
  cgimage, img_trim, ctindex=0, pos=pos, $
           minval=-10, maxval=200, /noerase, /rev
  cgplot, [0], [0], xr=[-1*imgsz-0.5,imgsz+0.5], $
          yr=[-1*imgsz-0.5,imgsz+0.5], $
          /nodata, /noerase, pos=pos, xs=1, ys=1, xcharsize=0.2, $
          ycharsize=0.2, xtickname=replicate(' ',50), $
          ytickname=replicate(' ',50)
  cgplot, [wfc3_g280_strct(ii).x0-x0],[wfc3_g280_strct(ii).y0-y0], $
          psym=16, color='red', /ov, symsize=0.3
  ;cgplot, [wfc3_g280_strct(ii).xguess-wfc3_g280_strct(ii).x0], $
  ;        [wfc3_g280_strct(ii).yguess-wfc3_g280_strct(ii).y0], $
  ;        psym=16, color='red', /ov, symsize=0.3
     
  cgtext, 0.121,0.59, 'Pixel centroid = ('+ $
          strcompress(string(round(wfc3_g280_strct(ii).x0)),/r)+','+$
          strcompress(string(round(wfc3_g280_strct(ii).y0)),/r)+')', $
          color='black', charsize=0.29, /norm
  
  ;; plot the 1d spectra for beam A and beam C + combined and sky
  wgda=where(wfc3_g280_strct(ii).wavea gt 1800 and wfc3_g280_strct(ii).wavea lt 1d4)
  wgdc=where(wfc3_g280_strct(ii).wavec gt 1800 and wfc3_g280_strct(ii).wavec lt 1d4)
  wgd=where(wfc3_g280_strct(ii).wave gt 1800 and wfc3_g280_strct(ii).wave lt 1d4)
  badpixa=wfc3_g280_strct(ii).badpixa(where(wfc3_g280_strct(ii).badpixa ne 0))
  badpixc=wfc3_g280_strct(ii).badpixc(where(wfc3_g280_strct(ii).badpixc ne 0))
  minfx=min(wfc3_g280_strct(ii).fluxa(wgda)*1d17)<$
        min(wfc3_g280_strct(ii).fluxc(wgdc)*1d17)
  minsg=median(wfc3_g280_strct(ii).flux_siga(wgda)*1d17)>$
        median(wfc3_g280_strct(ii).flux_sigc(wgdc)*1d17)
  yr0=minfx>(-1d*minsg)
  bin=total(histogram(wfc3_g280_strct(ii).fluxa, nbin=1000, loc=xbin),/cum)
  yr1=xbin(max(where(bin lt 0.99*n_elements(wfc3_g280_strct(ii).fluxa))))*1d17*1.5
  pos=[0.10,0.10,0.95,0.40]
  
  cgplot, [0], [0], xs=1, ys=1, /nodata, /noerase, $
          xr=wvrg, yr=[yr0,yr1], pos=pos,$
          xthick=thick, ythick=thick, charsize=1.2
  cgplot, wvrg, [0,0], line=0, color='black', /ov

  ;; beam A
  cgplot, wfc3_g280_strct(ii).wavea(wgda), wfc3_g280_strct(ii).fluxa(wgda)*1d17, $
          psym=10, /ov, color='red3', thick=thick
  cgplot, wfc3_g280_strct(ii).wavea(wgda), wfc3_g280_strct(ii).flux_siga(wgda)*1d17, $
          psym=10, /ov, color='red3', thick=thick, line=1
  tmprga=[badpixa(0)]
  for jj=1L,n_elements(badpixa)-1 do begin
     if badpixa(jj) eq badpixa(jj-1)+1 and jj ne n_elements(badpixa)-1 then begin
        tmprga=[tmprga,badpixa(jj)]
        continue
     endif else begin
        cgplot, wfc3_g280_strct(ii).wavea(tmprga), wfc3_g280_strct(ii).fluxa(tmprga)*1d17, $
                psym=10, /ov, color='red1', thick=thick
        cgplot, wfc3_g280_strct(ii).wavea(tmprga), wfc3_g280_strct(ii).fluxa(tmprga)*1d17, $
                psym=10, /ov, color='white', thick=thick*0.6, line=2
     endelse
  endfor
  
  ;; beam C
  cgplot, wfc3_g280_strct(ii).wavec(wgdc), wfc3_g280_strct(ii).fluxc(wgdc)*1d17, $
          psym=10, /ov, color='blu3', thick=thick
  cgplot, wfc3_g280_strct(ii).wavec(wgdc), wfc3_g280_strct(ii).flux_sigc(wgdc)*1d17, $
          psym=10, /ov, color='blu3', thick=thick, line=1
  tmprgc=[badpixc(0)]
  for jj=1L,n_elements(badpixc)-1 do begin
     if badpixc(jj) eq badpixc(jj-1)+1 and jj ne n_elements(badpixc)-1 then begin
        tmprgc=[tmprgc,badpixc(jj)]
        continue
     endif else begin
        cgplot, wfc3_g280_strct(ii).wavec(tmprgc), wfc3_g280_strct(ii).fluxc(tmprgc)*1d17, $
                psym=10, /ov, color='blu1', thick=thick
        cgplot, wfc3_g280_strct(ii).wavec(tmprgc), wfc3_g280_strct(ii).fluxc(tmprgc)*1d17, $
                psym=10, /ov, color='white', thick=thick*0.6, line=2
        tmprgc=[badpixc(jj)]
     endelse
  endfor
  
  ;;  Combined
  cgplot, wfc3_g280_strct(ii).wave(wgd), wfc3_g280_strct(ii).flux(wgd)*1d17, $
          psym=10, /ov, color='pur7', thick=thick
  cgplot, wfc3_g280_strct(ii).wave(wgd), wfc3_g280_strct(ii).sig(wgd)*1d17, $
          psym=10, /ov, color='pur7', thick=thick, line=1

  ;; sky
  cgplot, wfc3_g280_strct(ii).wave(wgd), wfc3_g280_strct(ii).sky(wgd)*1d17, $
            psym=10, /ov, color='darkgreen', thick=thick, line=2
  
  ;; legend
  al_legend, ['Beam A','Beam C','Combined','Sky'], psym=[15,15,15,15], $
             colors=['red3','blu3','pur7','darkgreen'], /right, charsize=0.8, $
             bthick=thick, outline_color='black', /top, $
             textcolors=['red3','blu3','pur7','darkgreen']
  
  ;; define some text
  cgtext, 0.05,0.10, 'Flux ('+textoidl('10^{-17} erg s^{-1} cm^{-2} ')+cgsymbol('angstrom')+$
          textoidl('^{-1})'), /norm, color='black', orient=90, charsize=1.
  cgtext, 0.03,0.68, 'Row', /norm, color='black', orient=90, charsize=1.2
  cgtext, 0.47,0.03, 'Wavelength ('+cgsymbol('angstrom')+')', color='black', charsize=1.2, /norm
  ps_end

end
