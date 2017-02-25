;+ 
; NAME:
; sdss_fndlls
;    Version 1.0
;
; PURPOSE:
;    Searches SDSS QSO spectra for LLS absorption.
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; RETURNS:
;   
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
;PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   2005 Written by JO
;-
;------------------------------------------------------------------------------
; help, sdss_findlls(template, 2.8, 'sdss3.fit.gz')
Function sdss_findlls,template,QSOtempwav0,SDSSspec, $
  btrial=BTRIAL,npoly=NPOLY,PLOTFIT=PLOTFIT, taucut=TAUCUT,$
  PLOTCHI=PLOTCHI,SPECPLOT=SPECPLOT,nosilent=NOSILENT,$
  ALLMOD=allmod, FSTRCT=fstrct, ZEM=zem, MOCK_DATA=mock_data

;state
llsstate = {  plate: 0L, $
              mjd: 0L, $
              fiber: 0L, $
              zem: 0.0, $
              zlls: 0.0, $
              taulls: 0.0, $
              blls: 0.0, $
              llsflg: 0, $
              zstart: 0.0, $
              zend: 0.0, $
              svchi: fltarr(100), $
              svz: fltarr(100), $
              uminusg: 0.0} 
           

if not keyword_set(ALLMOD) then begin
;    nnhi = 4
;    nh1 = xmrdfits('nhi16.5_18b30.fits',0)
    nnhi = 20
    nh1 = xmrdfits(getenv('XIDL_DIR')+'/SDSS/LLS/nhi16_19b30.fits',0) 
    npnh = n_elements(nh1) 
    allnh = fltarr(npnh,nnhi)
    for i=0L,nnhi-1 do allnh[*,i] = $
      xmrdfits(getenv('XIDL_DIR')+'/SDSS/LLS/nhi16_19b30.fits',i)
    
    ;; Create model image
    nztest = long((alog10(911.7633*6) - alog10(911.7833*4.2))/0.0001) + 11
    zllsarr = (10^(findgen(nztest)*0.0001 + alog10(4.2*911.7633)))/911.7633 $
              - 1.0
    wv_mod = 10^(2.7d + dindgen(npnh)*1e-4)
    mn = min(abs(wv_mod-911.7633),mmn)

    wv_sdss = 10^(3.58 + dindgen(3900)*1e-4)
    mn = min(abs(4.2*911.7633-wv_sdss),imn)
    shft = imn-mmn
    
    allmod = fltarr(npnh, nztest, nnhi)
    for rr=0l,nnhi-1 do begin
        for qq=0L,nztest-1 do begin
            allmod[*,qq,rr] = shift(allnh[*,rr], shft+qq)
            allmod[(npnh+qq+shft)<(npnh-1):*,qq,rr] = 1.
       endfor
    endfor
ENDIF
    
    
;;Defaults
IF NOT keyword_set(btrial) THEN btrial=[30.]
IF NOT keyword_set(npoly) THEN npoly=2
IF NOT keyword_set(taucut) THEN taucut=0.0

;;open template spectrum
tlogwav = findgen(n_elements(template))*0.0001 + QSOtempwav0

if not keyword_set(MOCK_DATA) then begin
    ;;open SDSS spectrum
    spspec=xmrdfits(SDSSspec,0,hdr,/silent)
    spec = spspec[*,0]
    specerr = spspec[*,2]
    speczem = sxpar(hdr,"Z")
    specwav0 = sxpar(hdr,"COEFF0")
;dwave=0.0001

    ;;read in some of the state variables
    llsstate.plate = sxpar(hdr,"PLATEID")
    llsstate.fiber = sxpar(hdr,"FIBERID")
    llsstate.mjd   = sxpar(hdr,"MJD")
    if not keyword_set(ZEM) then speczem = sxpar(hdr,"Z") else speczem = zem
    llsstate.zem   = speczem

    ;;get the colors
    color = sxpar(hdr,"MAG")
    colorsplit=strsplit(color)
    color_u=float(strmid(color,colorsplit[0],7))
    color_g=float(strmid(color,colorsplit[1],7))
    color_r=float(strmid(color,colorsplit[2],7))
    llsstate.uminusg = color_u - color_g
endif else begin
    spec = mock_data.spec
    specerr = mock_data.sig
    bad = where(specerr GT 999.,nbad)
    if nbad GT 0 then specerr[bad] = 0.
    specwav0 = mock_data.w0
    speczem = ZEM
endelse

logwav = findgen(n_elements(spec))*0.0001 + specwav0

;;;redshift the template
;tzlogwav = alog10((1.0 + speczem) * (10.0^tlogwav))


;;;min wavelength/z condition
llszmin = 10^specwav0/911.7633 - 1.0
IF (llszmin GT speczem+0.1) THEN BEGIN
    print,"No soup for you, zem lower than first pixel!!!"
    llsstate.llsflg=-1
    stop
    return,llsstate
ENDIF

;;;range in redshift/wavelegth for search
llswavmin = specwav0
llswavmax = alog10((1. + speczem)*911.7633)

;;;pick the template minimum pixel to extract from 
;pixmin = long((specwav0 - QSOtempwav0 - alog10(1. + speczem))/dwave)

;;;You can unweight the forest by dividing the forest by 0.7 if you include
;;;the forest in the regions which effect the slope of the template
;tempsec = extrac(template,pixmin,n_elements(spec))
;lyafrange=where(10^logwav LT (1.0 + speczem)*1215.67)
;tempsec(lyafrange)=tempsec(lyafrange) / 0.7

;;;Find the good pixels
;tempmask = (tempsec GT 0 AND specerr GT 0)
;good=where(tempmask)
;res=ladfit(logwav[good],spec[good]/tempsec[good])
;modelfit=poly(logwav,res)

;;;This part of the code slaps the slope on the SDSS template
;;;Here you can set the model matching range to either be where there are good
;;;pixels OR where there are good pixels AND wave > lya-emission

;good2=where(logwav GT llswavmax AND tempsec GT 0 and specerr GT 0)
;good2=where(logwav GT alog10((1.0 + speczem)*1215.67) $
;            AND tempsec GT 0 and specerr GT 0)
; res2=ladfit(logwav[good2],spec[good2]/tempsec[good2])
; modelfit2=poly(logwav,res2)


 ;nh1 = xmrdfits('nhi16.5_18b30.fits',0)
 sz = size(allmod, /dimensions)
 npnh = sz[0]
 nztest = sz[1]
 ;allnh = fltarr(npnh,4)
 ;for i=0L,3 do allnh[*,i] = xmrdfits('nhi16.5_18b30.fits',i)
 wv_mod = 10^(2.7d + dindgen(npnh)*1e-4)
 mn = min(abs(wv_mod-911.7633),mmn)

;loop it brotha, loop it
;FOR i=0,nztest-1 DO BEGIN
;    IF keyword_set(nosilent) THEN BEGIN
;        print,'Now doing bin '+string(i)+' of '+string(nztest)
;    ENDIF
;    bestchi=1.0E38
;    beststat2=1.0E38

 ;; Make redshifted array
 wv_sdss = 10^(3.58 + dindgen(3900)*1e-4)
 mn = min(abs(4.2*911.7633-10^logwav),imn)
; nztest = long((alog10(911.7633*6) - alog10(911.7833*4.2))/0.0001) + 11
; zllsarr = (10^(findgen(nztest)*0.0001 + alog10(4.2*911.7633)))/911.7633 - 1.0

; shft = imn-mmn


;fpx = round((specwav0 - 3.58)*1e4)
;IF (fpx LT 0) THEN fpx=0
lend = alog10((1.0 + speczem)*1220.0)
;epx = round((lend - 3.58)*1e4)
;normalize the slope, take the median, divide out the mean flux from
;an extrapolation to the kirkman 05 DA fit.

templ_spec = sdss_qsotempl(template, QSOtempwav0, logwav, spec, $
                           specerr, speczem, $
                           EPX=epx, FPX=fpx)
;specimg = (spec/modelfit2/tempsec)[0:epx-fpx] # replicate(1.,sz[1])
specimg = (spec[0:epx-fpx]/templ_spec)  # replicate(1.,sz[1])


medrange=where((10^logwav[0:epx-fpx])/(1.+speczem) GT 1050.0 AND $
               (10^logwav[0:epx-fpx])/(1.+speczem) LT 1170.0)
;medval=median(specimg[medrange,0])
;specimg=specimg/medval

;; new stuff, 2/28/07, finds the mean forest redshift between 1050 and 1170
;; then looks up the kirkman 05 DA, and divides that out
;medforz = median(10^logwav[medrange]/1215.67 -1.0)
;kirkman_da = 0.0062*(1.0 + medforz)^(2.75)
;dacorr = (1.0 - kirkman_da)

;specimg = specimg / dacorr

zllsarr = (10^(findgen(sz[1])*0.0001 + alog10(4.2*911.7633)))/911.7633 - 1.0

;stop
bestchi=1E38
;; Loop on NHI
for rr=0L,sz[2]-1 DO BEGIN
    chisquared=total(((specimg - allmod[fpx:epx,*,rr]))^2,1)
    redchi = chisquared / float(epx-fpx+1)
    mn = min(redchi,minchi)
;    minchi=where(redchi EQ min(redchi))
    llsstate.svchi[rr] = min(redchi)
    llsstate.svz[rr] = zllsarr[minchi]

    if (min(redchi) LT bestchi) THEN BEGIN
        savmod=allmod[fpx:epx,minchi,rr]
        chi2=redchi
        bestchi=min(redchi)
        savh1=rr
        zguess=zllsarr[minchi]
        bestchipix=minchi
    endif
endfor

;;;plot the chi-sq and a polynomial fit to the chi-sq
IF keyword_set(PLOTCHI) THEN BEGIN
    window,0
    plot,zllsarr,chi2,yrange=[0.95*min(chi2),1.05*max(chi2)],ystyle=1
    oplot,[speczem,speczem],[0.95*min(chi2),1.05*max(chi2)],linestyle=2
    oplot,[zguess,zguess],[0.95*min(chi2),1.05*max(chi2)],linestyle=3,color=255
    xyouts,zguess+0.001,0.9*max(chi2),'LLS',color=255
    xyouts,speczem+0.001,0.9*max(chi2),'zem'
ENDIF

IF keyword_set(SPECPLOT) THEN BEGIN
    x_specplot,SDSSspec,inflg=5,zin=zguess,/lls
ENDIF


;;;plot the data, template+lls, and smoothed data (for low snr)
IF keyword_set(PLOTFIT) THEN BEGIN
    clr = getcolor(/load)
    window,1
    plotmodel=savmod*templ_spec
    
    plot,10^logwav[0:epx-fpx],spec[0:epx-fpx], $
         position=[0.1,0.1,0.9,0.9], $
         xmin=specwav0,xstyle=1,yrange=[-0.5,max(plotmodel)],ystyle=1
    oplot,10^logwav[0:epx-fpx], plotmodel, color=clr.red
          
    ;;plot subwindows with lya and LLS region

    plotwavlls=(1.0 + zguess)*912.0
    
    plot,10^logwav[0:epx-fpx], $
         specimg[*,bestchipix], $
         xrange=[plotwavlls-100.,plotwavlls+100],xstyle=1, $
         position=[0.7,0.7,0.9,0.9],/noerase,charsize=0.75,ymin=0,ystyle=1
    oplot,10^logwav[0:epx-fpx],savmod,color=clr.red
    oplot,10^logwav,fltarr(n_elements(spec)),color=clr.green,linestyle=2
    
    plot,10^logwav,spec,xrange=[(1.0+zguess)*1215.67 - 100., $
                                (1.0+zguess)*1215.67 + 100.],xstyle=1, $
         position=[0.12,0.7,0.32,0.9],/noerase,charsize=0.75,ymin=0,ystyle=1
    oplot,[(1.0+zguess)*1215.67,(1.0+zguess)*1215.67],[0,100],color=clr.green
    
    plottitle1=SDSSspec
    plottitle2='zem='+strcompress(string(speczem))+$
               '  zlls='+strcompress(string(zguess),/remove_all) + $
               ' HI_idx=' + strcompress(string(savh1),/remove_all) 
    xyouts,0.1,0.95,plottitle1,/normal
    xyouts,0.1,0.92,plottitle2,/normal
ENDIF


  if arg_present(FSTRCT) then begin
      plotmodel=savmod*templ_spec
;      plotmodel=savmod*modelfit2[0:epx-fpx]*tempsec[0:epx-fpx]*medval
      fstrct = { $
                 zlls: zllsarr, $
                 chi2: chi2, $
                 zem: speczem, $
                 model: plotmodel, $
                 wv: 10^logwav[0:epx-fpx], $
                 spec: spec[0:epx-fpx], $
                 sig: specerr[0:epx-fpx], $
                 xmin: specwav0, $
                 speci: specimg[*,bestchipix], $
                 savmod: savmod $
               }
  endif
                 

  llsstate.llsflg = 1
  llsstate.zlls = zguess
  llsstate.taulls = savh1
  llsstate.blls = 30.
  llsstate.zstart = zllsarr[0]
  llsstate.zend = zllsarr[nztest - 11] ;because of 10 pix buffer before

  return,llsstate
END
