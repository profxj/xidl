;+ 
; NAME:
; fuse_cog
;  V1.1
;
; PURPOSE:
;    Calculate a COG solution from a FUSE structure file
;
; CALLING SEQUENCE:
;   
;  fuse_cog, strct_fil, cog_fil, [Nlmt, blmt], /CHICHK, PLTONLY=
;    NSTP=, BSTP=, PSFILE=, OUTFIL=, /EXACT, ZLBL=,DEBLEND=,UPLIM=,RMS=
;
; INPUTS:
;  strct_fil -- FITS file for the FUSE structure
;  cog_fil -- COG input file (lists redshift and transitions to use)
;  [Nlmt] -- Range of column densities to explore (2 element array)
;  [blmt] -- Range of Doppler parameters to explore (2 element array)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /CHICHK -- Plot Chi^2 image
;  PLTONLY --  4-element array of N,b values and error for a plot
;  NSTP -- Number of steps to search N space [default: 100L]
;  BSTP -- Number of steps to search b space [default: 100L]
;  /EXACT -- Calculate EW exactly (Spline is generally good enough)
;  ZLBL= -- Label for Plot giving redshift of the absorber (string)
;  UPLIM -- file (like cog_fil) of lines to color upper limits on
;           curve (1: excluded, 2: 2-sig upper limit, 4: blend upper limit
;  DEBLEND -- file with information of EW to remove for specific 
;             features because they are blended OR include upper limits
;  RMS -- add RMS to EW error in quadrature
;  LABEL -- indicate line across top (either /label or label='FeII')
;  ASYMERR -- calculate asymmetric errors by measuring the chi^2=1 ellipse
;
; OPTIONAL OUTPUTS:
;   OUTFIL -- File with best fit values and error
;   PSFILE -- File for postscript plot
;
; COMMENTS:
;
; EXAMPLES:
;  fuse_cog, '/u/xavier/FUSE/data/PKS0405-12/Analysis/pks0405_abslin.fits', $
;    '/u/xavier/FUSE/data/PKS0405-12/Analysis/COG/Input/pks0405_z0918.cog', $
;    PSFIL='Figures/z0918_cog.ps', PLTONLY=[14.52, 38.2, 0.04, 1.8]
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   11-Sep-2003 Written by JXP
;   30-Nov-2005 modifed so can handle inability to calculate error, KLC
;   12-Jan-2006 added UPLIM keyword, KLC
;    8-Mar-2006 added DEBLEND keyword, KLC
;   13-Sep-2006 added RMS keyword, KLC
;   22-Sep-2006 added LABEL keyword, KLC
;    2-Jan-2007 added ASYMERR keyword, KLC and PJ; modify plot labels
;    8-Jan-2007 corrected and modified UPLIM
;-
;------------------------------------------------------------------------------
function lblnam,fnam,elem
;; Trim line name down
nn = n_elements(fnam)
lygrk = strarr(nn)
nam = strtrim(fnam,2)

for ii=0,nn-1 do begin
    case nam[ii] of
        'HI 930': lygrk[ii] = '!9z!3'
        '930':  lygrk[ii] = '!9z!3'
        'HI 937': lygrk[ii] = '!9e!3'
        '937':  lygrk[ii] = '!9e!3'
        'HI 949': lygrk[ii] = '!9d!3'
        '949':  lygrk[ii] = '!9d!3'
        'HI 972': lygrk[ii] = '!9g!3'
        '972':  lygrk[ii] = '!9g!3'
        'HI 1025': lygrk[ii] = '!9b!3'
        '1025':  lygrk[ii] = '!9b!3'
        'HI 1215': lygrk[ii] = '!9a!3'
        '1215':  lygrk[ii] = '!9a!3'
        else: $
          if keyword_set(elem) and size(elem,/type) eq 7 then $
          lygrk[ii] = strmid(nam[ii],strpos(nam[ii],elem)+strlen(elem)+1,$
                             strlen(nam[ii])) $
        else lygrk[ii] = strmid(nam[ii],strpos(nam[ii],'HI')+3,$
                                strlen(nam[ii]))
    endcase
endfor

return,lygrk
end


pro fuse_cog, strct_fil, cog_fil, Nlmt, blmt, CHICHK=chichk, PLTONLY=pltonly, $
              NSTP=nstp, BSTP=bstp, PSFILE=psfile, OUTFIL=outfil, EXACT=exact,$
              ZLBL=zlbl, UPLIM=uplim, DEBLEND=deblend, RMS=rms, LABEL=label, $
              ASYMERR=asymerr

  common x_calccog_cmm, cog_bval, cog_tau, cog_strct

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
      'fuse_cog, strct, cog_fil, [Nlmt, blmt], /CHICHK, /PLTONLY,' + $
      'NSTP=, BSTP=,  PSFILE=, OUTFIL=, /EXACT, ZLBL= (v1.1)' 
    return
  endif 

  if not keyword_set( NSTP ) then nstp = 100L
  if not keyword_set( BSTP ) then bstp = 100L
  if not keyword_set(EXACT) then begin
      print, 'fuse_cog: Using spline interpolation for the COG.  '
      print, '      This assumes a Maxwellian profile'
      cogmax_fil = getenv('XIDL_DIR')+'/Spec/Analysis/cogmax_tab.fits'
      cog_strct = xmrdfits(cogmax_fil, 1, /silent)
  endif
  if not keyword_set( LSIZE ) then lsize = 2.7

  compile_opt strictarr
  resolve_routine, 'x_calccog', /COMPILE_FULL_FILE, /EITHER
  
;  Open structure
  strct = xmrdfits(strct_fil, 1, /silent)
  msk = lonarr(n_elements(strct)) - 1L
  
; Read cog files
  close, /all
  readcol, cog_fil, val, FORMAT='D'
  zabs = val[0]
  wv_lin = val[1:*]
  nlin = n_elements(wv_lin)
  
  ;; Read list
  for i=0L,nlin-1 do begin
      ;; Mask
      a = where(abs(strct.zabs-zabs) LT 0.0002 AND $
                abs(strct.wrest-wv_lin[i]) LT 0.003, na)
      if na NE 1 then begin
          print, 'fuse_cog: no line', wv_lin[i], na
          stop
      endif
      msk[a] = a[0]
  endfor
  
  ;; Grab structure subset 
  gd = where(msk GE 0L)
  strct = strct[gd]
  
  ;; Reduce the EW (and put in Ang not mA)
  redew = strct.EW[0] / strct.wrest / 1000.
  if keyword_set(RMS) then redsigew = $
    sqrt(rms^2+strct.sigEW[0]^2) / strct.wrest / 1000. $
  else redsigew = strct.sigEW[0] / strct.wrest / 1000.


  ;; Remove blend EW (added by KLC)
  ;if keyword_set(deblend) then begin
  if size(deblend,/type) eq 7 then begin
     readcol,deblend,blndlin,blndwrest,blndew,blndewerr,format='(d,d,d,d)'
      nblnd = n_elements(blndlin)
      for ii=0,nblnd-1 do begin
          gd = where(abs(strct.wrest-blndlin[ii]) le 0.0005,ngd)
          case ngd of
              0: stop,'fuse_cog: blend not listed ',blndlin[ii]
              1: begin
                  redew[gd] = redew[gd]-blndew[ii]/ blndwrest[ii] / 1000.
                  redsigew[gd] = sqrt(redsigew[gd]^2+$
                                      (blndewerr[ii] /blndwrest[ii]/1000.)^2)
              end 
              else: stop,'fuse_cog: multiple matching blends ',blndlin
          endcase 
      endfor 
  endif 
  ;; End deblending (KLC)
  

  ;; Read in upper limit lines and use those in setting plot bounds
  ;; (added by KLC)
  if keyword_set(uplim) then begin
      strct_ul = xmrdfits(strct_fil, 1, /silent)
      msk_ul = replicate(-1,n_elements(strct_ul))
      close,/all

      ;;Flags: 1 excluded point; 2 2-sigma upper limit; 4 blend upper limit
      readcol,uplim,wv_lin_ul, flg_ul,format='d9.4,i1'      
      nlin_ul = n_elements(wv_lin_ul)
      srt = sort(wv_lin_ul)
      wv_lin_ul = wv_lin_ul[srt]
      flg_ul = flg_ul[srt]

      for ii=0L,nlin_ul-1 do begin
          a = where(abs(strct_ul.zabs-zabs) lt 0.0003 and $
                    abs(strct_ul.wrest-wv_lin_ul[ii]) lt 0.003,na)
          if na ne 1 then stop,'fuse_cog: no line ',wv_lin_ul[ii], na
          msk_ul[a] = a[0]
      endfor

      gd = where(msk_ul ge 0L)
      strct_ul = strct_ul[gd]
      srt = sort(strct_ul.wrest)
      strct_ul = strct_ul[srt]

      getfnam, strct_ul.wrest, fval_ul, nam_ul
      flambda_ul = fval_ul*strct_ul.wrest
      srt = sort(flambda_ul)       ;(cm)

      print,''
      printcol, strct_ul[srt].wrest, alog10(flambda_ul[srt]*1.e-8)
      xplt_ul = alog10(fval_ul*strct_ul.wrest*1.e-8)
      yplt_ul = alog10(strct_ul.ew[0]/strct_ul.wrest/1000.)
      if keyword_set(RMS) then ysig_ul = sqrt(rms^2+strct_ul.sigew[0]^2)/$
        (alog(10)*strct_ul.ew[0]) $
      else ysig_ul = strct_ul.sigew[0]/(alog(10)*strct_ul.ew[0])

      ul = where(flg_ul eq 2,nul)
      if nul ne 0 then begin
          yplt_ul[ul] = $       ;2-sigma upper limit
            alog10(2*strct_ul[ul].sigew[0]/(1000.*strct_ul[ul].wrest)) 
          
          if keyword_set(deblend) then begin
              ;; Include upper limits
              redew = [redew,replicate(0.,nul)]
              ;redew = [redew,strct_ul[ul].sigew[0]/1000./strct_ul[ul].wrest]
              redsigew = [redsigew,strct_ul[ul].sigew[0]/1000./strct_ul[ul].wrest]
              nlin=nlin+nul
              strct = [strct,strct_ul[ul]]
          endif 
      endif 

      ll = where(flg_ul eq 4,nll)
      ;if nll ne 0 then $
      ;  yplt_ul[ll] = yplt_ul[ll] ;blend upper limit
      excl = where(flg_ul eq 1,nexcl)

  endif 
  ;; ended if uplim (KLC)


  ;; Grab fvalues
  getfnam, strct.wrest, fval, nam
  flambda = strct.wrest*fval
  srt = sort(flambda) ;(cm)
  printcol, strct[srt].wrest, alog10(flambda[srt]*1.e-8)
  
  if not keyword_set( PLTONLY ) then begin
      ;; CHISQ array
      chisq = dblarr(NSTP,BSTP)
      
; LOOP
      
      calcEW = dblarr(nlin)
      
      for ii=0L,NSTP-1 do begin
          cog_Nval = Nlmt[0] + (Nlmt[1]-Nlmt[0])*float(ii)/float(nstp)
          if ii MOD 10 EQ 0 then print, ii
          for jj=0L,BSTP-1 do begin
              cog_bval = blmt[0] + (blmt[1]-blmt[0])*float(jj)/float(bstp)
              
              ;; Get TAU
              tau = 1.497e-2*(flambda*1.e-8)*(10^cog_Nval)/(cog_bval*1e5) 
              for kk=0L, nlin-1 do $
                calcEW[kk] = x_calccog_cog(tau[kk], EXACT=exact)
              ;; Calculate chisq
              chisq[ii,jj] = total( ((calcEW-redew)/(redsigew))^2 )
              
          endfor
      endfor
      if keyword_set( CHICHK ) then xatv, chisq/float((nlin-2)>1), /block
      
      ;; Find min chisq
      mn = min(chisq, ichi)
      print, 'fuse_cog: min chi_red = ', mn/float((nlin-2)>1)
      ii = ichi MOD NSTP
      jj = (ichi/NSTP)
      Ngd = Nlmt[0] + (Nlmt[1]-Nlmt[0])*float(ii)/float(nstp)
      bgd = blmt[0] + (blmt[1]-blmt[0])*float(jj)/float(bstp)
      cog_bval = bgd

      ;; Calculate the probability array
      prob_arr = 1. - chisqr_pdf(chisq, 2)
      tot_prob = total(prob_arr)

      ;; Error in N
      sz = size(prob_arr)
      for i=0L, 2*NSTP/3 do begin
          if ii-i ge 0 and ii+i lt sz[1] then begin
              sub_prob = total( prob_arr[(ii-i)>0:ii+i,*] )
              if sub_prob GT 0.683*tot_prob then break
          endif else begin
              print,'fuse_cog: cannot determine sigN ',sub_prob
              sigN = -9.99
              goto,assume_sigN
          endelse 
      endfor
      if i EQ 0 then stop
      sigN = (Nlmt[1]-Nlmt[0])*float(i)/float(nstp)
      assume_sigN:
      
      ;; Error in b
      for j=0L, 2*BSTP/3 do begin
          if jj-j ge 0 and jj+j lt sz[2] then begin
              sub_prob = total( prob_arr[*,(jj-j)>0:jj+j] )
              if sub_prob GT 0.683*tot_prob then break
          endif else begin
              print,'fuse_cog: cannot determine sigb ',sub_prob
              sigb = -9.99
              goto,assume_sigb
          endelse
      endfor
      if j EQ 0 then stop
      sigb = (blmt[1]-blmt[0])*float(j)/float(bstp)
      assume_sigb:
      
      print, 'fuse_cog: N = ', Ngd, '+/-', sigN
      print, 'fuse_cog: b = ', bgd, '+/-', sigb

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Calcuate asymmetric errors at chi^2=1 ellipse
      ;; Patrik Jonsson, 22 Dec 2006
      if keyword_set(asymerr) then begin
          pjunkarr1= where(chisq lt mn+1)
          pjunkarri= pjunkarr1 mod nstp
          pjunkarrj = pjunkarr1 / nstp
          pimin = min(pjunkarri,max=pimax)
          pjmin = min(pjunkarrj,max=pjmax)
          if pimin eq 0 or pimax eq sz[1]-1 then $
            stop,"Badness!: Error ellipse not closed in N: ",pimin,pimax
          if pimax-pimin le 5 then $
            stop,"Badness!: Error ellipse too small in N"
          if pjmin eq 0 or pjmax eq sz[2]-1 then $
            stop,"Badness!: Error ellipse not closed in b: ",pjmin,pjmax
          if pjmax-pjmin le 5 then $
            stop,"Badness!: Error ellipse too small in b"
          
          sigbm = (blmt[1]-blmt[0])*float(-pjmin+jj)/float(bstp)
          sigNm = (Nlmt[1]-Nlmt[0])*float(-pimin+ii)/float(nstp)
          sigbp = (blmt[1]-blmt[0])*float(pjmax-jj)/float(bstp)
          sigNp = (Nlmt[1]-Nlmt[0])*float(pimax-ii)/float(nstp)

          print, 'fuse_cog: N = ', Ngd, '+', sigNp,'/-',signm
          print, 'fuse_cog: b = ', bgd, '+', sigbp,'/-',sigbm

          if keyword_set(chichk) then begin
              contour,chisq-mn,level=[0.01,0.5,1,1.5,2,2.5],$
                xtitle='logN index',ytitle='b index',/font
              stop,'fuse_cog: view chi^2 contours ',cog_fil
          endif 

          sigN = sigNp > sigNm
          sigb = sigbp > sigbm 
      endif                     ;asymerr
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      Naodm = total(strct.ncolm/strct.signcolm^2)/total(1./strct.signcolm^2)
      sigNaodm = sqrt(total(strct.signcolm^2))
      print, 'fuse_cog: Naodm = ',Naodm,'+/-',sigNaodm
      ;; Write to file
      if keyword_set( OUTFIL ) then begin
          openw, 22, outfil
          if not keyword_set(asymerr) then begin
              printf, 22, Ngd, sigN
              printf, 22, bgd, sigb
          endif else begin
              printf, 22, Ngd, sigNp, -sigNm
              printf, 22, bgd, sigbp, -sigbm
          endelse 
          printf, 22, mn/float((nlin-2)>1) ;;reduced chi^2
          close, 22
      endif
  endif else begin
          Ngd = pltonly[0]
          bgd = pltonly[1]
      if keyword_set(asymerr) then begin
          sigNp = pltonly[2]
          sigbp = pltonly[3]
          sigNm = pltonly[4]
          sigbm = pltonly[5]
          sigN = sigNp > sigNm
          sigb = sigbp > sigbm
      endif else begin
          sigN = pltonly[2]
          sigb = pltonly[3]
      endelse 
  endelse

  cog_bval=bgd

  if keyword_set(uplim) and keyword_set(deblend) then begin
      if nul gt 0 then begin
          nlin = nlin-nul
          redew = redew[0:nlin-1]
          redsigew = redsigew[0:nlin]
          strct = strct[0:nlin-1]
      endif 
  endif 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plot
  xplt = alog10(flambda*1e-8)  ; flambda (cm)
  yplt = alog10(redew)  ; log (EW/lambda)
  ysig = redsigew/(alog(10)*10^yplt)

  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs,/encapsulated

  xmn = min(xplt, MAX=xmx) - 0.1
   xmx = xmx + 0.1
  ymn = min(yplt, max=ymx)

  if keyword_set(uplim) then begin
      xmn = xmn < min(xplt_ul) - 0.1
      xmx = xmx > max(xplt_ul) + 0.1
      ymn = ymn < min(yplt_ul) 
      ;if nul ne 0 then ymn = ymn < min(yplt_ul[ul]-2*ysig_ul[ul])
      ymx = ymx > max(yplt_ul)
  endif 

  ymnx = [ymn-0.1, ymx+0.1]

  clr = getcolor(/load)
  plot, [xmn,xmx], ymnx, /nodata, background=clr.white, $ 
    color=clr.black, xtitle='log(f!9l!X)', $
    ytitle='log(W!Dr!N /!9l!X)', xstyle=1, ystyle=1, charsize=2.0, $
    xmargin=[8,2], ymargin=[4,1], xtickinterval=0.5
  oploterror, xplt, yplt, ysig, psym=1, color=clr.black, ERRCOLOR=clr.black

  if keyword_set(label) then begin
      ytmp = replicate(ymnx[1]-0.025*(ymnx[1]-ymnx[0]), n_elements(xplt))
      tmp = lblnam(nam,label)
      xyouts,xplt,ytmp,tmp,orientation=-90,color=clr.black,/font,$
        charsize=2.
  endif 


  ;; plot excluded points and upper limits of other transitions (Added
  ;; by KLC) 
  if keyword_set(uplim) then begin
      if nexcl ne 0 then $ ;Excluded (triang)
        oploterror, xplt_ul[excl], yplt_ul[excl], ysig_ul[excl], psym=5, $
        color=clr.cyan, ERRCOLOR=clr.cyan 
      plotsym,1,4,thick=2.5      ;down arrow, now psym=8
      if nul ne 0 then begin    ;2-sigma upper limits
          oplot,xplt_ul[ul],yplt_ul[ul],psym=8,color=clr.magenta
          oplot,xplt_ul[ul],yplt_ul[ul],psym=6,color=clr.magenta ;square
      endif 
      if nll ne 0 then begin    ;blend upper limits
          oplot,xplt_ul[ll],yplt_ul[ll],psym=8,color=clr.cyan
          oplot,xplt_ul[ll],yplt_ul[ll],psym=7,color=clr.cyan ;X
      endif 

      if keyword_set(label) then begin
          ytmp = replicate(ymnx[1]-0.025*(ymnx[1]-ymnx[0]), n_elements(xplt_ul))
          tmp = lblnam(nam_ul,label)
          xyouts,xplt_ul,ytmp,tmp,orientation=-90,color=clr.black,/font,$
            charsize=2.0
      endif
  endif
  ;;end /uplim (KLC)
  
  nplt = 100L
  cog = dblarr(nplt)

  ;; Good value
  xplt2 = (xmn + (xmx-xmn)*findgen(nplt)/float(nplt))
  tau = 1.497e-2*(10^xplt2)*(10^Ngd)/(bgd*1e5)
  for q=0L,nplt-1 do cog[q] = x_calccog_cog(tau[q])
  oplot, xplt2, alog10(cog), color=clr.red

  Nsv = Ngd
  bsv = bgd

  ;; Error edges (high N, high b)
  if not keyword_set(asymerr) then begin
      Ngd = Ngd + sigN
      bgd = bgd + sigb
  endif else begin
      Ngd = Ngd + sigNp
      bgd = bgd + sigbp
  endelse

  cog_bval=bgd
  tau = 1.497e-2*(10^xplt2)*(10^Ngd)/(bgd*1e5)
  for q=0L,nplt-1 do cog[q] = x_calccog_cog(tau[q])
  oplot, xplt2, alog10(cog), color=clr.green, linestyle=2

  ;; Error edges (low N, low b)
  if not keyword_set(asymerr) then begin
      Ngd = Ngd - 2*sigN
      bgd = bgd - 2*sigb
  endif else begin
      Ngd = Ngd - sigNp - sigNm
      bgd = bgd - sigbp - sigbm
  endelse
 
  cog_bval=bgd
  tau = 1.497e-2*(10^xplt2)*(10^Ngd)/(bgd*1e5)
  for q=0L,nplt-1 do cog[q] = x_calccog_cog(tau[q])
  oplot, xplt2, alog10(cog), color=clr.green, linestyle=2

  if keyword_set(outfil) then begin
      ;; Calculate predicted EW and sigEW
      ysig2 = fltarr(nlin)
      yplt2 = fltarr(nlin)

      ;; Lower 1-sigma limit (cog_bval set previously)
      tau = 1.497e-2*(10^xplt)*(10^Ngd)/(bgd*1e5)
      for q=0L,nlin-1 do ysig2[q] = x_calccog_cog(tau[q])

      ;; Predicted COG
      tau = 1.497e-2*(10^xplt)*(10^Nsv)/(bsv*1e5)
      cog_bval = bsv
      for q=0L,nlin-1 do yplt2[q] = x_calccog_cog(tau[q])
      
      ;; Difference (as estimate of error)
      ysig2 = yplt2 - ysig2

      ;; Convert COG to EW and sigEW (mA)
      yplt2 = strct.wrest*yplt2*1000.
      ysig2 = strct.wrest*ysig2*1000.

      ;; Append to file
      openu,22,outfil,/append
      if keyword_set(rms) then printf,22,'RMS = ',rms,' mA'

      srt = sort(strct.wrest)
      
      printf,22,'#f*wv','wrest','ew_M','sigew','ew_F','sigew','diff','err',$
        format='(a5,2x,a9,2x,a8,1x,a8,2x,a8,1x,a8,2x,a10,1x,a10)'
      diff = abs(strct.ew[0]-yplt2)
      if keyword_set(RMS) then differr = diff*$
        sqrt((sqrt(rms^2+strct.sigew[0]^2)/strct.ew[0])^2 + $
             (ysig2/yplt2)^2) $
      else differr = diff*sqrt((strct.sigew[0]/strct.ew[0])^2 + $
                               (ysig2/yplt2)^2)
      if keyword_set(RMS) then writecol,outfil,xplt[srt],strct[srt].wrest,$
        strct[srt].ew[0],sqrt(rms^2+strct[srt].sigew[0]^2),yplt2[srt],$
        ysig2[srt],diff[srt],differr[srt],$
        fmt='(f5.2,2x,f9.4,2x,f8.2,1x,f8.2,2x,f8.2,1x,f8.2,2x,f10.4,1x,f10.4)',$
        filnum=22 $
      else writecol,outfil,xplt[srt],strct[srt].wrest,strct[srt].ew[0],$
        strct[srt].sigew[0],yplt2[srt],ysig2[srt],diff[srt],differr[srt],$
        fmt='(f5.2,2x,f9.4,2x,f8.2,1x,f8.2,2x,f8.2,1x,f8.2,2x,f10.4,1x,f10.4)',$
        filnum=22

      if keyword_set(uplim) then begin
          ysig_ul = fltarr(nlin_ul)
          yplt_ul = fltarr(nlin_ul)

          ;; Predicted COG
          tau = 1.497e-2*(10^xplt_ul)*(10^Nsv)/(bsv*1e5)
          for q=0L,nlin_ul-1 do yplt_ul[q] = x_calccog_cog(tau[q])
      
          ;; Lower 1-sigma limit
          tau = 1.497e-2*(10^xplt_ul)*(10^Ngd)/(bgd*1e5)
          cog_bval = bgd
          for q=0L,nlin_ul-1 do ysig_ul[q] = x_calccog_cog(tau[q])

          ;; Difference as estimate of error
          ysig_ul = abs(yplt_ul - ysig_ul)
          
          ;; Convert COG to EW and sigEW (mA)
          yplt_ul = strct.wrest*yplt_ul*1000.
          ysig_ul = strct.wrest*ysig_ul*1000.

          srt = sort(strct_ul.wrest)
          diff = abs(strct_ul.ew[0]-yplt_ul)
          if keyword_set(RMS) then differr = diff*$
            sqrt((sqrt(rms^2+strct_ul.sigew[0]^2)/strct_ul.ew[0])^2 + $
                 (ysig_ul/yplt_ul)^2) $
          else differr = diff*sqrt((strct_ul.sigew[0]/strct_ul.ew[0])^2 + $
                                   (ysig_ul/yplt_ul)^2)

          printf,22,''
          printf,22,'# Upper limits'
          if keyword_set(rms) then writecol,outfil,xplt_ul[srt],$
            strct_ul[srt].wrest,strct_ul[srt].ew[0],$
            sqrt(rms^2+strct_ul[srt].sigew[0]^2),$
            yplt_ul[srt],ysig_ul[srt],diff[srt],differr[srt],$
            fmt='(f5.2,2x,f9.4,2x,f8.2,1x,f8.2,2x,f8.2,1x,f8.2,2x,f10.4,1x,f10.4)',$
            filnum=22 $
          else writecol,outfil,xplt_ul[srt],strct_ul[srt].wrest,$
            strct_ul[srt].ew[0],strct_ul[srt].sigew[0],$
            yplt_ul[srt],ysig_ul[srt],diff[srt],differr[srt],$
            fmt='(f5.2,2x,f9.4,2x,f8.2,1x,f8.2,2x,f8.2,1x,f8.2,2x,f10.4,1x,f10.4)',$
            filnum=22
      endif 

      close,22
  endif 


  ;; Label
  if keyword_set(ZLBL) then xyouts, 0.55, 0.37, zlbl, /normal, charsize=lsize
  if not keyword_set(asymerr) then begin
      xyouts, 0.55, 0.29, 'log N!DHI!N = '+string(Nsv,FORMAT='(f5.2)')+'!9 '+ $
        string("261B)+' !7'+$
               string(sigN,FORMAT='(f4.2)'), charsize=lsize, /normal
      xyouts, 0.55, 0.21, 'b = '+string(round(bsv),FORMAT='(i2)')+'!9 '+$
        string("261B)+ ' !7'+string(sigb,FORMAT='(f3.1)')+' km s!U-1!N', $
               charsize=lsize, /normal
  endif else begin
      xyouts, 0.55, 0.29, 'log N!DHI!N = '+string(Nsv,FORMAT='(f5.2)')+'!S!U+'+$
        string(sigNp,FORMAT='(f4.2)')+'!N!R!D-'+$
        string(sigNm,FORMAT='(f4.2)')+'!N', charsize=lsize, /normal
      if sigbp ge 10. or sigbm ge 10. then $
        xyouts, 0.55, 0.21, 'b = '+string(round(bsv),FORMAT='(i2)')+'!S!U+'+$
        string(sigbp,FORMAT='(f4.1)')+'!N!R!D-'+$
        string(sigbm,FORMAT='(f4.1)')+'!N km s!U-1!N', $
        charsize=lsize, /normal $
      else xyouts, 0.55, 0.21, 'b = '+string(round(bsv),FORMAT='(i2)')+'!S!U+'+$
        string(sigbp,FORMAT='(f3.1)')+'!N!R!D-'+$
        string(sigbm,FORMAT='(f3.1)')+'!N km s!U-1!N', $
        charsize=lsize, /normal
  endelse 

  if keyword_set( PSFILE ) then x_psclose

  return
end
