;+ 
; NAME:
; fuse_twocog
;  V1.1
;
; PURPOSE:
;    Given a list of DLA base files, fill up the structure ;
; CALLING SEQUENCE:
;   
;   lowzovi_prsdat, stucture, filename
;
; INPUTS:
;
; RETURNS:
;   structure      - IDL structure
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  LIST - File
;  ION - Input ionic column densities
;  NOELM - Supress inputting Elemental values
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   fuse_calccolm, struct, fil_instr
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   15-Sep-2003 Written by JXP
;-
;------------------------------------------------------------------------------
pro fuse_twocog_initcomm

common fuse_twocog_cmm, cog_b1, cog_N1, cog_b2, cog_N2, tcog_const, tcog_dv

return

end

;;;
function fuse_twocog_cog, tau, EXACT=exact

common fuse_twocog_cmm

  ;; Integrate (EXACT)
  ftau = qromo('fuse_twocog_ftau', -300., 300., /double)

  ; Other factors
  return, ftau / 3e5
end
  
;;;;
function fuse_twocog_ftau, x

common fuse_twocog_cmm

  stp1 = cog_N1 * exp(-x*x*cog_b1)
  stp2 = cog_N2 * exp(-(x-tcog_dv)*(x-tcog_dv)*cog_b1)
  ; Step 3
  stp3 = tcog_const*(stp1 + stp2)
  ; Integrand
  if stp3 LT -80.d then ftauint = 1. else $
    ftauint = (1.d - exp(stp3))
  return, ftauint
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro fuse_twocog, strct_fil, cog_fil, N1lmt, b1lmt, N2lmt, b2lmt, delv, $
                 CHICHK=chichk, PLTONLY=pltonly, $
                 NSTP=nstp, BSTP=bstp, PSFILE=psfile, OUTFIL=outfil, EXACT=exact

  common fuse_twocog_cmm

  if (N_params() LT 4) then begin 
    print,'Syntax - ' + $
             'fuse_twocog, strct, cog_fil (v1.0)' 
    return
  endif 

  if not keyword_set( NSTP ) then nstp = 5L
  if not keyword_set( BSTP ) then bstp = 5L
;  if not keyword_set(EXACT) then begin
;      print, 'fuse_cog: Using spline interpolation for the COG.  '
;      print, '      This assumes a Maxwellian profile'
;      cogmax_fil = getenv('XIDL_DIR')+'/Spec/Analysis/cogmax_tab.fits'
;      cog_strct = xmrdfits(cogmax_fil, 1, /silent)
;  endif

  tcog_dv = delv

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
      a = where(abs(strct.zabs-zabs) LT 0.001 AND $
                abs(strct.wrest-wv_lin[i]) LT 0.003, na)
      if na NE 1 then stop
      msk[a] = a[0]
  endfor
  
  ;; Grab structure subset
  gd = where(msk GE 0L)
  strct = strct[gd]
  
  ;; Reduce the EW (and deal with mA)
  redew = strct.EW[0] / (strct.zabs+1.) / strct.wrest / 1000.
  redsigew = strct.sigEW[0] / (strct.zabs+1.) / strct.wrest / 1000.
  
  ;; Grab fvalues
  getfnam, strct.wrest, fval
  flambda = strct.wrest*fval
  
  if not keyword_set( PLTONLY ) then begin
      ;; CHISQ array
      chisq = dblarr(NSTP,BSTP,NSTP,BSTP)
      
; LOOP
      
      calcEW = dblarr(nlin)
      
      for ii=0L,NSTP-1 do begin
          N1 = N1lmt[0] + (N1lmt[1]-N1lmt[0])*float(ii)/float(nstp)
          print, 'ii', ii
;          if ii MOD 10 EQ 0 then print, ii
          for jj=0L,BSTP-1 do begin
              b1 = b1lmt[0] + (b1lmt[1]-b1lmt[0])*float(jj)/float(bstp)
              for mm=0L,NSTP-1 do begin
                  N2 = N2lmt[0] + (N2lmt[1]-N2lmt[0])*float(mm)/float(nstp)
                  for nn=0L,BSTP-1 do begin
                      b2 = b2lmt[0] + (b2lmt[1]-b2lmt[0])*float(nn)/float(bstp)

                      ;; Reset N1, N2
                      cog_b1 = 1./b1^2 ;; Inverse to speed up calculation
                      cog_b2 = 1./b2^2 ;; Inverse to speed up calculation
                      cog_N1 = 10^N1 / (b1*1.e5)
                      cog_N2 = 10^N2 / (b2*1.e5)
                      
                      ;; Get REW
                      for kk=0L, nlin-1 do begin
                          tcog_const = -1.*1.497e-2*(flambda[kk]*1.e-8)
;                          xpnt = -100+ findgen(200)
;                          ypnt = fltarr(200)
;                          for tt=0L,200-1 do begin
;                              ypnt[tt] = fuse_twocog_ftau(xpnt[tt]) / 3e5
;                          endfor
;                          x_splot, xpnt, ypnt, /block
                          calcEW[kk] = fuse_twocog_cog() 
                      endfor
                      ;; Calculate chisq
                      chisq[ii,jj,mm,nn] = total( ((calcEW-redew)/(redsigew))^2 )
;      xplt = alog10(flambda*1e-8) ; flambda (cm)
;      yplt = alog10(redew)      ; log (EW/lambda)
;;      ysig = redsigew/(alog(10)*10^yplt)
;      xmn = min(xplt, MAX=xmx) - 0.1
;      xmx = xmx + 0.1
;      ymn = min(yplt, max=ymx)
;      ymnx = [ymn-0.1, ymx+0.1]
;      
;      clr = getcolor(/load)
;      plot, [xmn,xmx], ymnx, /nodata, background=clr.white, $ 
;        color=clr.black, xthick=3.0, ythick=3.0, xtitle='!17 log!d10!n(f!7k!X)', $
;        ytitle='log!d10!n(W/!7k!X)', xstyle=1, ystyle=1, charsize=1.9
;      oploterror, xplt, yplt, ysig, psym=1, color=clr.blue, ERRCOLOR=clr.blue
;      oplot, xplt, alog10(calcEW), color=clr.red
;      stop
      
                  endfor
              endfor
          endfor
      endfor
;      if keyword_set( CHICHK ) then xatv, chisq/float((nlin-2)>1), /block

      ;; Find min chisq
      mn = min(chisq, ichi)
      ii = ichi MOD NSTP
      jj = (ichi MOD (NSTP*BSTP))/NSTP
      mm = (ichi MOD (NSTP*BSTP*NSTP))/(NSTP*BSTP)
      nn = ichi / (NSTP*BSTP*NSTP)
      N1gd = N1lmt[0] + (N1lmt[1]-N1lmt[0])*float(ii)/float(nstp)
      b1gd = b1lmt[0] + (b1lmt[1]-b1lmt[0])*float(jj)/float(bstp)
      N2gd = N2lmt[0] + (N2lmt[1]-N2lmt[0])*float(mm)/float(nstp)
      b2gd = b2lmt[0] + (b2lmt[1]-b2lmt[0])*float(nn)/float(bstp)
;      cog_bval = bgd

      ;; Calculate the probability array
;      prob_arr = 1. - chisqr_pdf(chisq, 2)
;      tot_prob = total(prob_arr)

      ;; Error in N
;      for i=0L, 2*NSTP/3 do begin
;          sub_prob = total( prob_arr[ii-i:ii+i,*] )
;          if sub_prob GT 0.683*tot_prob then break
;      endfor
;      if i EQ 0 then stop
;      sigN = (Nlmt[1]-Nlmt[0])*float(i)/float(nstp)
      
      ;; Error in b
;      for j=0L, 2*BSTP/3 do begin
;          sub_prob = total( prob_arr[*,jj-j:jj+j] )
;          if sub_prob GT 0.683*tot_prob then break
;      endfor
;      if j EQ 0 then stop
;      sigb = (blmt[1]-blmt[0])*float(j)/float(bstp)
      
      print, 'fuse_cog: N1 = ', N1gd
      print, 'fuse_cog: N2 = ', N2gd
      print, 'fuse_cog: b1 = ', b1gd
      print, 'fuse_cog: b2 = ', b2gd
      ;; Write to file
;      if keyword_set( OUTFIL ) then begin
;          openw, 22, outfil
;          printf, 22, Ngd, sigN
;          printf, 22, bgd, sigb
;          close, 22
;      endif
  endif else begin
      stop
      Ngd = pltonly[0]
      bgd = pltonly[1]
      sigN = pltonly[2]
      sigb = pltonly[3]
  endelse

;  cos_bval=bgd

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plot
  xplt = alog10(flambda*1e-8)  ; flambda (cm)
  yplt = alog10(redew)  ; log (EW/lambda)
  ysig = redsigew/(alog(10)*10^yplt)

  if keyword_set( PSFILE ) then begin
      device, decompose=0
      ps_open, filename=PSFILE, /color, bpp=8
      !p.thick = 5
      !p.charthick = 3
  endif

  xmn = min(xplt, MAX=xmx) - 0.1
  xmx = xmx + 0.1
  ymn = min(yplt, max=ymx)
  ymnx = [ymn-0.1, ymx+0.1]

  clr = getcolor(/load)
  plot, [xmn,xmx], ymnx, /nodata, background=clr.white, $ 
    color=clr.black, xthick=3.0, ythick=3.0, xtitle='!17 log!d10!n(f!7k!X)', $
    ytitle='log!d10!n(W/!7k!X)', xstyle=1, ystyle=1, charsize=1.9
  oploterror, xplt, yplt, ysig, psym=1, color=clr.blue, ERRCOLOR=clr.blue
  
  nplt = 100L
  cog = dblarr(nplt)

  cog_b1 = 1./b1gd^2 ;; Inverse to speed up calculation
  cog_b2 = 1./b2gd^2 ;; Inverse to speed up calculation
  cog_N1 = 10^N1gd / (b1gd*1.e5)
  cog_N2 = 10^N2gd / (b2gd*1.e5)

  ;; Good value
  xplt2 = (xmn + (xmx-xmn)*findgen(nplt)/float(nplt))
  for kk=0L, nplt-1 do begin
      tcog_const = -1.*1.497e-2*(10^xplt2[kk])
      cog[kk] = fuse_twocog_cog() 
  endfor
  oplot, xplt2, alog10(cog), color=clr.red

  ;; Error edges (high N, high b)
;  Ngd = Ngd + 2*sigN
;  bgd = bgd + 2*sigb
;  cos_bval=bgd
;  tau = 1.497e-2*(10^xplt2)*(10^Ngd)/(bgd*1e5)
;  for q=0L,nplt-1 do cog[q] = x_calccog_cog(tau[q])
;  oplot, xplt2, alog10(cog), color=clr.green, linestyle=2


;; Error edges (high N, high b)
;  Ngd = Ngd - 4*sigN
;  bgd = bgd - 4*sigb
;  cos_bval=bgd
;  tau = 1.497e-2*(10^xplt2)*(10^Ngd)/(bgd*1e5)
;  for q=0L,nplt-1 do cog[q] = x_calccog_cog(tau[q])
;  oplot, xplt2, alog10(cog), color=clr.green, linestyle=2

  if keyword_set( PSFILE ) then begin
      ps_close, /noprint, /noid
      device, decompose=1
      !p.thick = 1
      !p.charthick = 1
  endif

  return
end
