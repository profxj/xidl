;+ 
; NAME:
; fuse_cog
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
;   11-Sep-2003 Written by JXP
;-
;------------------------------------------------------------------------------
pro fuse_cog, strct_fil, cog_fil, Nlmt, blmt, CHICHK=chichk, PLTONLY=pltonly, $
              NSTP=nstp, BSTP=bstp, PSFILE=psfile, OUTFIL=outfil, EXACT=exact

  common x_calccog_cmm, cog_bval, cog_tau, cog_strct

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'fuse_cog, strct, cog_fil (v1.0)' 
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
  redew = strct.EW[0] / (strct.zabs+1.) / strct.wrest / 1000.
  redsigew = strct.sigEW[0] / (strct.zabs+1.) / strct.wrest / 1000.
  
  ;; Grab fvalues
  getfnam, strct.wrest, fval
  flambda = strct.wrest*fval
  srt = sort(flambda)
  printcol, strct[srt].wrest, alog10(flambda[srt])
  
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
      for i=0L, 2*NSTP/3 do begin
          sub_prob = total( prob_arr[ii-i:ii+i,*] )
          if sub_prob GT 0.683*tot_prob then break
      endfor
      if i EQ 0 then stop
      sigN = (Nlmt[1]-Nlmt[0])*float(i)/float(nstp)
      
      ;; Error in b
      for j=0L, 2*BSTP/3 do begin
          sub_prob = total( prob_arr[*,(jj-j)>0:jj+j] )
          if sub_prob GT 0.683*tot_prob then break
      endfor
      if j EQ 0 then stop
      sigb = (blmt[1]-blmt[0])*float(j)/float(bstp)
      
      print, 'fuse_cog: N = ', Ngd, '+/-', sigN
      print, 'fuse_cog: b = ', bgd, '+/-', sigb
      ;; Write to file
      if keyword_set( OUTFIL ) then begin
          openw, 22, outfil
          printf, 22, Ngd, sigN
          printf, 22, bgd, sigb
          close, 22
      endif
  endif else begin
      Ngd = pltonly[0]
      bgd = pltonly[1]
      sigN = pltonly[2]
      sigb = pltonly[3]
  endelse

  cog_bval=bgd

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
    ytitle='log!d10!n(W/!7k!X)', xstyle=1, ystyle=1, charsize=2.2, $
    xticks=5, xmargin=[2,1], ymargin=[2,0]
  oploterror, xplt, yplt, ysig, psym=1, color=clr.blue, ERRCOLOR=clr.blue
  
  nplt = 100L
  cog = dblarr(nplt)

  ;; Good value
  xplt2 = (xmn + (xmx-xmn)*findgen(nplt)/float(nplt))
  tau = 1.497e-2*(10^xplt2)*(10^Ngd)/(bgd*1e5)
  for q=0L,nplt-1 do cog[q] = x_calccog_cog(tau[q])
  oplot, xplt2, alog10(cog), color=clr.red

  ;; Error edges (high N, high b)
  Ngd = Ngd + sigN
  bgd = bgd + sigb
  cog_bval=bgd
  tau = 1.497e-2*(10^xplt2)*(10^Ngd)/(bgd*1e5)
  for q=0L,nplt-1 do cog[q] = x_calccog_cog(tau[q])
  oplot, xplt2, alog10(cog), color=clr.green, linestyle=2

  ;; Error edges (high N, high b)
  Ngd = Ngd - 2*sigN
  bgd = bgd - 2*sigb
  cog_bval=bgd
  tau = 1.497e-2*(10^xplt2)*(10^Ngd)/(bgd*1e5)
  for q=0L,nplt-1 do cog[q] = x_calccog_cog(tau[q])
  oplot, xplt2, alog10(cog), color=clr.green, linestyle=2

  if keyword_set( PSFILE ) then begin
      ps_close, /noprint, /noid
      device, decompose=1
      !p.thick = 1
      !p.charthick = 1
  endif

  return
end
