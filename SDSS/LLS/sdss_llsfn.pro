;+ 
; NAME:
; sdss_llsfn
;
; PURPOSE:
;    Given a LLS struct and gz list, calculate f(N)
;
; CALLING SEQUENCE:
;   sdss_llsfn
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
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Nov-2004 Written by JXP
;-
pro sdss_llsfn, llslst, PSFILE=psfile, NMIN=nmin, $
                GZFIL=gzfil, DLALST=dlalst, ZBIN=zbin, $
                STRCT=strct, GZSTR=gzstr, $
                NBIN=nbin, stp=stp, SUBR=subr, NODBL=nodbl, OUTFIL=outfil

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'sdss_llsfn, LIST, STRCT=strct  [v1.0]'
    return
  endif 
  ;; Get structure if necessary

  if not keyword_set( CSZ ) then csz = 2.3
  if not keyword_set( NBIN ) then nbin = 4L
  if not keyword_set( STP ) then stp = 0.3
  if not keyword_set( NMIN ) then nmin = 19.
  if not keyword_set( NMAX ) then nmax = 20.3
  if not keyword_set( ZBIN ) then zbin = [1.7, 5.5]
  if not keyword_set( YRNG ) then yrng = [-22., -19.5]
  if not keyword_set( XRNG ) then xrng = [18.8, 20.4]
  if not keyword_set( YMRG ) then ymrg = [4,1]
  if not keyword_set(LBLSZ) then lblsz = 1.4
  
  compile_opt strictarr

  ;; 3 cases
  if not keyword_set(GZSTR) then gzstr = xmrdfits(gzfil, 1, /silent)

  ;; Read LLS
  for qq=0L,n_elements(llslst)-1 do begin
      lls_struct, tmp, llslst[qq], /noelm, ROOT=getenv('LLSTREE')+'/'
      if qq EQ 0 then all_lls = tmp else all_lls = [all_lls,tmp]
  endfor

  ;; SDSS
  i1 = sdss_llsstat(all_lls, gzstr) 

  i2 = where(all_lls[i1].zabs GE zbin[0] AND $
             all_lls[i1].zabs LE zbin[1] AND $
             all_lls[i1].NHI GE NMIN AND $
             all_lls[i1].NHI LT 20.3, nlls)
  if nlls NE 0 then begin
      indx = i1[i2]
      all_lls = all_lls[indx]
  endif else delvarx, all_lls

  zmean = mean(all_lls.zabs)
  print, 'Mean z = ', zmean

  ;; dX
  dX = sdss_llsdx( zbin[0], zbin[1], GZSTR=gzstr)
  print, 'dX = ', dX

  ;; f(N) plot
  bins = nmin + findgen(nbin)*stp
  yplt = dblarr(nbin)
  svex = dblarr(nbin)
  yerr1 = dblarr(nbin) 
  yerr2 = dblarr(nbin)
  xplt = dblarr(nbin)
  o3xplt = dblarr(nbin)
  xval = dblarr(nbin)
  x1 = dblarr(nbin)
  x2 = dblarr(nbin)
  svdn = dblarr(nbin)

  for qq=0L,nbin-1 do begin
      dex = abs(bins[qq]+stp-20.3)
      if dex LT stp then ex = dex else ex = 0.
      svex[qq] = ex
      xval[qq] = bins[qq] + (stp+ex)/2.
      gd = where(all_lls.NHI GE bins[qq] AND $
                 all_lls.NHI LT bins[qq]+stp+ex,ngd)
      dN = 10^(bins[qq]+stp+ex) - 10^bins[qq]
      svdN[qq] = dN
      case ngd of
          0: xplt[qq] = xval[qq]
          1: xplt[qq] = all_lls[gd].NHI
          else: begin
              xplt[qq] = alog10(mean(10^all_lls[gd].NHI))
              o3xplt[qq] = mean(all_lls[gd].NHI)
          end
      endcase
      if ngd NE 0 then begin
          yplt[qq] = alog10(ngd/dN/dX) 
          print, qq, ngd, bins[qq], dN, dX, yplt[qq]
          ;; Poisson 68.3% (1sigma)
          val = x_poisscl(double(ngd), sigma=1)
          ;; Fill
          yerr1[qq] = alog10(val[0]/dN/dX) 
          yerr2[qq] = alog10(val[1]/dN/dX) 
      endif else begin
          yplt[qq] = -99.
          val = x_poisscl(0., sigma=2)
          yerr1[qq] = alog10(val[0]/dN/dX) 
      endelse
  endfor

  ;; Output
  if keyword_set(OUTFIL) then $
    writecol, outfil, xplt, yplt, yerr1, yerr2

  ;; Functions
  ;; Single power law
  maxb = x_maxsngpow(10^all_lls.NHI, NMIN=10.^NMIN, SMM=smmnh, $
                    KSPROB=ksprob, NOISE=0.05, ERR=errs, CL=0.683, $
                    ARNG=[-5., -0.], NMAX=10^NMAX)
  k1 = -1.*float(nlls/dX) * (maxb+1.) / $
       ((10.^Nmin)^(1.+maxb) - (10.^Nmax)^(1.+maxb))
  sigk1 = x_poisscl(double(nlls), 0.8413)*k1/float(nlls)

  print, 'Maxb = ', maxb, errs
  print, 'k1 = ', k1
  if keyword_set(KSPROB) then print, 'Single KS = ', ksprob

  ;; Structure
  strct = { $
          xplt: xplt, $
          o3xplt: o3xplt, $
          xval: xval, $
          fN: yplt, $
          yerr1: yerr1, $
          yerr2: yerr2, $
          dN: svdN, $
          bins: bins, $
          sigfn1: yerr1, $
          sigfn2: yerr2, $
          svex: svex, $
          nlls: nlls, $
          zmean: zmean, $
          dX: dX, $
          k1: k1, $
          sigk1: sigk1, $
          a1: maxb, $
          siga1: errs $
          }
;            k2: k2, $
;            a2: val[1], $
;            Ng: val[0], $
;            gerr:gerr, $
;            k3: k3, $
;            Nd: vald[0], $
;            a3: vald[1], $
;            a4: vald[2], $
;           derr:derr $

  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  
  print, 'fig_fnall: All done!'

  return

end
