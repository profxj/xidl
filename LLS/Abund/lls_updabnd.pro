;+ 
; NAME:
; lls_updabnd
;  V1.2
;
; PURPOSE:
;    Given a list of LLS .dat files, calculate the abundances, modify
;    the .dat files and write everything out.  This is the main driver
;    of AODM analysis for LLS.
;
; CALLING SEQUENCE:
;   lls_updabd, dla, /FLG_PLT, /FLG_CII, /FIXZN, /FINE
;
; INPUTS:
;  list -- List of LLS .dat files.
;
; RETURNS:
;
; OUTPUTS:
;  Series of DLA files
;
; OPTIONAL KEYWORDS:
;  FLAG=  0=Default; 1=GRB; 2=DLA
;  /NOXH  --  Dont calcualte X/H values 
;  /SKIP -- Skip systems without a datafile provided
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   lls_updabd, list_fil
;   lls_updabd, 'blah', LLS=lls
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   01-Oct-2004 Written by JXP 
;-
;------------------------------------------------------------------------------
pro lls_updabnd, list, FLG_PLT=flg_plt, FLG_CII=flg_cii, FLAG=flag, $
                 PRINTONLY=printonly, ROOT=root, FILE=file, NOWRITE=nowrite, $
                 NOXH=noxh, LLS=lls, SKIP=skip, CORRO=corro

; parse_llslst -- Reads in DLA data to a structure

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
          'lls_updabd, list, /FLG_PLT, /FLG_CII, /printonly, FLAG=, /FILE ' + $
          '/NOWRITE, /NOXH, /SKIP [v1.1]'
    return
  endif 

  if not keyword_set(FLAG) then flag = 0
  if not keyword_set(ROOT) then root = ''
  sv_cldy_fil = ''

  ;; Parse
  if not keyword_set(PRINTONLY) then begin
     if not keyword_set(LLS) then $
        lls_struct, lls, list, ROOT=root, FILE=file, /NOFILLELM ;, FLAG=flag
      nlls = n_elements(lls)
  endif else nlls = 1L

  ;; Loop
  for nn=0L,nlls-1 do begin
      print, 'lls_updabnd: Analysing ', lls[nn].qso, lls[nn].zabs

      ;; Abundances
      for sys=0L,lls[nn].nsys-1 do begin
          ;; Full calc
          lls_allabnd, lls, nn, sys, FLG_PLT=flg_plt, ROOT=root, SKIP=skip

          ;; Ion List
          len = strlen(lls[nn].systems[sys].tab_fil)
          ION_fil = root+strmid(lls[nn].systems[sys].tab_fil,0,len-3)+'ion'
          close, 13
          openw, 13, ION_fil
          imx = 99 < ((size(lls.systems[sys].ion))[1]-1)
          for ww=2,imx do begin ;; Now excluding Hydrogen!  Should have all along (JXP: 5 Jan 2015)
              for qq=1,9 do begin
                  tidx = lls[nn].systems[sys].ion[ww].indx[qq]
                  for k=1,tidx do begin
                      if lls[nn].systems[sys].ion[ww].state[qq,k].flgclm NE -1 then begin
                          if ww NE 1  then begin
                              x_logclm, lls[nn].systems[sys].ion[ww].state[qq,k].clm, $
                                        lls[nn].systems[sys].ion[ww].state[qq,k].sigclm, ans, sig
                          endif else begin
                              ans =  lls[nn].systems[sys].ion[ww].state[qq,k].clm
                              sig =  lls[nn].systems[sys].ion[ww].state[qq,k].sigclm
                          endelse
                      endif else begin
                        ans = 0.
                        sig = 9.9
                      endelse      
                      if sig GT 9.9 then sig = 9.9
                      printf, 13,  $
                              lls[nn].systems[sys].ion[ww].state[qq,k].lambda, ans,sig, $
                              lls[nn].systems[sys].ion[ww].state[qq,k].flgclm, $
                              lls[nn].systems[sys].ion[ww].state[qq,k].flginst, $
                              format='(f10.4,1x,2f8.4,1x,i2,1x,i3)'
                  endfor
              endfor
          endfor
          close, 13

          ;;cccccccccccccccccccccccccccccc
          ;; VPFIT?
          ;;
          if strlen(strtrim(lls[nn].systems[sys].vpfil,2)) GT 0 then begin
              ism_sumvpfit, lls[nn].systems[sys].vpfil, vpion
              ;; Overwrite? (No limits)
              gd = where(vpion.state[*,0].sigclm GT 0., ngd)
              szvp = size(vpion.state[*,0],/dimensions)
              Zgd = gd/szvp[0]
              igd = gd mod szvp[0]
              for ss=0L,ngd-1 do begin
                  if lls[nn].systems[sys].ion[Zgd[ss]].state[igd[ss],0].flgclm $
                    GT 1 then $
                    print, 'lls_updabnd: Warning, the AODM says ' + $
                           'it is a limit, but we will take VPFIT.', $
                           Zgd[ss], igd[ss]
                  ;; Ions
                  lls[nn].systems[sys].ion[Zgd[ss]].state[igd[ss],0].flgclm = 1
                  lls[nn].systems[sys].ion[Zgd[ss]].state[igd[ss],0].clm = $
                    vpion[Zgd[ss]].state[igd[ss],0].clm
                  lls[nn].systems[sys].ion[Zgd[ss]].state[igd[ss],0].sigclm = $
                    vpion[Zgd[ss]].state[igd[ss],0].sigclm
                  ;; Elm
                  if Zgd[ss] NE 1 and igd[ss] NE 1 then begin
                      ;; flag needs to be odd in the .clm file 
                      mtz = where((lls[nn].systems[sys].ion[Zgd[ss]].state[*,1:*].flgclm $
                                   MOD 2)  GT 0, nmtz)
                      mti = where((lls[nn].systems[sys].ion[Zgd[ss]].state[igd[ss],1:*].flgclm $
                                   MOD 2)  GT 0, nmti)
                      if nmti NE 0 or (nmti + nmtz) EQ 0 then begin
                          lls[nn].systems[sys].elm[Zgd[ss]].flgclm = 1
                          lls[nn].systems[sys].elm[Zgd[ss]].clm = $
                            vpion[Zgd[ss]].state[igd[ss],0].clm
                          lls[nn].systems[sys].elm[Zgd[ss]].sigclm = $
                            vpion[Zgd[ss]].state[igd[ss],0].sigclm
                      endif
                  endif
              endfor
          endif
      

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Calculate XH using Cloudy corrections
          gdelm = where(lls[nn].systems[sys].elm.flgclm NE 0, ngdelm)
          if ngdelm NE 0 and not keyword_set(NOXH) then begin
              flg_IC = 0
              ;; Cloudy?
              fil = findfile(strtrim(lls[nn].cldyfil,2)+'*',count=nfil)
              if nfil EQ 0 then begin
                  print, 'lls_updabnd:  No cloudy file.  Assuming no IC'
              endif else begin
                  print, 'lls_updabnd:  Using cldy file ', fil[0]
                  ;; U
                  Ubest = lls[nn].systems[sys].U $
                          + [0., [-1,1]*lls[nn].systems[sys].Usig]
                  ;;
                  flg_IC = 1
                  ;; Read Cloudy (avoid reading twice)
                  if not keyword_set(cldy) or (not strmatch(sv_cldy_fil, strtrim(fil[0],2))) then $
                     cldy = xmrdfits(fil[0], 1, /silent)
                  sv_cldy_fil = strtrim(fil[0],2)
                  ;; Parse
                  mn = min(abs(cldy.NHI - lls[nn].systems[sys].NHI),imn)
                  NHIc = cldy[imn].NHI
                  nHc = -1.  ; Default
                  mn = min(abs(cldy.FeH - lls[nn].systems[sys].alphaH),imn)
                  FeHc = cldy[imn].FeH
                  gd = where(cldy.NHI EQ NHIc AND cldy.FeH EQ FeHc $
                             AND cldy.nH EQ nHc, ngd)
                  scldy = cldy[gd]
                  srt = sort(scldy.U)
                  scldy = scldy[srt]
              endelse
              ;; Calculate X/H with IC (as available)
              for tt=0L,ngdelm-1 do begin
                  ;; Z and ion
                  Zval = gdelm[tt]
                  mt = where((lls[nn].systems[sys].ion[Zval].state[*,1:*].flgclm $
                              MOD 2)  GT 0, nmt)
                  sz = size(lls[nn].systems[sys].ion[Zval].state[*,1:*].flgclm,/dim)
                  mt = mt MOD sz[0]
                  mt = mt[uniq(mt, sort(mt))]
                  nmt = n_elements(mt)
                  
                  if nmt NE 1 then begin
                     print, 'Zval = ', Zval, mt
                     stop ;; Cannot have 2 ions chosen for 1 element
                  endif
                  ival = mt[0]
;                  print, 'Ion: ', Zval, ival

                  ;; Flag and instrument
                  lls[nn].systems[sys].XH[Zval].flgclm = $
                    lls[nn].systems[sys].elm[Zval].flgclm 
                  lls[nn].systems[sys].XH[Zval].flginst = $
                    lls[nn].systems[sys].elm[Zval].flginst 
                  ;; Value and error
                  x_logclm, lls[nn].systems[sys].elm[Zval].clm, $
                            lls[nn].systems[sys].elm[Zval].sigclm, ans, sig
                  getabnd, nm, Zval, abnd, flag=1
                  case flg_IC of 
                      0: begin ; No IC -- Assuming low-ions 
                          lls[nn].systems[sys].XH[Zval].clm = $
                            ans - lls[nn].systems[sys].NHI + 12.0 - abnd
                          ;;MF. Correct sigNHI in NHIsig
                          ;;if lls[nn].systems[sys].sigNHI[1] GT 0. then $
                          ;;  sigNHI = mean(lls[nn].systems[sys].sigNHI) else $
                          ;;  sigNHI = lls[nn].systems[sys].sigNHI[0]
                          if lls[nn].systems[sys].NHIsig[1] GT 0. then $
                             sigNHI = mean(lls[nn].systems[sys].NHIsig) else $
                                sigNHI = lls[nn].systems[sys].NHIsig[0]
                          lls[nn].systems[sys].XH[Zval].sigclm = $
                             sqrt([sig,sig]^2 + sigNHI^2)
;                             sqrt(sig^2 + sigNHI^2)
                       end
                      1: begin  ; IC
                          pmod = scldy.X[Zval,ival]-scldy.X[1,1]
                          IC = interpol(pmod, scldy.U, Ubest)
                          ;; OI
                          if Zval EQ 8 and ival EQ 1 and $
                            NOT keyword_set(CORRO) then IC = [0., 0.2, 0.2]
                          ;; Value
                          lls[nn].systems[sys].XH[Zval].clm = $
                            ans - lls[nn].systems[sys].NHI + 12.0 - abnd - IC[0]
                          ;; Error
                          sigIC = abs(IC[1:2]-IC[0])
                          lls[nn].systems[sys].XH[Zval].sigclm = $
                            sqrt([sig,sig]^2 + lls[nn].systems[sys].NHIsig^2 $
                                + sigIC^2)
;                          if nn EQ 1 and Zval EQ 6 then stop
                      end
                      else: stop
                  endcase
              endfor
              ;; Set NH and x too
              if flg_IC GT 0 then begin
                  pmod = scldy.X[1,1]
                  IC = interpol(pmod, scldy.U, Ubest)
                  ;; NH
                  lls[nn].systems[sys].NH = $
                    lls[nn].systems[sys].NHI - IC[0]
;                  sigIC = total( abs(IC[1:2]-IC[0]) ) /2. 
                  sigIC = abs(IC[1:2]-IC[0])
                  if lls[nn].systems[sys].NHIsig[1] LE 1e-5 then $
                    lls[nn].systems[sys].NHIsig[1] = lls[nn].systems[sys].NHIsig[0] 
                  lls[nn].systems[sys].NHsig = $
                    sqrt(lls[nn].systems[sys].NHIsig^2 + sigIC^2)
                  ;; x
                  lls[nn].systems[sys].logx = IC[0]
                  lls[nn].systems[sys].sig_logx = abs(IC[1:2]-IC[0])
              endif
          endif
                  
          
          ;;cccccccccccccccccccccccccccccc
          ;; Output XH List
          ;;
          XH_fil = root+strmid(lls[nn].systems[sys].tab_fil,0,len-3)+'XH'
          openw, 13, XH_fil
          imx = 99 < ((size(lls.systems[sys].XH))[1]-1)
          for ww=1,imx do begin
              if lls[nn].systems[sys].XH[ww].flgclm NE 0 then begin
                  ans = lls[nn].systems[sys].XH[ww].clm
                  sig = lls[nn].systems[sys].XH[ww].sigclm
;                  x_logclm, lls[nn].systems[sys].XH[ww].clm, $
;                            lls[nn].systems[sys].XH[ww].sigclm, ans, sig
                  ;if SIG GT 9.9 then sig = 9.9
                  sig[where(sig GT 9.9,/null)] = 9.9
                  print, ww, ans, sig, lls[nn].systems[sys].XH[ww].flgclm, $
                         lls[nn].systems[sys].XH[ww].flginst, $
                         format='(i2,1x,f6.3,1x,f6.3,1x,f6.3,1x,i1,1x,i3)'
                  printf, 13, ww, ans, sig, lls[nn].systems[sys].XH[ww].flgclm, $
                          lls[nn].systems[sys].XH[ww].flginst, $
                          format='(i2,1x,f6.3,1x,f6.3,1x,f6.3,1x,i1,1x,i3)'
              endif 
          endfor
          close,13
          
          if not keyword_set(NOXH) then begin
             ;; Alpha
             lls_calcalpha, lls, nn, sys
             
             ;; Fe/H
             lls_calcfeh, lls, nn, sys
          endif
          ;stop

          ;;cccccccccccccccccccccccccccccc
          ;; Output Total List
          ;;
          XH_fil = root+strmid(lls[nn].systems[sys].tab_fil,0,len-3)+'all'
          openw, 13, XH_fil
          imx = 99 < ((size(lls.systems[sys].ion))[1]-1)
          for ww=1,imx do begin
              for qq=1,6 do begin
                  if lls[nn].systems[sys].ion[ww].state[qq].flgclm NE 0 then begin
                      x_logclm, lls[nn].systems[sys].ion[ww].state[qq].clm, $
                                lls[nn].systems[sys].ion[ww].state[qq].sigclm, $
                                ans, sig
                      if sig GT 9.9 then sig = 9.9
                      printf, 13,  ww, qq, ans, sig, $
                              lls[nn].systems[sys].ion[ww].state[qq].flgclm, $
                              lls[nn].systems[sys].ion[ww].state[qq].flginst, $
                              format='(i2,1x,i2,1x,f6.3,1x,f6.3,1x,i1,1x,i3)'
                  endif
              endfor
          endfor
          close,13

          ;;cccccccccccccccccccccccccccccc
          ;; DH
          flag = 1
          case flag of
              1: 
              else: begin
                  lls[nn].flg_DH = lls[nn].elm[99].flgclm 
                  if lls[nn].elm[99].flgclm NE 0 then begin
                      lls[nn].DH = lls[nn].elm[99].clm / 10^lls[nn].NHI * 1e5
                      lls[nn].sigDH = lls[nn].elm[99].sigclm / 10^lls[nn].NHI * 1e5
                  endif
              end
          endcase

          ;; Fe
;          dla_calcfeh, lls.systems, sys

          ;; CII* 
;          if keyword_set(flg_CII) then begin
;              if lls[nn].systems[sys].ion[6].state[6].flgclm GE 0 then begin
;                  lls[nn].flgCII = lls[nn].systems[sys].ion[6].state[6].flgclm + 1
;                  x_logclm, lls[nn].systems[sys].ion[6].state[6].clm, $
;                            lls[nn].systems[sys].ion[6].state[6].sigclm, ans, sig
;                  lls[nn].systems[sys].CII = ans
;                  lls[nn].systems[sys].sigCII = sig
;              endif else lls[nn].systems[sys].flgCII = 0
;          endif
      endfor

      ;; Update totals
      if lls[nn].nsys GT 0 then begin
          nsys = lls[nn].nsys
          ;; NHI, NH
          ;; Turned off by JXP On 09 Dec 2014
          ;NHI = total(10.^[lls[nn].systems[0:nsys-1].NHI])
          ;sigv = lls[nn].systems[0:nsys-1].NHIsig[0] * alog(10.) * $
          ;       10.d^lls[nn].systems[0:nsys-1].NHI
          ;sigv = sqrt(total(sigv^2))
          ;x_logclm, NHI, sigv, ans, sig
          ;lls[nn].NHI = ans
          ;lls[nn].sigNHI = sig
          
          ;; NH
          NH = total(10.^[lls[nn].systems[0:nsys-1].NH])
          sigv = lls[nn].systems[0:nsys-1].NHsig[0] * alog(10.) * $
                 10.d^lls[nn].systems[0:nsys-1].NH
          sigv = sqrt(total(sigv^2))
          x_logclm, NH, sigv, ans, sig
          lls[nn].NH = ans
          lls[nn].NHsig = sig

          ;; [M/H]
          ;; Marie's version (1/9/2015)
          good=where(lls[nn].systems[0:nsys-1].flg_alpha eq 1 or $
                     lls[nn].systems[0:nsys-1].flg_alpha eq 2,ngd)
          lower=where(lls[nn].systems[0:nsys-1].flg_alpha eq 2,nlo)
          upper=where(lls[nn].systems[0:nsys-1].flg_alpha eq 3,nup)	
          if ngd ne 0 then begin
             mtlsgd=10.^(lls[nn].systems[good].NH+lls[nn].systems[good].alphaH)
             if nlo gt 0 then lls[nn].flg_MH=2 else lls[nn].flg_MH=1
             lls[nn].MHave=alog10(total(mtlsgd))-lls[nn].NH
             ;; Error (crude)
             mx=max(mtlsgd,imx)
             sig=sqrt(lls[nn].systems[good[imx]].NHsig^2+$
                      lls[nn].systems[good[imx]].sig_alphaH^2)
             lls[nn].MHsig=[mean(sig[0,*]),mean(sig[1,*])]
;          endif else lls[nn].flg_MH=3
          endif else begin
             mtlsup=10.^(lls[nn].systems[upper].NH+lls[nn].systems[upper].alphaH)
             lls[nn].MHave=alog10(total(mtlsup))-lls[nn].NH
             lls[nn].flg_MH=3
          endelse
          
          ;; Previous code
          ;gda = where(lls[nn].systems[0:nsys-1].flg_alpha EQ 1 OR $
          ;           lls[nn].systems[0:nsys-1].flg_alpha EQ 5 OR $
          ;           lls[nn].systems[0:nsys-1].flg_alpha EQ 4, ngda)
          ;if ngda NE 0 then begin
          ;    lls[nn].flg_MH = 1
          ;    ;; Value
          ;    mtls = 10.^(lls[nn].systems[gda].NH + lls[nn].systems[gda].alphaH)
          ;    lls[nn].MHave = alog10( total(mtls) ) - lls[nn].NH
          ;    ;; Error (crude)
          ;    mx = max(mtls, imx)
          ;    sig = sqrt( lls[nn].systems[gda[imx]].NHsig^2 + $
          ;                lls[nn].systems[gda[imx]].sig_alphaH^2)
          ;    lls[nn].MHsig = mean(sig)
          ;endif
      endif
  endfor
      
  ;; Write
  if not keyword_set(NOWRITE) then lls_writestr, lls;, flag=flag
  print, 'lls_updabd: All done'

  return
end
