;+ 
; NAME:
; lls_ioncorr
;  V1.2
;
; PURPOSE:
;    Given a list of LLS .dat files each with a Cloudy ioncorr file
;    and a U value, calculate the IC and apply to the column densities
;    to derive [Fe/H] and [alpha/H]
;
; CALLING SEQUENCE:
;   lls_ioncorr, list
;
; INPUTS:
;  list -- ASCII file which is a list of LLS .dat files.
;
; RETURNS:
;
; OUTPUTS:
;  Series of DLA files
;
; OPTIONAL KEYWORDS:
;  FLAG=  0=Default; 1=GRB; 2=DLA
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   lls_ioncorr, list_fil
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   20-Jul-2008 Written by JXP 
;-
;------------------------------------------------------------------------------
pro lls_ioncorr, list, FLG_PLT=flg_plt, FLG_CII=flg_cii, FLAG=flag, $
                 PRINTONLY=printonly, ROOT=root, FILE=file, NOWRITE=nowrite

; parse_llslst -- Reads in DLA data to a structure

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
          'lls_ioncorr, list, /FLG_PLT, /FLG_CII, /printonly, FLAG=, /FILE ' + $
          '/NOWRITE [v1.1]'
    return
  endif 

  if not keyword_set(FLAG) then flag = 0
  if not keyword_set(ROOT) then root = ''

  ;; Read the list
  readcol, list, list_nm, format='A'
  nlls = n_elements(list_nm)

  ;; Loop
  for nn=0L,nlls-1 do begin
      print, 'lls_updabnd: Analysing ', list_nm[nn]

      ;; Read the column densities
      lls_struct, lls, list_nm[nn], /FILE, ROOT=root, /ION, /FINE, /VPFIT

      ;; Ion Corrections
      for sys=0L,lls[nn].nsys-1 do begin
          ;; Full calc
          lls_allabnd, lls, nn, sys, FLG_PLT=flg_plt, ROOT=root

          ;; Ion List
          len = strlen(lls[nn].systems[sys].tab_fil)
          ION_fil = root+strmid(lls[nn].systems[sys].tab_fil,0,len-3)+'ion'
          close, 13
          openw, 13, ION_fil
          imx = 99 < ((size(lls.systems[sys].ion))[1]-1)
          for ww=1,imx do begin
              for qq=1,6 do begin
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
                      endif
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
          ;; Output Element List
          ;;
          XH_fil = root+strmid(lls[nn].systems[sys].tab_fil,0,len-3)+'XH'
          openw, 13, XH_fil
          imx = 99 < ((size(lls.systems[sys].elm))[1]-1)
          for ww=1,imx do begin
              if lls[nn].systems[sys].elm[ww].flgclm NE 0 then begin
                  x_logclm, lls[nn].systems[sys].elm[ww].clm, $
                            lls[nn].systems[sys].elm[ww].sigclm, ans, sig
                  if SIG GT 9.9 then sig = 9.9
                  print, ww, ans, sig, lls[nn].systems[sys].elm[ww].flgclm, $
                         lls[nn].systems[sys].elm[ww].flginst, $
                         format='(i2,1x,f6.3,1x,f6.3,1x,i1,1x,i3)'
                  printf, 13, ww, ans, sig, lls[nn].systems[sys].elm[ww].flgclm, $
                          lls[nn].systems[sys].elm[ww].flginst, $
                          format='(i2,1x,f6.3,1x,f6.3,1x,i1,1x,i3)'
              endif 
          endfor
          close,13

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
                    GT 1 then begin
                      print, 'lls_updabnd: Taking the AODM limit, not VPFIT', $
                             Zgd[ss], igd[ss]
                  endif else begin
                      lls[nn].systems[sys].ion[Zgd[ss]].state[igd[ss],0].flgclm = 1
                      lls[nn].systems[sys].ion[Zgd[ss]].state[igd[ss],0].clm = $
                        vpion[Zgd[ss]].state[igd[ss],0].clm
                      lls[nn].systems[sys].ion[Zgd[ss]].state[igd[ss],0].sigclm = $
                        vpion[Zgd[ss]].state[igd[ss],0].sigclm
                  endelse
              endfor
          endif
      
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
          ;; Alpha
;          dla_calcalpha, lls.systems, sys
      endfor
  endfor
      
  ;; Write
  if not keyword_set(NOWRITE) then lls_writestr, lls;, flag=flag
  print, 'lls_updabd: All done'

  return
end
