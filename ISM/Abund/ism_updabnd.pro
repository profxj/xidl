;+ 
; NAME:
; ism_updabnd
;  V1.2
;
; PURPOSE:
;    Given a list of ISM .dat files, calculate the abundances, modifty
;    the .dat files and write everything out.  This is the main driver
;    of AODM analysis for ISM abundances.
;
; CALLING SEQUENCE:
;   ism_updabd, ism
;
; INPUTS:
;  list -- List of DLA .dat files.
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
;   ism_updabd, ism
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   01-Oct-2004 Written by JXP 
;-
;------------------------------------------------------------------------------
pro ism_updabnd, list, FLG_PLT=flg_plt, FLG_CII=flg_cii, FLAG=flag

; parse_ismlst -- Reads in DLA data to a structure

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'ism_updabd, list, /FLG_PLT, /FLG_CII, FLAG= [v1.1]'
    return
  endif 

  if not keyword_set(FLAG) then flag = 0

  ;; Parse
  parse_ismlst, ism, list, /noelm, FLAG=flag
  nism = n_elements(ism)

  ;; Loop
  for nn=0L,nism-1 do begin
      ;; Abundances
      ism_allabnd, ism, nn, FLG_PLT=flg_plt

      ;; Ion List
      len = strlen(ism[nn].tab_fil)
      ION_fil = strmid(ism[nn].tab_fil,0,len-3)+'ion'
      close, 13
      openw, 13, ION_fil
      imx = 99 < ((size(ism.ion))[1]-1)
      for ww=1,imx do begin
          for qq=1,6 do begin
              tidx = ism[nn].ion[ww].indx[qq]
              for k=1,tidx do begin
                  if ism[nn].ion[ww].state[qq,k].flgclm NE -1 then begin
                      if ww NE 1  then begin
                          x_logclm, ism[nn].ion[ww].state[qq,k].clm, $
                            ism[nn].ion[ww].state[qq,k].sigclm, ans, sig
                      endif else begin
                          ans =  ism[nn].ion[ww].state[qq,k].clm
                          sig =  ism[nn].ion[ww].state[qq,k].sigclm
                      endelse
                  endif
                  if sig GT 9.9 then sig = 9.9
                  printf, 13,  $
                    ism[nn].ion[ww].state[qq,k].lambda, ans,sig, $
                    ism[nn].ion[ww].state[qq,k].flgclm, $
                    ism[nn].ion[ww].state[qq,k].flginst, $
                    format='(f10.4,1x,2f8.4,1x,i2,1x,i3)'
              endfor
          endfor
      endfor
      close, 13

      ;;cccccccccccccccccccccccccccccc
      ;; Output Element List
      ;;
      XH_fil = strmid(ism[nn].tab_fil,0,len-3)+'XH'
      openw, 13, XH_fil
      imx = 99 < ((size(ism.elm))[1]-1)
      for ww=1,imx do begin
          if ism[nn].elm[ww].flgclm NE 0 then begin
              x_logclm, ism[nn].elm[ww].clm, ism[nn].elm[ww].sigclm, ans, sig
              if SIG GT 9.9 then sig = 9.9
              print, ww, ans, sig, ism[nn].elm[ww].flgclm, $
                ism[nn].elm[ww].flginst, $
                format='(i2,1x,f6.3,1x,f6.3,1x,i1,1x,i3)'
              printf, 13, ww, ans, sig, ism[nn].elm[ww].flgclm, $
                ism[nn].elm[ww].flginst, $
                format='(i2,1x,f6.3,1x,f6.3,1x,i1,1x,i3)'
          endif 
      endfor
      close,13
      
      ;;cccccccccccccccccccccccccccccc
      ;; Output Total List
      ;;
      XH_fil = strmid(ism[nn].tab_fil,0,len-3)+'all'
      openw, 13, XH_fil
      imx = 99 < ((size(ism.ion))[1]-1)
      for ww=1,imx do begin
          for qq=1,6 do begin
              if ism[nn].ion[ww].state[qq].flgclm NE 0 then begin
                  x_logclm, ism[nn].ion[ww].state[qq].clm, $
                    ism[nn].ion[ww].state[qq].sigclm, ans, sig
                  if sig GT 9.9 then sig = 9.9
                  printf, 13,  ww, qq, ans, sig, ism[nn].ion[ww].state[qq].flgclm, $
                    ism[nn].ion[ww].state[qq].flginst, $
                    format='(i2,1x,i2,1x,f6.3,1x,f6.3,1x,i1,1x,i3)'
              endif
          endfor
      endfor
      close,13

      ;;cccccccccccccccccccccccccccccc
      ;; DH
      case flag of
          1: 
          else: begin
              ism[nn].flgD = ism[nn].elm[99].flgclm 
              if ism[nn].elm[99].flgclm NE 0 then begin
                  ism[nn].DH = ism[nn].elm[99].clm / 10^ism[nn].NHI * 1e5
                  ism[nn].sigDH = ism[nn].elm[99].sigclm / 10^ism[nn].NHI * 1e5
              endif
          end
      endcase

      ;; Fe
      dla_calcfeh, ism, nn

      ;; Zn
      Z= 30
      ism[nn].flgZn = ism[nn].elm[Z].flgclm 
      if ism[nn].flgZn NE 0 then begin
          dla_calcabnd, ism, nn, Z, 1, ans, sig
          ism[nn].ZnH = ans
          ism[nn].sigZnH = sig
      endif

      ;; CII* 
      if keyword_set(flg_CII) then begin
          if ism[nn].ion[6].state[6].flgclm GE 0 then begin
              ism[nn].flgCII = ism[nn].ion[6].state[6].flgclm + 1
              x_logclm, ism[nn].ion[6].state[6].clm, $
                ism[nn].ion[6].state[6].sigclm, ans, sig
              ism[nn].CII = ans
              ism[nn].sigCII = sig
          endif else ism[nn].flgCII = 0
      endif
      ;; Alpha
      dla_calcalpha, ism, nn
  endfor

  ;; Write
  ism_writestr, ism, flag=flag
  print, 'ism_updabd: All done'

  return
end
