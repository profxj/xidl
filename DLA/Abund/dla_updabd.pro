;+ 
; NAME:
; dla_updabd
;  V1.1
;
; PURPOSE:
;    Given a list of DLA .dat files, calculate the abundances, modifty
;    the .dat files and write everything out.  This is the main driver
;    of AODM analysis for DLA.
;
; CALLING SEQUENCE:
;   dla_updabd, dla, /FLG_PLT, /FLG_CII, /FIXZN, /FINE
;
; INPUTS:
;  list -- List of DLA .dat files.
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /FIXZN -- Correct for MgI and CrII blending (MSDLA driven)
;  /FINE  -- Deal with fine-structure transitions (e.g. SiII*)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   dla_updabd, '~/Lists/all_mtl.lst'
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   01-Oct-2004 Written by JXP 
;-
;------------------------------------------------------------------------------
pro dla_updabd, list, FLG_PLT=flg_plt, FLG_CII=flg_cii, FIXZN=fixzn, $
                FINE=fine, EW=ew

; parse_dlalst -- Reads in DLA data to a structure

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'dla_updabd, list, /FLG_PLT, /FLG_CII, /FIXZN, /FINE [v1.1]'
    return
  endif 

  ;; Fix Zn check
  if keyword_set(FIXZN) then begin
      print, 'dla_updabd: Warning you are about to fiddle with Zn!'
      stop
  endif


  close, /all
  ;; Parse
  parse_dlalst, dla, list, /NORELM
  ndla = n_elements(dla)

  ;; Loop
  for nn=0L,ndla-1 do begin
      ;; Abund file?
      if strlen(strtrim(dla[nn].abndfil)) EQ 0 then begin
          print, 'dla_updabd:  No ABUND file for ', dla[nn].dlafil
          print, 'dla_updabd:  Skipping...'
          continue
      endif
      print, 'dla_updabd:  Processing ', dla[nn].dlafil

      ;; Abundances
      dla_allabnd, dla, nn, FLG_PLT=flg_plt

      ;; Fix Zn?
      if keyword_set(FIXZN) then $
        ism_corrzn, dla, nn

      ;; Fine structure?
      if keyword_set(FINE) then begin
          tmp = dla[nn]
          x_fineabnd, tmp, /NOLOG
          dla[nn] = tmp
      endif
        
      ;; Ion List
      len = strlen(dla[nn].tab_fil)
      ION_fil = strmid(dla[nn].tab_fil,0,len-3)+'ion'
      close, 13
      openw, 13, ION_fil
      for ww=1,99 do begin
          for qq=1,6 do begin
              tidx = dla[nn].ion[ww].indx[qq]
              for k=1,tidx do begin
                  if dla[nn].ion[ww].state[qq,k].flgclm NE -1 then begin
                      if ww NE 1  then begin
                          x_logclm, dla[nn].ion[ww].state[qq,k].clm, $
                            dla[nn].ion[ww].state[qq,k].sigclm, ans, sig
                      endif else begin
                          ans =  dla[nn].ion[ww].state[qq,k].clm
                          sig =  dla[nn].ion[ww].state[qq,k].sigclm
                      endelse
                  endif
                  if sig GT 9.9 then sig = 9.9
                  printf, 13,  $
                    dla[nn].ion[ww].state[qq,k].lambda, ans,sig, $
                    dla[nn].ion[ww].state[qq,k].flgclm, $
                    dla[nn].ion[ww].state[qq,k].flginst, $
                    format='(f10.4,1x,2f8.4,1x,i2,1x,i3)'
              endfor
          endfor
      endfor
      close, 13

      ;;cccccccccccccccccccccccccccccc
      ;; Output Element List
      ;;
      XH_fil = strmid(dla[nn].tab_fil,0,len-3)+'XH'
      openw, 13, XH_fil
      for ww=1,99 do begin
          if dla[nn].elm[ww].flgclm NE 0 then begin
              x_logclm, dla[nn].elm[ww].clm, dla[nn].elm[ww].sigclm, ans, sig
              if SIG GT 9.9 then sig = 9.9
              print, ww, ans, sig, dla[nn].elm[ww].flgclm, dla[nn].elm[ww].flginst, $
                format='(i2,1x,f6.3,1x,f6.3,1x,i1,1x,i3)'
              printf, 13, ww, ans, sig, dla[nn].elm[ww].flgclm, $
                dla[nn].elm[ww].flginst, $
                format='(i2,1x,f6.3,1x,f6.3,1x,i1,1x,i3)'
          endif 
      endfor
      close,13
      
      ;;cccccccccccccccccccccccccccccc
      ;; Output Total List
      ;;

      XH_fil = strmid(dla[nn].tab_fil,0,len-3)+'all'
      openw, 13, XH_fil
      dla[nn].ion[1].state[1].flgclm = 1 
      dla[nn].ion[1].state[1].clm =  10.^dla[nn].NHI
      dla[nn].ion[1].state[1].sigclm =  10.^(dla[nn].NHI+dla[nn].sigNHI[0]) - 10.^dla[nn].NHI
      dla[nn].ion[1].state[1].flginst = 99
      for ww=1,99 do begin
          for qq=1,6 do begin
              if dla[nn].ion[ww].state[qq].flgclm NE 0 then begin
                  x_logclm, dla[nn].ion[ww].state[qq].clm, $
                    dla[nn].ion[ww].state[qq].sigclm, ans, sig
                  if sig GT 9.9 then sig = 9.9
                  printf, 13,  ww, qq, ans, sig, dla[nn].ion[ww].state[qq].flgclm, $
                    dla[nn].ion[ww].state[qq].flginst, $
                    format='(i2,1x,i2,1x,f6.3,1x,f6.3,1x,i1,1x,i3)'
              endif
          endfor
      endfor
      close,13

      ;;cccccccccccccccccccccccccccccc
      ;; Fe
      dla_calcfeh, dla, nn, /nohsig, FINE=fine

      ;; Zn
      Z= 30
      dla[nn].flgZn = dla[nn].elm[Z].flgclm 
      if dla[nn].flgZn NE 0 then begin
          dla_calcabnd, dla, nn, Z, 1, ans, sig, /nosigy
          dla[nn].ZnH = ans
          dla[nn].sigZnH = sig
      endif
 
     ;; CII*  ;RAJ changed this to write out cii* properly
      if keyword_set(flg_CII) then begin
          if dla[nn].ion[6].state[6,1].flgclm GE 0 then begin
              dla[nn].flgCII = dla[nn].ion[6].state[6,1].flgclm + 1
              x_logclm, dla[nn].ion[6].state[6,1].clm, $
                dla[nn].ion[6].state[6,1].sigclm, ans, sig
              dla[nn].CII = ans
              dla[nn].sigCII = sig
          endif else dla[nn].flgCII = 0
      endif 

    ;; Alpha
      dla_calcalpha, dla, nn, /NOHSIG
      ;; Metallicity
      dla_calcmtl, dla, nn

      ;; VPFIT
      fil = file_search(strtrim(dla[nn].vpfit_fil,2), count=nfil)
      if nfil EQ 1 then dla_abund_vpfit, dla, nn
  endfor


  ;; Write
  dla_writestr, dla
  print, 'dla_updabd: All done'

  ;; EW?
  if keyword_set(EW) then dla_updew, list

  return
end
