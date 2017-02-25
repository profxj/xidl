;+ 
; NAME:
; ism_corrzn
;  V1.2
;
; PURPOSE:
;    Corrects Zn columns based on MgI and CrII values
;    Primarily for ESI or other low resolution data
;
; CALLING SEQUENCE:
;   ism_corrzn, dla, qq, ZN_STR=, /Si
;
; INPUTS:
;   dla -- DLA structure
;   qq  -- Index of the structure
;
; RETURNS:
;
; OUTPUTS:
;  Series of DLA files
;
; OPTIONAL KEYWORDS:
;  /SI  --  Measure SiII 1808 EW too (MSDLA paper)
;
; OPTIONAL OUTPUTS:
;  ZN_STR=  -- Structure containing the info
;
; COMMENTS:
;
; EXAMPLES:
;   dla_allabd, dla
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   28-Mar-2006 Written by JXP
;-
;------------------------------------------------------------------------------
pro ism_corrzn, dla, qq, ZN_STR=zn_str, SI=si

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'ism_corrzn, dla, qq, ZN_STR=, /SI [v1.1]'
    return
  endif 


  ;; Create the structrue
  c = x_constants()
  ndla = n_elements(dla)

  zn_str = { $
          lin_ew: fltarr(2,2), $
          zn_ew: fltarr(2,2), $
          zn_colm: fltarr(3,2), $
          si_ew: fltarr(2), $  ; For MSDLA paper only
          cr_ew: fltarr(2), $
          mg_ew: fltarr(2,2) $ ; 2026, 2852
        }

  dumc = ''
  dumf = 0.0
  dumf2 = 0.0
  dumd = 0.d
  dumd2 = 0.d
  dumi = 0

  ;; Measure EWs
  openr, lun, strtrim(dla[qq].abndfil,2), /get_lun
  readf, lun, dumc
  readf, lun, dumi
  if dla[qq].ndfil MOD 2 GT 0 then readf, lun, dumc
  if dla[qq].ndfil MOD 4 GT 1 then readf, lun, dumc
  if dla[qq].ndfil MOD 8 GT 3 then readf, lun, dumc
  if dla[qq].ndfil MOD 16 GT 7 then readf, lun, dumc
  readf, lun, dumd
  readf, lun, dumc
  readf, lun, dumd, dumd2
  readf, lun, nhand
  for i=0L,nhand-1 do begin
      readf, lun, Znmb
      readf, lun, dumf, dumf2, dumi
  endfor
  
  flg_zn = 0
  while not eof(lun) do begin
      readf, lun, dumi
      ionflg = dumi
      readf, lun, dumd, dumf, dumf2, dumi
      
      ;; 2026
      if abs(dumd-2026.) LT 0.5 and flg_zn EQ 0 then begin
          flg_zn = 1
          
          ;; Read in data
          mdflg = dumi mod 64
          case mdflg of
              1: dat = x_readspec(dla[qq].Hfil, /struct, /autofsig)
              2: dat = x_readspec(dla[qq].Efil, /struct, /autofsig)
              4: dat = x_readspec(dla[qq].Ufil, /struct, /autofsig)
              8: dat = x_readspec(dla[qq].Xfil, /struct, /autofsig)
          endcase
          
          ;; Measure Rest EW
          x_pixminmax, dat.wv, dumd, dla[qq].zabs, dumf, dumf2, $
                       PIXMIN=pmin, PIXMAX=pmax
          ew = x_calcew(dat.wv, dat.fx, [pmin,pmax], $
                        dat.sig, sigew, /fpix)
          rew = ew / (1. + dla[qq].zabs)
          rsig = sigew / (1. + dla[qq].zabs)
          
          ;; Save
          zn_str.lin_ew[0,0] = rew * 1000. ; mA
          zn_str.lin_ew[0,1] = rsig * 1000. ; mA
      endif

      ;; 1808  (MSDLA paper only)
      if abs(dumd-1808.) LT 0.5 and keyword_set(SI) then begin
          
          ;; Read in data
          mdflg = dumi mod 64
          case mdflg of
              1: dat = x_readspec(dla[qq].Hfil, /struct, /autofsig)
              2: dat = x_readspec(dla[qq].Efil, /struct, /autofsig)
              4: dat = x_readspec(dla[qq].Ufil, /struct, /autofsig)
              8: dat = x_readspec(dla[qq].Xfil, /struct, /autofsig)
          endcase
          
          ;; Measure Rest EW
          x_pixminmax, dat.wv, dumd, dla[qq].zabs, dumf, dumf2, $
                       PIXMIN=pmin, PIXMAX=pmax
          ew = x_calcew(dat.wv, dat.fx, [pmin,pmax], $
                        dat.sig, sigew, /fpix)
          rew = ew / (1. + dla[qq].zabs)
          rsig = sigew / (1. + dla[qq].zabs)
          
          ;; Save
          zn_str.si_ew[0] = rew * 1000. ; mA
          zn_str.si_ew[1] = rsig * 1000. ; mA
      endif

      ;; 2062
      flg_2062 = 0
      if abs(dumd-2062.) LT 0.5 and flg_2062 EQ 0 then begin
          flg_2062 = 1
          
          ;; Read in data
          mdflg = dumi mod 64
          case mdflg of
              1: dat = x_readspec(dla[qq].Hfil, /struct, /autofsig)
              2: dat = x_readspec(dla[qq].Efil, /struct, /autofsig)
              4: dat = x_readspec(dla[qq].Ufil, /struct, /autofsig)
              8: dat = x_readspec(dla[qq].Xfil, /struct, /autofsig)
          endcase
          
          ;; Measure Rest EW
          x_pixminmax, dat.wv, dumd, dla[qq].zabs, dumf, dumf2, $
                       PIXMIN=pmin, PIXMAX=pmax
          ew = x_calcew(dat.wv, dat.fx, [pmin,pmax], $
                        dat.sig, sigew, /fpix)
          rew = ew / (1. + dla[qq].zabs)
          rsig = sigew / (1. + dla[qq].zabs)
          
          ;; Save
          zn_str.lin_ew[1,0] = rew * 1000. ; mA
          zn_str.lin_ew[1,1] = rsig * 1000. ; mA
      endif

      ;; 2852
      if abs(dumd-2853.) LT 0.5 then begin
          flg_zn = 1
          
          ;; Read in data
          mdflg = dumi mod 64
          case mdflg of
              1: dat = x_readspec(dla[qq].Hfil, /struct, /autofsig)
              2: dat = x_readspec(dla[qq].Efil, /struct, /autofsig)
              4: dat = x_readspec(dla[qq].Ufil, /struct, /autofsig)
              8: dat = x_readspec(dla[qq].Xfil, /struct, /autofsig)
          endcase
          
          ;; Measure Rest EW
          x_pixminmax, dat.wv, dumd, dla[qq].zabs, dumf, dumf2, $
                       PIXMIN=pmin, PIXMAX=pmax
          ew = x_calcew(dat.wv, dat.fx, [pmin,pmax], $
                        dat.sig, sigew, /fpix)
          rew = ew / (1. + dla[qq].zabs)
          rsig = sigew / (1. + dla[qq].zabs)
          
          ;; Save 2852
          zn_str.mg_ew[1,0] = rew * 1000. ; mA
          zn_str.mg_ew[1,1] = rsig * 1000. ; mA
      endif
      
  endwhile
  
  close, lun
  
  ;; MgI 2026 EW
  if dla[qq].ion[12].state[1].flgclm EQ 2 then mgoff = 0.1 else mgoff=0.
  nval = dla[qq].ion[12].state[1].clm+mgoff
  sigv = dla[qq].ion[12].state[1].sigclm
  ewv = ew_to_colm([2026.4768,2026.4768], [nval, sigv], $
                   /RVRS, /silent)
  zn_str.mg_ew[0,0] = ewv[0]
  case dla[qq].ion[12].state[1].flgclm of
      0: zn_str.mg_ew[0,1] = 0.
      1: zn_str.mg_ew[0,1] = ewv[1] 
      2: zn_str.mg_ew[0,1] = ewv[1]
      3: zn_str.mg_ew[0,1] = -9
      else: stop
  endcase
  
  ;; ZnII 2026 Colm
  if zn_str.lin_ew[0,0] GT 0. then begin
      ;; Value
      if zn_str.lin_ew[0,0] GT 3.*zn_str.lin_ew[0,1] then begin
          zn_str.zn_ew[0,0] = zn_str.lin_ew[0,0] $
            - (zn_str.mg_ew[0,0] > 0.)
          zn_str.zn_ew[0,1] = sqrt(zn_str.lin_ew[0,1]^2 $
                                       + zn_str.mg_ew[0,1]^2)
          colv = ew_to_colm([2026.136,2026.136], $
                            zn_str.zn_ew[0,*], /sile)
          zn_str.zn_colm[0,0] = colv[0]
          ;; Flg corresponding to 2026.
          mt = where(abs(dla[qq].ion[30].state[2,1:*].lambda-2026.2) LT 0.5, $
                     nmt)
          if nmt EQ 1 then flg2026 = dla[qq].ion[30].state[2,1+mt[0]].flgclm $
          else flg2026 = 0
          ;; Limit?
          if colv[0] GT 3.*colv[1] AND flg2026 NE 5 then begin
              zn_str.zn_colm[0,0] = colv[0]
              zn_str.zn_colm[0,1] = colv[1] 
          endif else begin
              if flg2026 EQ 5 then zn_str.zn_colm[0,0] = colv[0] + 2*colv[1] $
              else zn_str.zn_colm[0,0] = 3*colv[1] 
              zn_str.zn_colm[0,1] = -9
          end
      endif else begin ;; Limit
          zn_str.zn_ew[0,0] = 3*zn_str.lin_ew[0,1] 
          colv = ew_to_colm([2026.136], $
                            [zn_str.zn_ew[0,0]], /sile)
          zn_str.zn_colm[0,0] = colv[0]
          zn_str.zn_colm[0,1] = -9
      endelse
  endif
  
  
  ;; CrII 2062 EW
  if dla[qq].ion[24].state[2].flgclm EQ 1 then begin
      nval= dla[qq].ion[24].state[2].clm
      sigv= dla[qq].ion[24].state[2].sigclm
      ewv = ew_to_colm([2062.234,2062.234], [nval, sigv], $
                       /RVRS, /silent)
      zn_str.cr_ew[0] = ewv[0]
      zn_str.cr_ew[1] = ewv[1]
  endif
  
  ;; ZnII 2062 Colm
  if zn_str.lin_ew[1,0] GT 0. then begin
      ;; Value
      if zn_str.lin_ew[1,0] GT 3.*zn_str.lin_ew[1,1] then begin
          zn_str.zn_ew[1,0] = zn_str.lin_ew[1,0] $
            - (zn_str.cr_ew[0] > 0.)
          zn_str.zn_ew[1,1] = sqrt(zn_str.lin_ew[1,1]^2 $
                                       + zn_str.cr_ew[1]^2)
          colv = ew_to_colm([2062.664,2062.664], $
                            zn_str.zn_ew[1,*], /sile)
          ;; Flg corresponding to 2062.
          mt = where(abs(dla[qq].ion[30].state[2,1:*].lambda-2062.2) LT 0.5, $
                     nmt)
          if nmt EQ 1 then flg2062 = dla[qq].ion[30].state[2,1+mt[0]].flgclm $
          else flg2062 = 0
          ;; Limit?
          if colv[0] GT 3.*colv[1] or flg2062 EQ 5 then begin
              zn_str.zn_colm[1,0] = colv[0]
              zn_str.zn_colm[1,1] = colv[1] 
          endif else begin
              if flg2062 EQ 5 then zn_str.zn_colm[1,0] = colv[0] + 2*colv[1] $
              else zn_str.zn_colm[1,0] = 4*colv[1] 
              zn_str.zn_colm[1,1] = -9
          end
      endif else begin ;; Limit
          zn_str.zn_ew[1,0] = 3*zn_str.lin_ew[1,1] 
          colv = ew_to_colm([2062.664],$
                            [zn_str.zn_ew[1,0]], /sile)
          zn_str.zn_colm[1,0] = colv[0]
          zn_str.zn_colm[1,1] = -9
      endelse
  endif

  ;; Combine
  gd = where(zn_str.zn_colm[0:1,1] GT 0., ngd)
  case ngd of 
      0: begin                  ; Limits
          val = where(zn_str.zn_colm[0:1,0] GT 0., nval)
          if nval NE 0 then begin
              zn_str.zn_colm[2,0] = min(zn_str.zn_colm[val,0])
              zn_str.zn_colm[2,1] = -9
          endif
      end
      1: begin                  ; One value
          zn_str.zn_colm[2,0] = zn_str.zn_colm[gd,0]
          zn_str.zn_colm[2,1] = zn_str.zn_colm[gd,1]
      end
      else: begin               ; Weighted mean
          zn_str.zn_colm[2,0] = $
            total( zn_str.zn_colm[gd,0] / zn_str.zn_colm[gd,1]^2)/$
            total( 1. / zn_str.zn_colm[gd,1]^2)
          zn_str.zn_colm[2,1] = $
            sqrt(1. / total( 1. / zn_str.zn_colm[gd,1]^2))
          
      endelse
  endcase

  ;; Update Zn
  nonz = where(zn_str.zn_colm[0:1,1] NE 0., nnz)
  dla[qq].ion[30].indx[2] = nnz

  for jj=0L,nnz-1 do begin
      idx = nonz[jj]
      ;; Fill in
      case idx of
          0: dla[qq].ion[30].state[2,jj+1].lambda = 2026.136
          1: dla[qq].ion[30].state[2,jj+1].lambda = 2062.664
          else: stop
      endcase
      dla[qq].ion[30].state[2,jj+1].clm = zn_str.zn_colm[idx,0]
      dla[qq].ion[30].state[2,jj+1].sigclm = zn_str.zn_colm[idx,1]
      if zn_str.zn_colm[idx,1] GT 0. then $
        dla[qq].ion[30].state[2,jj+1].flgclm = 1 $
      else dla[qq].ion[30].state[2,jj+1].flgclm = 5
      if dla[qq].ion[30].state[2,jj+1].flginst EQ 0 then $
        dla[qq].ion[30].state[2,jj+1].flginst = dla[qq].ion[30].state[2,jj].flginst
  endfor
  ;; Final value
  if nnz NE 0 then begin
      ;; Ion
      dla[qq].ion[30].state[2].clm = zn_str.zn_colm[2,0]
      dla[qq].ion[30].state[2].sigclm = zn_str.zn_colm[2,1]
      if zn_str.zn_colm[2,1] GT 0. then $
        dla[qq].ion[30].state[2].flgclm = 1 $
      else dla[qq].ion[30].state[2].flgclm = 3
      ;; Element
      dla[qq].elm[30].clm = zn_str.zn_colm[2,0]
      dla[qq].elm[30].sigclm = zn_str.zn_colm[2,1]
      dla[qq].elm[30].flginst = dla[qq].ion[30].state[2,1].flginst
      if zn_str.zn_colm[2,1] GT 0. then $
        dla[qq].elm[30].flgclm = 1 $
      else dla[qq].elm[30].flgclm = 3
  endif

end
  
