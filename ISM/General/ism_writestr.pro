;+ 
; NAME:
; ism_writestr
;  V1.1
;
; PURPOSE:
;    Given a ISM structre write the files out
;
; CALLING SEQUENCE:
;   ism_writestr, ism
;
; INPUTS:
;  ism --  ISM structure
;
; RETURNS:
;
; OUTPUTS:
;  Series of DLA files
;
; OPTIONAL KEYWORDS:
;  FLAG -- 0=Default; 1=GRB
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   ism_writestr, ism
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   01-Oct-2004 Written by JXP
;-
;------------------------------------------------------------------------------
pro ism_writestr, ism, FLAG=flag

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'ism_writestr, ism, FLAG= [v1.1]'
    return
  endif 
  if not keyword_set( FLAG ) then flag = 0

;
  close, /all
  ;;
  nism = n_elements(ism)

  for i=0,nism-1 do begin
      openw, 1, ism[i].ismfil
;		openr, 1, strmid(fil,0,xlc(fil))
      ;; Name
      case flag of 
          1: printf, 1,  ism[i].target, '! GRB name', format='(a60,a10)'
          else: printf, 1,  ism[i].target, '! Target name', format='(a60,a13)'
      endcase
      printf, 1,  ism[i].targ_ra, '! RA (2000)', format='(a19,41x,a11)' 
      printf, 1,  ism[i].targ_dec, '! DEC (2000)', format='(a19,41x,a12)' 
      ;; Velocity
      case flag of 
          1: stop
          else: printf, 1,  ism[i].targ_vel, '! Target vel', $
            format='(f11.6,49x,a12)' 
      endcase
      printf, 1,  ism[i].flg_mag, '! flg_mag', format='(i2,58x,a9)' 
      printf, 1,  ism[i].targ_mag, '! Targ Mag', format='(f9.6,51x,a10)' 
      printf, 1,  ism[i].zabs, '! zabs', format='(f12.9,48x,a6)' 
      printf, 1,  ism[i].NHI, '! NHI', format='(f6.3,54x,a5)' 
      printf, 1,  ism[i].sigNHI[*], '! sig(NHI)', format='(2f6.3,48x,a10)' 
      printf, 1,  ism[i].abndfil, '! Abund file', format='(a60,a12)' 
      case flag of 
          1:
          else: begin
              ;; D
              if ism[i].flgD NE 0 then begin
                  printf, 1, ism[i].flgD, '! flg_D', format='(i3,57x,a7)' 
                  printf, 1, ism[i].DH, '! (D/H)x1e5', format='(f7.3,53x,a11)' 
                  printf, 1, ism[i].sigDH,'! sig(D/H)',format='(f7.3,53x,a10)' 
              endif else begin
                  printf, 1, '0', '! flg_D', format='(a1,59x,a7)'
                  printf, 1, '! [D/H]', format='(60x,a7)'
                  printf, 1, '! sig(D)', format='(60x,a8)' 
              endelse
          end
      endcase
      ;; Fe
      if ism[i].flgFe NE 0 then begin
          printf, 1, ism[i].flgFe, '! flg_Fe', format='(i3,57x,a8)' 
          printf, 1, ism[i].FeH, '! [Fe/H]', format='(f7.3,53x,a8)' 
          printf, 1, ism[i].sigFeH, '! sig(Fe)', format='(f7.3,53x,a9)' 
      endif else begin
          printf, 1, '0', '! flg_Fe', format='(a1,59x,a8)'
          printf, 1, '! [Fe/H]', format='(60x,a8)'
          printf, 1, '! sig(Fe)', format='(60x,a9)' 
      endelse
      ;; Zn
      if ism[i].flgZn NE 0 then begin
          printf, 1, ism[i].flgZn,'! flg_Zn', format='(i2,58x,a8)'
          printf, 1, ism[i].ZnH,'! [Zn/H]', format='(f7.3,53x,a8)' 
          printf, 1, ism[i].sigZnH,'! sig(Zn)', format='(f7.3,53x,a9)' 
      endif else begin
          printf, 1, '0', '! flg_Zn', format='(a1,59x,a8)'
          printf, 1, '! [Zn/H]', format='(60x,a8)'
          printf, 1, '! sig(Zn)', format='(60x,a9)' 
      endelse
      ;; Alpha
      if ism[i].flgAlpha NE 0 then begin
          printf, 1, ism[i].flgAlpha,'! flg_Alpha', format='(i2,58x,a11)'
          printf, 1, ism[i].Alpha,'! [Alpha/H]', format='(f7.3,53x,a11)' 
          printf, 1, ism[i].sigAlpha,'! sig(Alpha)', format='(f7.3,53x,a12)' 
      endif else begin
          printf, 1, '0', '! flg_Alpha', format='(a1,59x,a11)' 
          printf, 1, '! [Alpha/H]', format='(60x,a11)'
          printf, 1, '! sig(Alpha)', format='(60x,a12)' 
      endelse
      ;;    Low
      if ism[i].flglw NE 0 then begin
          printf, 1, ism[i].flglw, '! flg_low', format='(i2,58x,a9)'
          printf, 1, ism[i].lwfil, '! hi res file', format='(a60,a13)' 
          printf, 1, ism[i].lwwav, '! low_wav', format='(f9.4,51x,a9)' 
          printf, 1, ism[i].lwvmn, ism[i].lwvmx, $
            '! low_vmn,vmx (2f7)', format='(2f7.1,46x,a19)' 
          printf, 1, ism[i].lwfvel, '! fdelv', format='(f7.2,53x,a7)'
          printf, 1, ism[i].lwfmm, '! fmm', format='(f7.2,53x,a5)' 
          printf, 1, ism[i].lwfedg, '! fedg', format='(f7.2,53x,a6)' 
          printf, 1, ism[i].lwftpk, '! ftpk', format='(f7.2,53x,a6)' 
      endif else begin
          printf, 1, '0', '! flg_low', format='(a1,59x,a9)' 
          printf, 1, ism[i].lwfil, '! hi res file', format='(a60,a13)'
          printf, 1, '! low_wav', format='(60x,a9)' 
          printf, 1, '! low_vmn,vmx (2f7)', format='(60x,a19)'
          printf, 1, '! fdelv', format='(60x,a7)' 
          printf, 1, '! fmm', format='(60x,a5)' 
          printf, 1, '! fedg', format='(60x,a6)' 
          printf, 1, '! ftpk', format='(60x,a6)' 
      endelse
      ;;  CII*
      if ism[i].flgCII NE 0  then begin
          printf, 1, ism[i].flgCII,'! flg_CII', format='(i2,58x,a9)'   
          printf, 1, ism[i].CII,'! N(CII*)', format='(f6.3,54x,a9)' 
          printf, 1, ism[i].sigCII,'! sig(CII*)', format='(f5.3,55x,a11)'
      endif else begin
          printf, 1, '0', '! flg_CII', format='(a1,59x,a9)'
          printf, 1, '! N(CII*)', format='(60x,a9)'
          printf, 1, '! sig(CII*)', format='(60x,a11)' 
      endelse
      ;; CIV
      if ism[i].flgciv NE 0 then begin
          printf, 1, ism[i].flgciv, '! flg_civ', format='(i2,58x,a9)' 
          printf, 1, ism[i].civfil, '! civ hi res file', format='(a60,a17)' 
          printf, 1, ism[i].civwav, '! civ_wav', format='(f9.4,51x,a9)' 
          printf, 1, ism[i].civvmn, ism[i].civvmx, $
            '! civ_vmn,vmx (2f7)', format='(2f7.1,46x,a19)' 
          printf, 1, ism[i].civfvel, '! civ fdelv', format='(f7.2,53x,a11)' 
          printf, 1, ism[i].civfmm, '! civ fmm', format='(f7.2,53x,a9)' 
          printf, 1, ism[i].civfedg, '! civ fedg', format='(f7.2,53x,a10)' 
          printf, 1, ism[i].civftpk, '! civ ftpk', format='(f7.2,53x,a10)' 
      endif else begin
          printf, 1, '0', '! flg_civ', format='(a1,59x,a9)' 
          printf, 1, '! civ hi res file', format='(60x,a17)' 
          printf, 1, '! civ_wav', format='(60x,a9)' 
          printf, 1, '! civ_vmn,vmx (2f7)', format='(60x,a19)'
          printf, 1, '! civ fdelv', format='(60x,a11)'
          printf, 1, '! civ fmm', format='(60x,a9)'
          printf, 1, '! civ fedg', format='(60x,a10)' 
          printf, 1, '! civlow ftpk', format='(60x,a13)' 
      endelse
      case flag of
          1:
          else: begin
              ;; Cross-correlation
              if ism[i].flglw EQ 1 AND ism[i].flgciv EQ 1 then begin
                  printf, 1, ism[i].civlwfdv, '! civ fdv', format='(f7.2,53x,a9)'
                  printf, 1, ism[i].civlwfrto, '! civ frto', format='(f7.2,53x,a10)'
                  printf, 1, ism[i].civlwfnmm, '! civ fnmm', format='(f7.2,53x,a10)' 
                  printf, 1, ism[i].civlwftvm, '! civ ftvm', format='(f7.2,53x,a10)' 
              endif else begin
                  printf, 1, '! civlow fdv', format='(60x,a12)' 
                  printf, 1, '! civlow frto', format='(60x,a13)'
                  printf, 1, '! civlow fnmm', format='(60x,a13)' 
                  printf, 1, '! civlow ftvm', format='(60x,a13)' 
              endelse
          end
      endcase
      ;; E(B-V)
      printf, 1, ism[i].ebv, '! E(B-V)',format='(f5.3,55x,a8)'
      ;; Imaging, QSO survey
      case flag of 
          1:
          9: begin
              ;; Imaging stuff
              printf, 1, ism[i].ffilt, '! Filt 124816',format='(i2,58x,a13)' 
              printf, 1, ism[i].fslit, '! Slit 0n1s2y',format='(i2,58x,a13)'
              ;; QSO Survey
              printf, 1, ism[i].srvy, '! QSO Survey',format='(i2,58x,a12)'
              printf, 1, ism[i].srvy_mag,  '! Survey Mag',format='(f5.2,55x,a12)'
          end
          else:
      endcase
      ;; References
      printf, 1, ism[i].ref, '! References', format='(a60,a12)'
      ;; Metallicity
      if ism[i].flgmtl NE 0  then begin
          printf, 1, ism[i].flgmtl,'! flg_mtl', format='(i4,56x,a9)'
          printf, 1, ism[i].mtl,'! [M/H]', format='(f7.3,53x,a7)'
          printf, 1, ism[i].sigmtl,'! sig([M/H])', format='(f7.3,53x,a12)'
      endif else begin
          printf, 1, '0', '! flg_mtl', format='(a1,59x,a9)'
          printf, 1, '! [M/H]', format='(60x,a7)'
          printf, 1, '! sig([M/H])', format='(60x,a12)'
      endelse
      close, 1
  endfor
  close, /all

  print, 'ism_writestr: All done'

  return
end
