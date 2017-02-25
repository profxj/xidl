;+ 
; NAME:
; dla_writestr
;  V1.2
;
; PURPOSE:
;    Given a DLA structre write the files out with the .dat format
;
; CALLING SEQUENCE:
;   dla_writestr, dla
;
; INPUTS:
;  dla -- DLA structure
;
; RETURNS:
;
; OUTPUTS:
;  Series of DLA files
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   dla_writestr, dla
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   01-Oct-2004 Written by JXP
;   2006   Added metallicity and SDSS tags
;-
;------------------------------------------------------------------------------
pro dla_writestr, dla

; parse_dlalst -- Reads in DLA data to a structure

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'dla_writestr, dla, [v1.2]'
    return
  endif 

;
  close, /all
  ;;
  ndla = n_elements(dla)

  for i=0,ndla-1 do begin
      openw, 1, dla[i].dlafil
;		openr, 1, strmid(fil,0,xlc(fil))
      padnm = x_padstr(dla[i].qso,60L,/trim)
      printf, 1,  padnm, '! QSO name', format='(a60,a10)'
      padra = x_padstr(dla[i].qso_ra,19L,/trim)
      printf, 1,  padra, '! RA (2000)', format='(a19,41x,a11)' 
      paddec = x_padstr(dla[i].qso_dec,19L,/trim)
      printf, 1,  paddec, '! DEC (2000)', format='(a19,41x,a12)' 
      printf, 1,  dla[i].qso_zem, '! QSO zem', format='(f9.6,51x,a9)' 
      printf, 1,  dla[i].flg_QSOmag, '! flg_QSOmag', format='(i2,58x,a12)' 
      printf, 1,  dla[i].qso_mag, '! QSO Mag', format='(f9.6,51x,a9)' 
      printf, 1,  dla[i].zabs, '! zabs', format='(f12.9,48x,a6)' 
      printf, 1,  dla[i].NHI, '! NHI', format='(f6.3,54x,a5)' 
      printf, 1,  dla[i].sigNHI[*], '! sig(NHI)', format='(2f6.3,48x,a10)' 
      printf, 1,  dla[i].abndfil, '! Abund file', format='(a60,a12)' 
      ;; Fe
      if dla[i].flgFe NE 0 then begin
          printf, 1, dla[i].flgFe, '! flg_Fe', format='(i3,57x,a8)' 
          printf, 1, dla[i].FeH, '! [Fe/H]', format='(f7.3,53x,a8)' 
          printf, 1, dla[i].sigFeH, '! sig(Fe)', format='(f7.3,53x,a9)' 
      endif else begin
          printf, 1, '0', '! flg_Fe', format='(a1,59x,a8)'
          printf, 1, '! [Fe/H]', format='(60x,a8)'
          printf, 1, '! sig(Fe)', format='(60x,a9)' 
      endelse
      ;; Zn
      if dla[i].flgZn NE 0 then begin
          printf, 1, dla[i].flgZn,'! flg_Zn', format='(i2,58x,a8)'
          printf, 1, dla[i].ZnH,'! [Zn/H]', format='(f7.3,53x,a8)' 
          printf, 1, dla[i].sigZnH,'! sig(Zn)', format='(f7.3,53x,a9)' 
      endif else begin
          printf, 1, '0', '! flg_Zn', format='(a1,59x,a8)'
          printf, 1, '! [Zn/H]', format='(60x,a8)'
          printf, 1, '! sig(Zn)', format='(60x,a9)' 
      endelse
      ;; Alpha
      if dla[i].flgAlpha NE 0 then begin
          printf, 1, dla[i].flgAlpha,'! flg_Alpha', format='(i2,58x,a11)'
          printf, 1, dla[i].Alpha,'! [Alpha/H]', format='(f7.3,53x,a11)' 
          printf, 1, dla[i].sigAlpha,'! sig(Alpha)', format='(f7.3,53x,a12)' 
      endif else begin
          printf, 1, '0', '! flg_Alpha', format='(a1,59x,a11)' 
          printf, 1, '! [Alpha/H]', format='(60x,a11)'
          printf, 1, '! sig(Alpha)', format='(60x,a12)' 
      endelse
      ;;    Low
      if dla[i].flglw NE 0 then begin
          printf, 1, dla[i].flglw, '! flg_low', format='(i2,58x,a9)'
          printf, 1, dla[i].lwfil, '! hi res file', format='(a60,a13)' 
          printf, 1, dla[i].lwwav, '! low_wav', format='(f9.4,51x,a9)' 
          printf, 1, dla[i].lwvmn, dla[i].lwvmx, $
            '! low_vmn,vmx (2f7)', format='(2f7.1,46x,a19)' 
          printf, 1, dla[i].lwfvel, '! fdelv', format='(f7.2,53x,a7)'
          printf, 1, dla[i].lwfmm, '! fmm', format='(f7.2,53x,a5)' 
          printf, 1, dla[i].lwfedg, '! fedg', format='(f7.2,53x,a6)' 
          printf, 1, dla[i].lwftpk, '! ftpk', format='(f7.2,53x,a6)' 
      endif else begin
          printf, 1, '0', '! flg_low', format='(a1,59x,a9)' 
          printf, 1, dla[i].lwfil, '! hi res file', format='(a60,a13)'
          printf, 1, '! low_wav', format='(60x,a9)' 
          printf, 1, '! low_vmn,vmx (2f7)', format='(60x,a19)'
          printf, 1, '! fdelv', format='(60x,a7)' 
          printf, 1, '! fmm', format='(60x,a5)' 
          printf, 1, '! fedg', format='(60x,a6)' 
          printf, 1, '! ftpk', format='(60x,a6)' 
      endelse
      ;;  CII*
      if dla[i].flgCII NE 0  then begin
          printf, 1, dla[i].flgCII,'! flg_CII', format='(i2,58x,a9)'   
          printf, 1, dla[i].CII,'! N(CII*)', format='(f6.3,54x,a9)' 
          printf, 1, dla[i].sigCII,'! sig(CII*)', format='(f5.3,55x,a11)'
      endif else begin
          printf, 1, '0', '! flg_CII', format='(a1,59x,a9)'
          printf, 1, '! N(CII*)', format='(60x,a9)'
          printf, 1, '! sig(CII*)', format='(60x,a11)' 
      endelse
      ;; CIV
      if dla[i].flgciv NE 0 then begin
          printf, 1, dla[i].flgciv, '! flg_civ', format='(i2,58x,a9)' 
          printf, 1, dla[i].civfil, '! civ hi res file', format='(a60,a17)' 
          printf, 1, dla[i].civwav, '! civ_wav', format='(f9.4,51x,a9)' 
          printf, 1, dla[i].civvmn, dla[i].civvmx, $
            '! civ_vmn,vmx (2f7)', format='(2f7.1,46x,a19)' 
          printf, 1, dla[i].civfvel, '! civ fdelv', format='(f7.2,53x,a11)' 
          printf, 1, dla[i].civfmm, '! civ fmm', format='(f7.2,53x,a9)' 
          printf, 1, dla[i].civfedg, '! civ fedg', format='(f7.2,53x,a10)' 
          printf, 1, dla[i].civftpk, '! civ ftpk', format='(f7.2,53x,a10)' 
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
      ;; Cross-correlation
      if dla[i].flglw EQ 1 AND dla[i].flgciv EQ 1 then begin
          printf, 1, dla[i].civlwfdv, '! civ fdv', format='(f7.2,53x,a9)'
          printf, 1, dla[i].civlwfrto, '! civ frto', format='(f7.2,53x,a10)'
          printf, 1, dla[i].civlwfnmm, '! civ fnmm', format='(f7.2,53x,a10)' 
          printf, 1, dla[i].civlwftvm, '! civ ftvm', format='(f7.2,53x,a10)' 
      endif else begin
          printf, 1, '! civlow fdv', format='(60x,a12)' 
          printf, 1, '! civlow frto', format='(60x,a13)'
          printf, 1, '! civlow fnmm', format='(60x,a13)' 
          printf, 1, '! civlow ftvm', format='(60x,a13)' 
      endelse
      ;; Imaging stuff
      printf, 1, dla[i].qso_ebv, '! E(B-V)',format='(f5.3,55x,a8)'
      printf, 1, dla[i].ffilt, '! Filt 124816',format='(i2,58x,a13)' 
      printf, 1, dla[i].fslit, '! Slit 0n1s2y',format='(i2,58x,a13)'
      ;; QSO Survey
      printf, 1, dla[i].srvy, '! QSO Survey',format='(i2,58x,a12)'
      printf, 1, dla[i].srvy_mag,  '! Survey Mag',format='(f5.2,55x,a12)'
      ;; References
      printf, 1, x_padstr(dla[i].ref,60L,/trim), $
        '! References', format='(a60,a12)'
      ;; Metallicity
      if dla[i].flgmtl NE 0  then begin
          printf, 1, dla[i].flgmtl,'! flg_mtl', format='(i4,56x,a9)'
          printf, 1, dla[i].mtl,'! [M/H]', format='(f7.3,53x,a7)'
          printf, 1, dla[i].sigmtl,'! sig([M/H])', format='(f7.3,53x,a12)'
      endif else begin
          printf, 1, '0', '! flg_mtl', format='(a1,59x,a9)'
          printf, 1, '! [M/H]', format='(60x,a7)'
          printf, 1, '! sig([M/H])', format='(60x,a12)'
      endelse
      ;; SDSS plate, fiberid
      if (dla[i].srvy / 10 EQ 6) OR dla[i].sdss_plate GT 0 then begin
          printf, 1, dla[i].sdss_plate,dla[i].sdss_fibid,dla[i].sdss_mjd, $
            '! SDSS plt,fib,mjd', format='(i5,1x,i5,1x,i6,42x,a18)'
      endif else begin
          printf, 1, '    0,    0,     0', '! SDSS plt,fib,mjd', $
            format='(a18,42x,a18)'
      endelse
      ;; VPFIT fil
      printf, 1, dla[i].vpfit_fil, '! VPFIT FILE', format='(a60,a12)' 
      ;;  CI
      if dla[i].flg_CI NE 0  then begin
          printf, 1, dla[i].flg_CI,'! flg_CI', format='(i2,58x,a8)'   
          printf, 1, dla[i].CI,'! N(CI)', format='(f6.3,54x,a7)' 
          printf, 1, dla[i].sig_CI,'! sig(CI)', format='(f5.3,55x,a9)'
      endif else begin
          printf, 1, '0', '! flg_CI', format='(a1,59x,a8)'
          printf, 1, '! N(CI)', format='(60x,a7)'
          printf, 1, '! sig(CI)', format='(60x,a9)' 
      endelse
      ;;  H2
      if dla[i].flg_H2 NE 0  then begin
          printf, 1, dla[i].flg_H2,'! flg_H2', format='(i2,58x,a8)'   
          printf, 1, dla[i].H2,'! N(H2)', format='(f6.3,54x,a7)' 
          printf, 1, dla[i].sig_H2,'! sig(H2)', format='(f5.3,55x,a9)'
      endif else begin
          printf, 1, '0', '! flg_H2', format='(a1,59x,a8)'
          printf, 1, '! N(H2)', format='(60x,a7)'
          printf, 1, '! sig(H2)', format='(60x,a9)' 
      endelse
      close, 1
  endfor
  close, /all

  print, 'dla_writestr: All done'

  return
end
