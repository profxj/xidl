;+ 
; NAME:
; lls_writestr
;  V1.2
;
; PURPOSE:
;    Given a LLS structre write the files out
;
; CALLING SEQUENCE:
;   lls_writestr, lls
;
; INPUTS:
;  lls --  Array of LLS structures
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
;   lls_writestr, lls
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   01-Oct-2004 Written by JXP
;-
;------------------------------------------------------------------------------
pro lls_writestr, lls

; parse_llslst -- Reads in DLA data to a structure

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'lls_writestr, lls, [v1.1]'
    return
  endif 

;
  close, /all
  ;;
  nlls = n_elements(lls)

  for i=0,nlls-1 do begin
      openw, 1, lls[i].llsfil
;		openr, 1, strmid(fil,0,xlc(fil))
      padnm = x_padstr(lls[i].qso,60L,/trim)
      printf, 1,  padnm, '! QSO name', format='(a60,a10)'
      padra = x_padstr(lls[i].qso_ra,19L,/trim)
      printf, 1,  padra, '! RA (2000)', format='(a19,41x,a11)' 
      paddec = x_padstr(lls[i].qso_dec,19L,/trim)
      printf, 1,  paddec, '! DEC (2000)', format='(a19,41x,a12)' 
      printf, 1,  lls[i].qso_zem, '! QSO zem', format='(f9.6,51x,a9)' 
      printf, 1,  lls[i].flg_mag, '! flg_QSOmag', format='(i2,58x,a12)' 
      printf, 1,  lls[i].qso_mag, '! QSO Mag', format='(f9.6,51x,a9)' 
      printf, 1,  lls[i].srvy, '! QSO Survey', format='(i2,58x,a12)'
      printf, 1,  lls[i].srvy_mag, '! Survey Mag', format='(f9.5,51x,a9)' 
      padref = x_padstr(lls[i].ref,60L,/trim)
      printf, 1,  padref, '! References', format='(a60,a10)'
      printf, 1, lls[i].sdss_plate,lls[i].sdss_fibid,lls[i].sdss_mjd, $
        '! SDSS plt,fib,mjd', format='(i5,1x,i5,1x,i6,42x,a18)'
      printf, 1,  lls[i].zabs, '! zabs', format='(f9.5,51x,a6)' 
      printf, 1,  lls[i].nhi, '! NHI tot', format='(f9.4,51x,a9)' 
      printf, 1,  lls[i].signhi[0], lls[i].signhi[1], '! NHI sig', format='(f9.4,f9.4,42x,a9)' 
      printf, 1,  lls[i].nh, '! NH tot', format='(f9.4,50x,a9)' 
      printf, 1,  lls[i].nhsig, '! NH sig', format='(f9.4,f9.4,41x,a9)' 
      printf, 1, lls[i].vmn, lls[i].vmx, $
        '! vmn,vmx (2f7)', format='(2f7.1,42x,a19)' 
      printf, 1,  lls[i].fdelv, '! fdelv', format='(f9.2,49x,a9)' 
      printf, 1,  lls[i].fmm, '! fmm', format='(f9.2,47x,a9)' 
      printf, 1,  lls[i].fedg, '! fedg', format='(f9.2,48x,a9)' 
      printf, 1,  lls[i].ftpk, '! ftpk', format='(f9.2,48x,a9)' 
      printf, 1,  lls[i].flg_mh, '! flg_M/H', format='(i2,55x,a12)'
      printf, 1,  lls[i].mhave, '! [M/H] ave', format='(f9.2,51x,a11)' 
      ;printf, 1,  lls[i].mhsig, '! sig[M/H] ', format='(f9.2,51x,a11)' 
      printf, 1,  lls[i].mhsig, '! sig[M/H] ', format='(f9.2,1x,f9.2,41x,a11)' 
      printf, 1,  lls[i].flg_dh, '! flg_D/H', format='(i2,55x,a12)'
      printf, 1,  lls[i].dh, '!  [D/H]', format='(f9.6,50x,a9)' 
      printf, 1,  lls[i].nsys, '! N subsys', format='(i2,56x,a12)'
      padcldy = x_padstr(lls[i].cldyfil,60L,/trim)
      printf, 1,  padcldy, '! Cloudy Grid File', format='(a60,a18)' 
      for j=0,lls[i].nsys-1 do begin
          nm = lls[i].systems[j].name
          padsys = x_padstr(lls[i].systems[j].name,60L,/trim)
          printf, 1,  padsys, '! System '+nm, format='(a60,a10)'
          printf, 1,  lls[i].systems[j].zabs, '! '+nm+' zabs', format='(f9.6,50x,a9)' 
          printf, 1,  lls[i].systems[j].nhi, '! '+nm+' NHI', format='(f9.4,49x,a9)' 
          printf, 1,  lls[i].systems[j].nhisig, '! '+nm+' NHIsig', $
                  format='(f9.4,1x,f9.4,41x,a10)' 
          printf, 1,  lls[i].systems[j].nh, '! '+nm+' NH', format='(f9.4,48x,a9)' 
          printf, 1,  lls[i].systems[j].nhsig, '! '+nm+' NHsig', $
                  format='(f9.4,1x,f9.4,40x,a10)' 
          printf, 1,  lls[i].systems[j].logx, '! '+nm+' log x', format='(f9.4,47x,a13)' 
          printf, 1,  lls[i].systems[j].sig_logx, '! '+nm+' sigx', $
                  format='(f9.4,1x,f9.4,39x,a10)' 
          printf, 1,  lls[i].systems[j].b, '! '+nm+' b', format='(f9.4,47x,a9)' 
          printf, 1,  lls[i].systems[j].bsig, '! '+nm+' bsig', format='(f9.4,50x,a9)' 
          padabnd = x_padstr(lls[i].systems[j].abndfil,60L,/trim)
          printf, 1,  padabnd, '! '+nm+' Abund file', format='(a60,a14)' 
          printf, 1,  lls[i].systems[j].U, '! '+nm+' U', format='(f9.2,47x,a9)' 
          printf, 1,  lls[i].systems[j].Usig, '! '+nm+' Usig', $
                  format='(f9.4,1x,f9.4,39x,a10)' 
          printf, 1,  lls[i].systems[j].flg_low, '! '+nm+' flg_low', format='(i2,57x,a12)'
          printf, 1,  lls[i].systems[j].flg_alpha, '! '+nm+' flg_alpha', format='(i2,58x,a13)'
          printf, 1,  lls[i].systems[j].alphaH, '! '+nm+' [alpha/H]', format='(f9.6,51x,a13)' 
          printf, 1,  lls[i].systems[j].sig_alphaH, '! '+nm+' sig[a/H]', format='(f9.6,1x,f9.6,41x,a12)' 
          printf, 1,  lls[i].systems[j].flgfe,  '! '+nm+' flg_Fe', format='(i2,56x,a12)'
          printf, 1,  lls[i].systems[j].feh, '! '+nm+' [Fe/H]', format='(f9.6,51x,a10)' 
          printf, 1,  lls[i].systems[j].sig_feh, '! '+nm+' sig[Fe/H]', format='(f9.6,1x,f9.6,41x,a13)' 
          ;; VPFIT
          padvpfil = x_padstr(lls[i].systems[j].vpfil,60L,/trim)
          printf, 1,  padvpfil, '! '+nm+' VPFIT file', format='(a60,a14)' 
      endfor

      close, 1
  endfor
  close, /all

  print, 'lls_writestr: All done'

  return
end
