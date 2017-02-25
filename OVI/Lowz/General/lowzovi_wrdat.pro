;+ 
; NAME:
; lowzovi_wrdat
;  V1.1
;  
; PURPOSE:
;  Writes to an ASCII file a summary of the galaxy data for a
;  list of fields.
;
; CALLING SEQUENCE:
;   lowzovi_wrdat, stucture, LIST=
;
; INPUTS:
;
; RETURNS:
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
;   lowzovi_wrdat, struct, '/u/xavier/DLA/Lists/tot_dla.lst'
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Oct-2003 Added metallicity sturcture
;-
;------------------------------------------------------------------------------
pro lowzovi_wrdat, supstrc, list

; lowzovi_wrdat -- Writes OVI data to files from a structure

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'lowzovi_wrdat, struct, list [v1.1]' 
    return
  endif 

;

  fil=''
  dumc = ''
  dumr = 0.0
  dumd = double(0.0)
  dumr2 = 0.0
  dumi = 0

; Parse the Base File
  readcol, list, listnms, format='A'
  nfield = n_elements(listnms)
;
  close, 1
  for i=0,nfield-1 do begin
      fil = listnms[i]
      openw, 1, fil
      printf, 1,  format='(a15,45x,a10)', supstrc[i].qso, '! QSO name'
      printf, 1,  format='(a15,45x,a11)', supstrc[i].qso_ra, '! RA (2000)'
      printf, 1,  format='(a15,45x,a12)', supstrc[i].qso_dec, '! DEC (2000)'
      printf, 1,  format='(f9.6,51x,a9)', supstrc[i].qso_zem, '! QSO zem'
      printf, 1,  format='(f9.5,51x,a7)', supstrc[i].qso_vmag, '! QSO V'
      printf, 1,  format='(f10.5,50x,a15)', supstrc[i].qso_uv, '! QSO UV (e-15)'
      ;; Galaxy
      printf, 1,  format='(i3,57x,a12)', supstrc[i].flg_gal, '! Galaxy flg'
      printf, 1,  format='(a60,a13)', supstrc[i].gal_fil, '! Galaxy file'
      ;; R limit + galaxies
      printf, 1,  format='(f9.3,51x,a9)', supstrc[i].R_limit, '! R Limit'
      printf, 1,  format='(i3,57x,a13)', supstrc[i].N_gal[0], '! N (z>0.005)'
      printf, 1,  format='(i3,57x,a18)', supstrc[i].N_gal[1], '! N (0.11>z>0.005)'
      printf, 1,  format='(i3,57x,a17)', supstrc[i].N_gal[2], '! N (qso>z>0.005)'
      printf, 1,  format='(i3,57x,a11)', supstrc[i].N_gal[3], '! N (z~qso)'
      ;; Completeness
      printf, 1,  format='(i3,57x,a18)', supstrc[i].complete[0,0], $
        '! N (R<19.5; p<5'')'
      printf, 1,  format='(i3,57x,a18)', supstrc[i].complete[0,1], $
        '! % (R<19.5; p<5'')'
      printf, 1,  format='(i3,57x,a19)', supstrc[i].complete[1,0], $
        '! N (R<19.5; p<10'')'
      printf, 1,  format='(i3,57x,a19)', supstrc[i].complete[1,1], $
        '! % (R<19.5; p<10'')'
      printf, 1,  format='(i3,57x,a18)', supstrc[i].complete[2,0], $
        '! N (R<20.0; p<5'')'
      printf, 1,  format='(i3,57x,a18)', supstrc[i].complete[2,1], $
        '! % (R<20.0; p<5'')'
      printf, 1,  format='(i3,57x,a19)', supstrc[i].complete[3,0], $
        '! N (R<20.0; p<10'')'
      printf, 1,  format='(i3,57x,a19)', supstrc[i].complete[3,1], $
        '! % (R<20.0; p<10'')'
      ;; FUSE
      printf, 1,  format='(i3,57x,a10)', supstrc[i].flg_fuse, '! FUSE flg'
      printf, 1,  format='(f11.2,49x,a10)', supstrc[i].fuse_exp, '! FUSE exp'
      printf, 1,  format='(f6.2,54x,a10)', supstrc[i].fuse_snr, '! FUSE S/N'

      ;; STIS
      printf, 1,  format='(i3,57x,a10)', supstrc[i].flg_stis, '! STIS flg'
      printf, 1,  format='(a60,a13)', supstrc[i].stis_comm[0], '! STIS Line 1'
      printf, 1,  format='(a60,a13)', supstrc[i].stis_comm[1], '! STIS Line 2'
      printf, 1,  format='(a60,a13)', supstrc[i].stis_comm[2], '! STIS Line 3'
      printf, 1,  format='(a60,a13)', supstrc[i].stis_comm[3], '! STIS Line 4'

      ;; GHRS
      printf, 1,  format='(i3,57x,a10)', supstrc[i].flg_ghrs, '! GHRS flg'
      printf, 1,  format='(a60,a13)', supstrc[i].ghrs_comm[0], '! GHRS Line 1'
      printf, 1,  format='(a60,a13)', supstrc[i].ghrs_comm[1], '! GHRS Line 2'
      printf, 1,  format='(a60,a13)', supstrc[i].ghrs_comm[2], '! GHRS Line 3'
      printf, 1,  format='(a60,a13)', supstrc[i].ghrs_comm[3], '! GHRS Line 4'

      ;; CLOSE
      close, 1


  endfor

  return
end
