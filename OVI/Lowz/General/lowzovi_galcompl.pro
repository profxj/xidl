;+ 
; NAME:
; lowzovi_galcompl   
;   Version 1.0
;
; PURPOSE:
;    Launches a cw_field and grabs input from the user
;
; CALLING SEQUENCE:
;   
;   string = lowzovi_galcompl(title)
;
; INPUTS:
;   title - Title
;
; RETURNS:
;   string - String
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
;   string = lowzovi_galcompl( 'Enter filename: ')
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   30-Sep-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro lowzovi_galcompl, list

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'lowzovi_galcompl, list [v1.0]'
    return
  endif 

;  Optional Keywords
  if not keyword_set( XOFFSET ) then xoffset = 200

; Parse List
  lowzovi_prsdat, struct, list, listnms=listnms
  nfield = n_elements(listnms)

; Loop on good files

  gd_gal = where(struct.flg_gal NE 0, ngd)

  for qq=0L,ngd-1 do begin
      ;; Open the file
      gal_str = xmrdfits(struct[gd_gal[qq]].gal_fil, 1, /silent)

   ; Count the galaxies
      good_z = where(gal_str.flg_anly MOD 8 GT 3 AND $
                     gal_str.obj_id EQ 'a', nz)

      ;; z > 0.005
      all_z = where(gal_str[good_z].z GT 0.005, nall)
      struct[gd_gal[qq]].N_gal[0] = nall
      ;; 0.11 > z > 0.005
      fuse_z = where(gal_str[good_z].z GT 0.005 AND $
                     gal_str[good_z].z LT 0.11, nfuse)
      struct[gd_gal[qq]].N_gal[1] = nfuse
      ;; z_qso > z > 0.005
      lqso_z = where(gal_str[good_z].z GT 0.005 AND $
                    gal_str[good_z].z LT (struct[gd_gal[qq]].qso_zem+0.01), $
                    nlqso)
      struct[gd_gal[qq]].N_gal[2] = nlqso
      ;; z~z_qso
      qso_z = where(abs(gal_str[good_z].z $
                        - struct[gd_gal[qq]].qso_zem) LT 0.01, nqso)
      struct[gd_gal[qq]].N_gal[3] = nqso
      print, 'Ngal = ', struct[gd_gal[qq]].N_gal

   ; Completeness
      all_surv = where(gal_str.flg_survey MOD 2 EQ 1 AND $
                       gal_str.obj_id EQ 'a', nsurv)
      impact = sqrt(gal_str[all_surv].ra^2 + $
             gal_str[all_surv].dec^2)  ;; Arcseconds
      cut195 = where(gal_str[all_surv].mag[1] LT 19.5, n195)
      ;; R < 19.5, p < 5'
      all_imp = where(impact[cut195] LT 300., nall)
      good_z = where(gal_str[all_surv[cut195[all_imp]]].flg_anly MOD 8 GT 3,ngd)
      ngd = ngd > 0
      struct[gd_gal[qq]].complete[0,0] = ngd
      struct[gd_gal[qq]].complete[0,1] = round(100*float(ngd)/float(nall))
      ;; R < 19.5, p < 10'
      all_imp = where(impact[cut195] LT 600., nall)
      good_z = where(gal_str[all_surv[cut195[all_imp]]].flg_anly MOD 8 GT 3,ngd)
      ngd = ngd > 0
      struct[gd_gal[qq]].complete[1,0] = ngd
      struct[gd_gal[qq]].complete[1,1] = round(100*float(ngd)/float(nall))
      if struct[gd_gal[qq]].R_limit GT 19.5 then begin
          cut200 = where(gal_str[all_surv].mag[1] LT 20.0, n200)
          ;; R < 20.0, p < 5'
          all_imp = where(impact[cut200] LT 300., nall)
          good_z = where(gal_str[all_surv[cut200[all_imp]]].flg_anly $
                         MOD 8 GT 3,ngd)
          ngd = ngd > 0
          struct[gd_gal[qq]].complete[2,0] = ngd
          struct[gd_gal[qq]].complete[2,1] = round(100*float(ngd)/float(nall))
          ;; R < 19.5, p < 10'
          all_imp = where(impact[cut200] LT 600., nall)
          good_z = where(gal_str[all_surv[cut200[all_imp]]].flg_anly $
                         MOD 8 GT 3,ngd)
          ngd = ngd > 0
          struct[gd_gal[qq]].complete[3,0] = ngd
          struct[gd_gal[qq]].complete[3,1] = round(100*float(ngd)/float(nall))
      endif
      print, 'Completeness: '
      print, struct[gd_gal[qq]].complete

  endfor
  ;; Write-out
  lowzovi_wrdat, struct, LIST=list

  return
end
