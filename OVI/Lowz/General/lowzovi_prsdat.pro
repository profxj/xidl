;+ 
; NAME:
; lowzovi_prsdat
;  V1.1
;
; PURPOSE:
;    
; CALLING SEQUENCE:
;   lowzovi_prsdat, stucture, list, LISTNMS=
;
; INPUTS:
;  [list]  -- List of OVI files [default: tot_ovi.lst]
;
; RETURNS:
;   structure      - IDL structure
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;  LISTNMS=  -- Array of list names
;
; COMMENTS:
;
; EXAMPLES:
;   lowzovi_prsdat, struct
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Oct-2003 Added metallicity sturcture
;-
;------------------------------------------------------------------------------
pro lowzovi_prsdat, supstrc, list, listnms=listnms

; lowzovi_prsdat -- Reads in Galaxy data to a structure for the
;   OVI project

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'lowzovi_prsdat, struct, [filename], [v1.1]' 
    return
  endif 

;


  if not keyword_set( list ) then begin
      print, 'Using tot_ovi.lst as the list'
      list = '/u/xavier/OVI/Lists/tot_ovi.lst'
  endif

  fil=''
  dumc = ''
  dumr = 0.0
  dumd = double(0.0)
  dumr2 = 0.0
  dumi = 0

; Parse the Base File
  tmp = { lowzovidat }
  readcol, list, listnms, format='A'
  nfield = n_elements(listnms)
  supstrc = replicate(tmp,nfield)
;
  close, 1
  for i=0,nfield-1 do begin
      fil = listnms[i]
      openr, 1, fil
;		openr, 1, strmid(fil,0,xlc(fil))
      readf, 1,  format='(a15)', dumc
      supstrc[i].qso = dumc
      readf, 1,  format='(a15)', dumc
      supstrc[i].qso_ra = dumc
      readf, 1,  format='(a15)', dumc
      supstrc[i].qso_dec = dumc
      readf, 1,  format='(f9.6)', dumr
      supstrc[i].qso_zem = dumr
      readf, 1,  format='(f9.5)', dumr
      supstrc[i].qso_vmag = dumr
      readf, 1,  format='(f10.5)', dumr
      supstrc[i].qso_uv = dumr
      ; Galaxy file
      readf, 1,  format='(i3)', dumi
      supstrc[i].flg_gal = dumi   
      readf, 1,  format='(a60)', dumc
      supstrc[i].gal_fil = strtrim(dumc,2)
      ; R limit
      readf, 1,  format='(f9.3)', dumr
      supstrc[i].R_limit = dumr
      ; N galaxies
      readf, 1,  format='(i3)', dumi  ; z>0.005
      supstrc[i].N_gal[0] = dumi   
      readf, 1,  format='(i3)', dumi  ; 0.11>z>0.005
      supstrc[i].N_gal[1] = dumi   
      readf, 1,  format='(i3)', dumi  ; z_qso>z>0.005
      supstrc[i].N_gal[2] = dumi   
      readf, 1,  format='(i3)', dumi  ; z~z_qso
      supstrc[i].N_gal[3] = dumi   
      ; Completeness
      readf, 1,  format='(i3)', dumi  ; R<19.5; p<5'
      supstrc[i].complete[0,0] = dumi
      readf, 1,  format='(i3)', dumi  
      supstrc[i].complete[0,1] = dumi
      readf, 1,  format='(i3)', dumi  ; R<19.5; p<10'
      supstrc[i].complete[1,0] = dumi
      readf, 1,  format='(i3)', dumi  
      supstrc[i].complete[1,1] = dumi
      readf, 1,  format='(i3)', dumi  ; R<20.0; p<5'
      supstrc[i].complete[2,0] = dumi
      readf, 1,  format='(i3)', dumi  
      supstrc[i].complete[2,1] = dumi
      readf, 1,  format='(i3)', dumi  ; R<20.0; p<10'
      supstrc[i].complete[3,0] = dumi
      readf, 1,  format='(i3)', dumi  
      supstrc[i].complete[3,1] = dumi
      ; FUSE
      readf, 1,  format='(i3)', dumi
      supstrc[i].flg_fuse = dumi   
      readf, 1,  format='(f11.2)', dumr
      supstrc[i].fuse_exp = dumr
      readf, 1,  format='(f6.2)', dumr
      supstrc[i].fuse_snr = dumr
      ; STIS
      readf, 1,  format='(i3)', dumi
      supstrc[i].flg_stis = dumi   
      readf, 1,  format='(a60)', dumc
      supstrc[i].stis_comm[0] = dumc
      readf, 1,  format='(a60)', dumc
      supstrc[i].stis_comm[1] = dumc
      readf, 1,  format='(a60)', dumc
      supstrc[i].stis_comm[2] = dumc
      readf, 1,  format='(a60)', dumc
      supstrc[i].stis_comm[3] = dumc
      ; GHRS
      readf, 1,  format='(i3)', dumi
      supstrc[i].flg_ghrs = dumi   
      readf, 1,  format='(a60)', dumc
      supstrc[i].ghrs_comm[0] = dumc
      readf, 1,  format='(a60)', dumc
      supstrc[i].ghrs_comm[1] = dumc
      readf, 1,  format='(a60)', dumc
      supstrc[i].ghrs_comm[2] = dumc
      readf, 1,  format='(a60)', dumc
      supstrc[i].ghrs_comm[3] = dumc
      ;; CLOSE
      close, 1
  endfor

  return
end
