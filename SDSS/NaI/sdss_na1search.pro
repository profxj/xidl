;+ 
; NAME:
; sdss_na1search   
;    Version 1.0
;
; PURPOSE:
;    Given an array of strings, return an array of unique values +
;    number
;
; CALLING SEQUENCE:
;   
; uniq = x_chkfil(strings, COUNT=count)
;
; INPUTS:
;   strings - Array of strings
;
; RETURNS:
;   uniq  - Array of unique members
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   sort - Sort the output
;
; OPTIONAL OUTPUTS:
;  COUNT - number of unique strings
;
; COMMENTS:
;
; EXAMPLES:
;   flg = x_chkfil( lbls, count=count)
;
; sdss_na1search, 'Lists/dr3_abs.lst', 'sdss_DR3_NaI.fits',
;  SDSSSUM='dr3_qso.fits', DATR='DR3_QSO/'
;
;PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   07-May-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro sdss_na1search, filename, outfil, SDSSSUM=sdsssum, ZMIN=zmin, $
                    RMAX=rmax,LZ=lz, DATR=datr, INIFIL=inifil

  if  N_params() LT 2 then begin 
      print,'Syntax - ' + $
        'sdss_na1search, filename, outfil (v1.0)'
      return
  endif 
  
  if not keyword_set( ZMAX ) then zmax = 3.9
  if not keyword_set( ZMIN ) then zmin = 0.
  if not keyword_set( RMAX ) then rmax = 23.
  if not keyword_set( SDSSPATH ) then sdsspath = 'SDSSPATH'   


;; List
  
  readcol, filename, list, FORMAT='A'
  absdir = list[0]
  if keyword_set( DATR ) then begin
      list = getenv(SDSSPATH)+DATR+list[1:*] 
  endif else begin
      list = getenv(SDSSPATH)+list[1:*] 
  endelse
  npth = strlen( getenv(SDSSPATH) ) 

  ;; Summary table
  if not keyword_set( SDSSSUM ) then stop
  sdsstab = xmrdfits(sdsssum, 1, /silent)
  if n_elements(list) NE n_elements(SDSSTAB) then stop
  idx = lindgen(n_elements(list))
  nlst = n_elements(list)
  ;; Parse
  if tag_exist(sdsstab, 'PSF_R') then $
    gd = where(sdsstab.z GT zmin AND sdsstab.PSF_R LT Rmax AND $
               idx LE nlst and sdsstab.z LT zmax, ngd) $
  else $
    gd = where(sdsstab.z GT zmin AND sdsstab.MAG_R LT Rmax AND $
               idx LE nlst and sdsstab.z LT zmax, ngd) 
  if ngd EQ 0 then stop else list = list[gd]
  nfil = n_elements(list)


  ;;  Make the structure
  tmp = {qalcharstrct}
  qalstrct= replicate(tmp, nfil)

  flg_na1 = 0L

  if keyword_set( INIFIL ) then begin
      tmp2 = xmrdfits(inifil,1,/silent)
      b = where(tmp2.z_qso GT 0., nb)
      istrt = nb - (nb MOD 500)
      qalstrct = tmp2
      id = max(where(qalstrct[0:istrt-1].ndla2 GT 0))
      allna1 = xmrdfits(inifil,2,/silent)
      ;; Match
      mtch = where( abs(allna1.ra-qalstrct[id].ra) LT 1e-4 AND $
                    abs(allna1.dec-qalstrct[id].dec) LT 1e-4 AND $
                    abs(allna1.zabs- $
                        qalstrct[id].dla_zabs2[qalstrct[id].ndla2-1]) $ 
                    LT 1e-3, nm)
      if nm NE 1 then stop
      allna1 = allna1[0:mtch]
  endif else ISTRT=0L

;  for qq=0L, 10-1 do begin
  for qq=ISTRT, nfil-1 do begin
    
      print, 'lines_strct: Analysing file ', list[qq], qq
      ;; Magnitude
      if tag_exist(sdsstab, 'PSF_R') then mag = sdsstab[gd[qq]].PSF_R $
      else mag = sdsstab[gd[qq]].MAG_R
      
      ;; Grab names
;      fname = '/Users/prochter'+strmid(list[qq],6)
      parse_sdss, list[qq], flux, wave, SIG=sig, ZQSO=zqso, HEAD=head
;      if keyword_set(LZ) then begin
;        if zqso GE 0.45 then continue
;      endif
      nm = sxpar(head, 'NAME')
      fib= sxpar(head, 'FIBERID')
      qalstrct[qq].qso_name=strtrim(nm,2)+'-'+strtrim(fib,2)
      RA = sxpar(head, 'RAOBJ')
      qalstrct[qq].RA= RA
      DEC = sxpar(head, 'DECOBJ')
      qalstrct[qq].DEC= DEC
;      mag = sxpar(head, 'MAG_R')
      qalstrct[qq].qso_mag=mag
      qalstrct[qq].z_qso=zqso
      
      absfil = absdir+qalstrct[qq].qso_name+'.fits'
      if x_chkfil(absfil+'*', /silent) EQ 0 then stop
      abslin = xmrdfits(absfil, 1, /silent)
      sdss_fndna1, abslin, ZEM=zqso, QSTRCT=strct
      
      ;; fill up the structure
      qalstrct[qq].start_wave = strct.start_wave
      qalstrct[qq].dla_quality[0] = strct.dla_quality[0]
      qalstrct[qq].ndla2 = strct.ndla2
      qalstrct[qq].DLA_zabs2 = strct.dla_zabs2
      qalstrct[qq].dla_hits = strct.dla_hits
      
      ;; file name
      qalstrct[qq].file_name=list[qq]

      if strct.ndla2 NE 0 then conti = xmrdfits(absfil, 0, /silent)

      ;; Check for duplicate
      flg_dup = 0
      if flg_na1 NE 0 then begin
          a = where( abs(qalstrct[qq].ra-allna1.ra) LT 0.001 AND $
                     abs(qalstrct[qq].dec-allna1.dec) LT 0.001 AND $
                     abs(qalstrct[qq].z_qso-allna1.z_qso) LT 0.1, na)
          if na NE 0 then flg_dup = 1 
      endif
          
      if flg_dup EQ 1 then begin
          b = where(strlen(allna1[a].sdss_obs) NE 0, nb)
          if nb EQ 0 then stop
          allna1[a].sdss_obs[b[0]] = list[qq]
      endif else begin
          ;; Fill up EW
          for kk=0L,strct.ndla2-1 do begin
              sdss_ewna1, wave, flux/conti, sig/conti, na1str, $
                ZABS=strct.dla_zabs2[kk]
              
              stop
              na1str.qso_name = qalstrct[qq].qso_name
              na1str.ra = qalstrct[qq].ra
              na1str.dec = qalstrct[qq].dec
              na1str.rmag = qalstrct[qq].qso_mag
              na1str.z_qso = qalstrct[qq].z_qso
              ;; Strip off path
              na1str.sdss_obs[0] = strmid(list[qq],npth)
              na1str.zabs = qalstrct[qq].dla_zabs2[kk]
              
              ;; Structure
              if flg_na1 EQ 0 then begin
                  allna1 = na1str 
                  flg_na1 = 1
              endif else allna1 = [allna1, na1str]
          endfor
      endelse
    ;; Output
    if (qq + 1) MOD 500 EQ 0 then begin
        mwrfits, qalstrct, outfil, /create
        mwrfits, allna1, outfil
    endif
  endfor

  ;; Output
  mwrfits, qalstrct, outfil, /create
  mwrfits, allna1, outfil

  return

end
