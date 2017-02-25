;+ 
; NAME:
; sdss_mgiisearch   
;    Version 1.0
;
; PURPOSE:
;    Given a list of SDSS quasars, search these spectra for MgII
;    absorbers and output a structure of detections and EWs
;
; CALLING SEQUENCE:
;
; INPUTS:
;   filename - File containing a (long) list of SDSS quasar data files
;              (1d spectra, .fit)
;   SDSSSUM= - Summary file of QSO properties [required]
;   INIFIL=  - File containing the MgII structure which was created
;              previously. These will not be redone.
;
; RETURNS:
;
; OUTPUTS:
;  outfil  - Name for outputted structure of MgII absorbers
;
; OPTIONAL KEYWORDS:
; /LZ -- Restrict to low redshift (z = 0.45)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  sdss_mgiisearch, 'Lists/dr6_qso.list', 'DR6_QSO/', 'dr6_qso.fits', 'mgii_search.fits'
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   07-May-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro sdss_mgiisearch, filename, drroot, sdsssum, outfil, $
                     ZMIN=zmin, RMAX=rmax,LZ=lz, INIFIL=inifil

  if  N_params() LT 4 then begin 
      print,'Syntax - ' + $
        'sdss_mgiisearch, filename, drroot, sdsssum, outfil (v1.0)'
      return
  endif 
  
  if not keyword_set( ZMIN ) then zmin = 0.35
  if not keyword_set( RMAX ) then rmax = 19.5


;; List
  
  readcol, filename, list, FORMAT='A'
  absdir = getenv('SDSSPATH')+DRROOT+list[0]
  list = list[1:*]
  ;; Summary table
  if not keyword_set( SDSSSUM ) then stop
  sdsstab = xmrdfits(sdsssum, 1, /silent)
  ;; Parse
  gd = where(sdsstab.z GT zmin AND sdsstab.psf_R LT Rmax, ngd)
  if ngd EQ 0 then stop else list = list[gd]

  nfil = n_elements(list)


  ;;  Make the structure
  tmp = {qalcharstrct}
  qalstrct= replicate(tmp, nfil)

  flg_mgii = 0L

  if keyword_set( INIFIL ) then begin
      tmp2 = xmrdfits(inifil,1,/silent)
      allmgii = xmrdfits(inifil,2,/silent)
      b = where(tmp2.z_qso GT 0., nb)
      istrt = nb - (nb MOD 500)
      qalstrct = tmp2
      flg_mgii = 1L
  endif else ISTRT=0L


  print, 'Nfil = ', nfil

;  for qq=0L, 10-1 do begin
  for qq=ISTRT, nfil-1 do begin
    
      print, 'lines_strct: Analysing file ', list[qq], qq
      
      ;; Grab names
      fname = getenv('SDSSPATH')+DRROOT+list[qq]
      parse_sdss, fname, flux, wave, SIG=sig, ZQSO=zqso, HEAD=head
      if keyword_set(LZ) then begin
        if zqso GE 0.45 then continue
      endif
      nm = sxpar(head, 'NAME')
      fib= sxpar(head, 'FIBERID')
      plate= sxpar(head, 'PLATEID')
      qalstrct[qq].fiberid=fib
      qalstrct[qq].plate=plate
      qalstrct[qq].qso_name=strtrim(nm,2)+'-'+strtrim(fib,2)
      RA = sxpar(head, 'RAOBJ')
      qalstrct[qq].RA= RA
      DEC = sxpar(head, 'DECOBJ')
      qalstrct[qq].DEC= DEC
      mag = sxpar(head, 'MAG_R')
      qalstrct[qq].qso_mag=mag
      qalstrct[qq].z_qso=zqso
      
      absfil = absdir+qalstrct[qq].qso_name+'.fits'
      if x_chkfil(absfil+'*', /silent) EQ 0 then stop
      abslin = xmrdfits(absfil, 1, /silent)
      sdss_mgii, abslin, ZEM=zqso, QSTRCT=strct
      
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
      if flg_mgii NE 0 then begin
          a = where( abs(qalstrct[qq].ra-allmgii.ra) LT 0.001 AND $
                     abs(qalstrct[qq].dec-allmgii.dec) LT 0.001 AND $
                     abs(qalstrct[qq].z_qso-allmgii.z_qso) LT 0.1, na)
          if na NE 0 then flg_dup = 1 
      endif
          
      if flg_dup EQ 1 then begin
          b = where(strlen(allmgii[a].sdss_obs) NE 0, nb)
          if nb EQ 0 then stop
          allmgii[a].sdss_obs[b[0]] = list[qq]
      endif else begin
          ;; Fill up EW
          for kk=0L,strct.ndla2-1 do begin
              sdss_ewmgii, wave, flux, sig, conti, mgiistr, $
                ZABS=strct.dla_zabs2[kk]
              
              mgiistr.qso_name = qalstrct[qq].qso_name
              mgiistr.ra = qalstrct[qq].ra
              mgiistr.dec = qalstrct[qq].dec
              mgiistr.plate = qalstrct[qq].plate
              mgiistr.fiber = qalstrct[qq].fiberid
              mgiistr.rmag = qalstrct[qq].qso_mag
              mgiistr.z_qso = qalstrct[qq].z_qso
              mgiistr.sdss_obs[0] = list[qq]
              mgiistr.zabs = qalstrct[qq].dla_zabs2[kk]
              
              ;; Structure
              if flg_mgii EQ 0 then begin
                  allmgii = mgiistr 
                  flg_mgii = 1
              endif else allmgii = [allmgii, mgiistr]
          endfor
      endelse
    ;; Output
    if (qq + 1) MOD 500 EQ 0 then begin
        mwrfits, qalstrct, outfil, /create
        mwrfits, allmgii, outfil
    endif

  endfor

  ;; Output
  mwrfits, qalstrct, outfil, /create
  mwrfits, allmgii, outfil

  return

end
