;+ 
; NAME:
; sdss_search
;    Version 1.1
;
; PURPOSE:
;    Main program which drives the search for LLS, DLA and metal-line systems
;
; CALLING SEQUENCE:
;  sdss_search, filename, outfil, ESIF=, SDSSSUM=, ZMIN=, $
;                RMAX=, SDSSPATH=
;
; INPUTS:
;  filename -- ASCII file listing all SDSS quasars
;
; RETURNS:
;
; OUTPUTS:
;  outfil   -- Name of QAL file
;
; OPTIONAL KEYWORDS:
;  ZMIN=  -- Minimum redshift to search
;  RMAX=  -- Maximum magnitude to include [default: 21.]
;  INIFIL= -- Name of QAL file (good to restart after a 'bomb')
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   flg = x_chkfil( lbls, count=count)
;
; sdss_search, 'Lists/dr3_abs.lst', 'sdss_DR3_QAL.fits',
;  SDSSSUM='dr3_qso.fits'
;
; PROCEDURES/FUNCTIONS CALLED:
;  sdss_dla
;  sdss_metals
;
; REVISION HISTORY:
;   07-May-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro sdss_search, filename, outfil, ESIF=esif, SDSSSUM=sdsssum, ZMIN=zmin, $
                 RMAX=rmax, SDSSPATH=sdsspath, DATR=datr, INIFIL=inifil, $
                 NOMTL=nomtl, ZMAX=zmax, NODLA=nodla

  if  N_params() LT 2 then begin 
      print,'Syntax - ' + $
        'sdss_search, filename, outfil, ESIF=, SDSSSUM=, ZMIN=, RMAX= [v1.1]'
      return
  endif 
  
  if not keyword_set( ZMIN ) then zmin = 1.6
  if not keyword_set( ZMAX ) then zmax = 10.
  if not keyword_set( RMAX ) then rmax = 23.
  if not keyword_set( SDSSPATH ) then sdsspath = 'SDSSPATH'   

;; List
  
  if keyword_set( ESIF ) then begin
      readcol, filename, list, zem, sigv, FORMAT='A,F'
  endif else begin
      readcol, filename, list, FORMAT='A'
      absdir = list[0]
      if keyword_set( DATR ) then list = getenv(SDSSPATH)+DATR+list[1:*] $
      else list = getenv(SDSSPATH)+list[1:*] 
      ;; Summary table
      if not keyword_set( SDSSSUM ) then stop
      sdsstab = xmrdfits(sdsssum, 1, /silent)
      if n_elements(list) NE n_elements(SDSSTAB) then stop
      ;; 
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
  endelse	

  nfil = n_elements(list)

;; ESIF
if keyword_set( ESIF ) then readcol, esif, fil_esif, format='A'

;  Make the structure
tmp = {qalcharstrct}
qalstrct= replicate(tmp, nfil)
print, 'nfil = ', nfil

if keyword_set( INIFIL ) then begin
    tmp2 = xmrdfits(inifil,1,/silent)
    b = where(tmp2.z_qso GT 0., nb)
    istrt = nb - (nb MOD 500)
    qalstrct = tmp2
endif else ISTRT=0L

; Loop
for qq=ISTRT, nfil-1 do begin
    
    ;; Lya
    if not keyword_set( NODLA ) then begin
        if keyword_set( ESIF ) then begin
            sdss_dla, zem=zem[qq], ESI=fil_esif[qq], sigv=sigv[qq], STRCT=strct
            qalstrct[qq] = strct
        endif else begin
            ;; Magntiude
            if tag_exist(sdsstab, 'PSF_R') then mag = sdsstab[gd[qq]].PSF_R $
            else mag = sdsstab[gd[qq]].MAG_R
            ;; Search for DLA
            sdss_dla, SDSS=list[qq], STRCT=strct, MAG=mag
            ;; Add in plate, fiber, mjd
            strct.plate = sdsstab[gd[qq]].plate
            strct.fiberid = sdsstab[gd[qq]].fiberid
            qalstrct[qq] = strct
        endelse
    endif else begin
        qalstrct[qq] = {qalcharstrct}
        head = xheadfits(list[qq],/silent)
        name = sxpar(head, 'NAME')
        fiberid=strtrim(sxpar(head,'FIBERID'),2)
        qalstrct[qq].qso_name=name+'-'+fiberid
        qalstrct[qq].plate = sdsstab[gd[qq]].plate
        qalstrct[qq].fiberid = sdsstab[gd[qq]].fiberid
        qalstrct[qq].ra = sdsstab[gd[qq]].ra
        qalstrct[qq].dec = sdsstab[gd[qq]].dec
    endelse
        
        
    print, list[qq], qq

    ;; Look for metal lines
    flg_mtl = 0
    if not keyword_set( NOMTL ) then begin
        if keyword_set( ESIF) then begin
            findlines, list[qq], zem[qq], abswave, QSTRCT=strct 
            flg_mtl = 1
        endif else begin
;            parse_sdss, list[qq], flux, wave, SIG=sig, ZQSO=zqso
            zqso = sdsstab[gd[qq]].z
            absfil = absdir+qalstrct[qq].qso_name+'.fits'
            if x_chkfil(absfil+'*', /silent) EQ 0 then print, 'None' else begin
                abslin = xmrdfits(absfil, 1, /silent)
                sdss_metals, abslin, zqso, QSTRCT=strct, ZMIN=zmin
                flg_mtl = 1
            endelse
        endelse
    endif

    ;; fill up the structure
    if flg_mtl NE 0 then begin
        qalstrct[qq].ndla2 = strct.ndla2
        qalstrct[qq].dla_zabs2 = strct.DLA_zabs2
        qalstrct[qq].dla_score2 = strct.DLA_score2
        qalstrct[qq].dla_hits = strct.DLA_hits
    endif
     
;; compare z values and give quality
    sdss_compare, qalstrct[qq], FSTRCT=strct
    
;; fill structure with z and quality
    qalstrct[qq].dla_z=strct.dla_z
    qalstrct[qq].dla_quality=strct.dla_quality
    
    ipos = strpos(list[qq], 'DR')
    qalstrct[qq].file_name=strmid(list[qq],ipos)

    fndz = where(qalstrct[qq].dla_z NE 0., nz)
    if nz NE 0 then print, 'Found z: ', qalstrct[qq].dla_z[fndz]

    ;; Output
    if (qq + 1) MOD 500 EQ 0 then mwrfits, qalstrct, outfil, /create
endfor

if keyword_set(ESIF) then qalstrct.qso_name = ' '

;; Output
mwrfits, qalstrct, outfil, /create
mwrfits, findgen(100), 'dum.fil', /create

print, 'sdss_search: All done!!!'

return

end
