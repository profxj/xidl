;+ 
; NAME:
; x_chkfil   
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
;lines_strct, '/data4/SDSS/DR1_QSO/dr1_fit.lst', 'sdss_DR1_QAL.fits',
;  SDSSSUM='/data4/SDSS/DR1_QSO/summ_z04_R195.fit'
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   07-May-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro sdss_search, filename, outfil, ESIF=esif, SDSSSUM=sdsssum, ZMIN=zmin, RMAX=rmax

  if  N_params() LT 2 then begin 
      print,'Syntax - ' + $
        'find_lines, filename, outfil (v1.0)'
      return
  endif 
  
  if not keyword_set( ZMIN ) then zmin = 1.6
  if not keyword_set( RMAX ) then rmax = 19.5

;; List
  
  if keyword_set( ESIF ) then begin
      readcol, filename, list, zem, sigv, FORMAT='A,F'
  endif else begin
      readcol, filename, list, FORMAT='A'
      absdir = list[0]
      list = list[1:*]
      ;; Summary table
      if not keyword_set( SDSSSUM ) then stop
      sdsstab = xmrdfits(sdsssum, 1, /silent)
      ;; Parse
      gd = where(sdsstab.z GT zmin AND sdsstab.MAG_R LT Rmax, ngd)
      if ngd EQ 0 then stop else list = list[gd]
  endelse	

  nfil = n_elements(list)

;; ESIF
if keyword_set( ESIF ) then readcol, esif, fil_esif, format='A'

;  Make the structure
tmp = {qalcharstrct}
qalstrct= replicate(tmp, nfil)

; Loop

for qq=0L, nfil-1 do begin
    
    ;; Lya
    if keyword_set( ESIF ) then begin
        sdss_dla, zem=zem[qq], ESI=fil_esif[qq], sigv=sigv[qq], STRCT=strct
        qalstrct[qq] = strct
    endif else begin
        sdss_dla, SDSS=list[qq], STRCT=strct
        qalstrct[qq] = strct
    endelse

    print, 'lines_strct: Analysing file ', list[qq], qq

    ;; Look for metal lines
    if keyword_set( ESIF) then $
      findlines, list[qq], zem[qq], abswave, QSTRCT=strct $ 
    else begin
        parse_sdss, list[qq], flux, wave, SIG=sig, ZQSO=zqso
        absfil = absdir+qalstrct[qq].qso_name+'.fits'
        if x_chkfil(absfil+'*', /silent) EQ 0 then stop
        abslin = xmrdfits(absfil, 1, /silent)
        sdss_metals, wave, 0., abslin, ZEM=zqso, QSTRCT=strct
    endelse

    
    
    ;; fill up the structure
    qalstrct[qq].ndla2 = strct.ndla2
    qalstrct[qq].dla_zabs2 = strct.DLA_zabs2
    qalstrct[qq].dla_score2 = strct.DLA_score2
    qalstrct[qq].dla_hits = strct.DLA_hits
     
;; compare z values and give quality
    sdss_compare, qalstrct[qq], FSTRCT=strct
    
;; fill structure with z and quality
    qalstrct[qq].dla_z=strct.dla_z
    qalstrct[qq].dla_quality=strct.dla_quality
    
    qalstrct[qq].file_name=list[qq]

    fndz = where(qalstrct[qq].dla_z NE 0., nz)
    if nz NE 0 then print, 'Found z: ', qalstrct[qq].dla_z[fndz]

endfor

if keyword_set(ESIF) then qalstrct.qso_name = ' '

;; Output
mwrfits, qalstrct, outfil, /create

return

end
