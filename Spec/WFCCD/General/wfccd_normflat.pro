;+ 
; NAME:
; wfccd_normflat   
;   Version 1.0
;
; PURPOSE:
;    Flattens the flat
;
; CALLING SEQUENCE:
;   
;   wfccd_normflat, struct, mask_id
;
; INPUTS:
;   wfccd -- wfccd_strct defining the images of interest
;
; RETURNS:
;
; OUTPUTS:
;   flat - fits file in the dir Flats named 'Flat_##.fits'
;                 where ## is the mask_id value
;   VAR  - Variance in the flat (in electrons)
;
; OPTIONAL KEYWORDS:
;   SVOV - save ov files
;   NOFITS - No FITS output
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_normflat, nght1_strct, mask_id, /RMOLD
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_normflat, wfccd, mask_id, RMOLD=rmold, FLAT=flat, SLITSTR=slitstr, $
                    VAR=var

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'wfccd_normflat, struct, mask_id [v1.0]'
      return
  endif 
  
;  Optional Keywords
  
; Subsets
  obj = where(wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0 AND $
              wfccd.mask_id EQ mask_id, nobj)
  all = where(wfccd.type NE 'FLT' AND wfccd.flg_anly NE 0 AND $
              wfccd.mask_id EQ mask_id, nall)
  flts = where(wfccd.type EQ 'FLT' AND wfccd.flg_anly NE 0 AND $
              wfccd.mask_id EQ mask_id, nflt)
  if nflt EQ 0 then return

;  Filenames

  flat_fil = strtrim(wfccd[flts[0]].flat_fil,2)
  ; Add the N
  i = strpos(flat_fil, '_')
  outfil = strmid(flat_fil,0,i)+'N'+strmid(flat_fil,i)

; Input flat
  if not keyword_set( FLAT ) then flat = xmrdfits(flat_fil,0,/silent)
  if not keyword_set( VAR ) then var = xmrdfits(flat_fil,1,/silent)

; Slitstrct
  if not keyword_set( SLITSTR ) then $
    slitstr = xmrdfits(wfccd[obj[0]].slit_fil,1, STRUCTYP='mslitstrct', $
                       /silent)

; Call x_normflat
  print, 'wfccd_normflat: Normalizing the flat: ', flat_fil
  nflat = x_normflat(flat, slitstr, var, nrmvar, /silent, OUTNRM=nrm, /ORIG)

; Write to file
  mwrfits, nflat, outfil, /create
  mwrfits, nrmvar, outfil
  mwrfits, nrm, outfil
  ;; COMPRESS
  spawn, 'gzip -f '+outfil

; Update the Structure (except the Flats)
  wfccd[all].flat_fil = outfil

; Alld one
  print, 'wfccd_normflat: All done!'
  
  return
end
