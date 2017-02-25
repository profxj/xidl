;+ 
; NAME:
; wfccd_mkspecfits
;    Version 1.0
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   
;   indx = wfccd_getobjnm(wffspec, obj_nm)
;
; INPUTS:
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
;   indx = wfccd_getobjnm(wffspec, obj_nm)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_mkspecfits, wffil, ID, OUTDIR=outdir, ROOT=root

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'wfccd_mkspefits, wffil, ID, OUTDIR=, ROOT= [v1.0]'
    return
  endif 

  ;; Keywords
  if not keyword_set(OUTDIR) then outdir = './'
  if not keyword_set(ROOT) then root = 'PKS0405-123_G'

  pks0405 = xmrdfits(wffil, 1, /silent)
  obj = strtrim(pks0405.id,2)+ strtrim(pks0405.obj_id,2)

  nid = n_elements(id)
  for qq=0L,nid-1 do begin
      ;; Match?
      mtc = where(strtrim(id[qq],2) EQ obj, nmt)
      if nmt EQ 0 then begin
          print, 'wfccd_mkspecfits:  No match to Object ', id[qq]
          continue
      endif

      ;; 
      mtc = mtc[0]
      gfil = where(strlen(strtrim(pks0405[mtc].fspec_fil,2)) GT 0, ngdf)
      if ngdf EQ 0 then begin
          print, 'wfccd_mkspecfits:  I believe no spectrum exists for ', $
                 id[qq]
          continue
      endif
      fil = strtrim(pks0405[mtc].fspec_fil[gfil[0]],2)
      wfccd_wrfspec, wffspec, fil, /read
      indx = x_getobjnm(wffspec, id[qq])

      ;; Data
      spec_wv = wffspec[indx].wave[0:wffspec[indx].npix-1]
      spec_fx = wffspec[indx].fx[0:wffspec[indx].npix-1]
      spec_sig = fltarr(n_elements(spec_wv))
      a = where(wffspec[indx].var[0:wffspec[indx].npix-1] GT 0.)
      spec_sig[a] = float(sqrt(wffspec[indx].var[a]))

      ;; Outfil
      outfil = OUTDIR+ROOT+id[qq]+'.fits'
      mwrfits, spec_fx, outfil, /create
      mwrfits, spec_sig, outfil
      mwrfits, spec_wv, outfil

  endfor

  return
end
