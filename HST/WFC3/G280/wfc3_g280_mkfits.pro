;+ 
; NAME:
;  wfc3_g280_mkfits
;
; PURPOSE:
;   Writes the structure into a FITS file
;
; CALLING SEQUENCE:
;   
;  wfc3_g280_strct, FITSDIR=
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;   Creates a FITS file for each object in the structure
;
; OPTIONAL KEYWORDS:
;  FITSDIR= -- Directory containing the FITS file
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
; wfc3_g280_mkfits, wfc3_g280_strct, FITSDIR=fitsdir
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   10-Jun-2016 Written by MN
;------------------------------------------------------------------------------


pro wfc3_g280_mkfits, wfc3_g280_strct, FITSDIR=fitsdir

  ;; first create the fits name
  fitsname=fitsdir+wfc3_g280_strct.name+'_g280F.fits'
  
  ;; the structure to save
  gd=lindgen(wfc3_g280_strct.cnt)
  outstrct={FLAM:wfc3_g280_strct.flux(gd),$
            FLAM_SIG:wfc3_g280_strct.sig(gd),$
            SKY:wfc3_g280_strct.sky(gd),$
            WAVE:wfc3_g280_strct.wave(gd) $
           }
  gda=lindgen(wfc3_g280_strct.cnta)
  outstrcta={FLAM:reverse(wfc3_g280_strct.fluxa(gda)),$
             FLAM_SIG:reverse(wfc3_g280_strct.flux_siga(gda)),$
             SKY:reverse(wfc3_g280_strct.skya(gda)),$
             WAVE:reverse(wfc3_g280_strct.wavea(gda)) $
            }
  gdc=lindgen(wfc3_g280_strct.cntc)
  outstrctc={FLAM:wfc3_g280_strct.fluxc(gdc),$
             FLAM_SIG:wfc3_g280_strct.flux_sigc(gdc),$
             SKY:wfc3_g280_strct.skyc(gdc),$
             WAVE:wfc3_g280_strct.wavec(gdc) $
            }
     
  ;; the header
  sxaddpar,head,'QSO_NAME',wfc3_g280_strct.name,'Name of QSO'
  sxaddpar,head,'RA',wfc3_g280_strct.ra,'Right Ascension (decimal)'
  sxaddpar,head,'DEC',wfc3_g280_strct.dec,'Declination (decimal)'
  sxaddpar,head,'EXPTIME',wfc3_g280_strct.exptime,'Exposure Time (s)'

  ;; write the fits file
  mwrfits, outstrct, fitsname, head, /create, /silent
  mwrfits, outstrcta, fitsname, head, /silent
  mwrfits, outstrctc, fitsname, head, /silent
 
end
