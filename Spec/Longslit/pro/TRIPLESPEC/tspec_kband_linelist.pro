;+
; NAME:
;   tspec_kband_linelist
;
; PURPOSE:
;   Generate a linelist for TripleSpec K-band (last order)
;   Here is what I have done:
;       1.  Used R=2700 for TripleSpec and a 41km/s dispersion
;       2.  Monkeyed with the BlackBody until it looked like the data
;       3.  Scaled the OH and H2O lines too
;       4.  Am using the BB spectrum as the continuum /CONTI_BB
;       5.  Generate a K-band linelist with nearir_modelsky_linelist
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   ;
; REVISION HISTORY:
;   
;-
;------------------------------------------------------------------------------
;; tspec_kband_linelist, /redo
pro tspec_kband_linelist, SCALE=scale, REDO=redo, T_BB=t_bb, SCL_BB=scl_bb, $
                          SCL_OH=scl_oh, scl_h2o=scl_h2o

  if not keyword_set(T_BB) then T_BB  = 305.
  if not keyword_set(SCL_BB) then SCL_BB  = 4e2
  if not keyword_set(SCL_OH) then SCL_OH = 3e2
  if not keyword_set(SCL_H2O) then SCL_H2O = 1e3

  ;; Read archived sky
  restore, '../../calib/linelists/tspec_wave_archive.sav'
  wave = x_calcfit(dindgen(2048L), fitstr=calib[4])

  ;; Generate model
  outfil = 'TSPEC_Kband_modelsky.fits'
  if x_chkfil(outfil) and not keyword_set(REDO) then begin
     print, 'tspec_kband_linelist: Using existing file'
     print, 'tspec_kband_linelist: But you may still need to generate the linelist!'
  endif else $
     nearir_modelsky_linelist, 2700., 'TSPEC_Kband_linelist.lst', outfil, $
                               dlam=41.1, flgd=1, wvmnx=[1.9, 2.5], T_BB=T_BB, $
                               SCL_BB=scl_bb, scl_oh=scl_oh, scl_h2o=scl_h2o, /CONTI_BB ;, /nowrite

  ;; Read
  model_fx = x_readspec(outfil, wav=model_wv)
  ;x_splot, model_wv, model_fx, /blo

  ;; Plot
  x_splot, wave, archive_arc[*,4], xtwo=model_wv, ytwo=model_fx, /blo
  stop

  return
end
