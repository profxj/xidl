;+ 
; NAME:
; wfccd_fixarc
;    Version 1.0
;
; PURPOSE:
;    Solves arc solutions for a given mask
;      Designed to do only 1 at a time
;
; CALLING SEQUENCE:
;   
;   wfccd_fixarc, wfccd, WFARC=
;
; INPUTS:
;   wfstrct     - WFCCD structure
;
; RETURNS:
;
; OUTPUTS:
;   wfarc      -  WFCCD arc structure (fits file)
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_fixarc, wfstrct, mask_id, exp_id, slit
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   12-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------


pro wfccd_fixarc, wfccd, mask_id, exp_id, slit, DEBUG=debug, CENT=cent

;
  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
      'wfccd_fixarc, wfccd, mask_id, exp_id, slit [v1.0]'
    return
  endif 

  obj_id = (where(wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0 AND $
              wfccd.mask_id EQ mask_id, nobj))[exp_id]
  arcfil = wfccd[obj_id].arc_fil

  ;; WFARC
  ipos = strpos(wfccd[obj_id].arc_fil, 'R_')
  wfarc_fil = 'Arcs/ArcS_'+strmid(wfccd[obj_id].arc_fil,ipos+2)
  wfccd_readastrct, wfarc_fil, wfarc

;  Read in the line list
  linelist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/wfccdB_HeNe.lst'

;  Set slit structure
  slitstr = xmrdfits(wfccd[obj_id].slit_fil,1, STRUCTYP='mslitstrct', /silent)

;  Read the Arc
  arc = xmrdfits(wfccd[obj_id].arc_fil, /silent)

; Find slits
  gd = where(slitstr.flg_anly NE 0, ngd)
  nslit = n_elements(slitstr)


  if not keyword_set( SILENT ) then print, 'wfccd_fixarc: Fixing slit ', slit
  gdi = slit 
;  gdi = (where(gd EQ slit))[0]

  ;; Center of slit
  if not keyword_set( CENT ) then begin
      cent = total(slitstr[gd[gdi]].yedg_flt)/2.
      ;; 0 Slit kludge
      if gdi EQ 0 AND cent GT 1850L then cent = slitstr[gd[gdi]].yedg_flt[0]+5
      if gdi EQ (ngd-1) AND cent LT 100L then cent = slitstr[gd[gdi]].yedg_flt[1]-5
      wfarc[gd[gdi]].cent = cent
  endif

  ;; Take center 5 rows
  spec = djs_median(arc[*,cent-2:cent+2], 2)
  npix = n_elements(spec)

  if keyword_set( DEBUG ) then stop

  ;; Run x_identify
  x_identify, spec, calib, WAVE=wave, LINELIST=linelist, /redblue

  ;; Save fitstr
  wfarc[gd[gdi]].wave = wave
  wfarc[gd[gdi]].spec = temporary(spec)
  wfarc[gd[gdi]].fit = temporary(calib)
  wfarc[gd[gdi]].flg_anly = 1


  ;; Output fits image
  wfccd_writeastrct, wfarc, wfarc_fil
  if not keyword_set( SILENT ) then print, 'wfccd_fixarc: All done'

  print, 'wfccd_fixarc: Now do --  wfccd_chkarc, wfccd, '+strtrim(mask_id,2)+'L, '+$
    strtrim(exp_id,2)+'L'
  print, 'wfccd_fixarc: And    --  wfccd_arcimg, wfccd, '+strtrim(mask_id,2)+'L, '+$
    strtrim(exp_id,2)+'L, /CLOBBER'
  print, 'wfccd_fixarc: And    --  wfccd_cleanarc, wfccd, '+strtrim(mask_id,2)+'L, '+$ 
strtrim(exp_id,2)+'L'

  print, 'wfccd_fixarc: Also reprocess as necessary!!'
  
  return
end
