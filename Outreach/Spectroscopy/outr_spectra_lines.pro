;+ 
; NAME:
; outr_spectra_lines
;   Version 1.1
;
; PURPOSE:
;    GUI for plotting a spectrum and doing quick analysis 
;
; CALLING SEQUENCE:
;   
;   outr_spectratool_ flux, [ysin], XSIZE=, YSIZE=, TITLE=, WAVE=, LLIST=,
;           /QAL, /GAL, ZIN=, /BLOCK, /NRM, /LLS
;
; INPUTS:
;   flux  - Flux array (or FITS file)
;   [ysin]  - Sigma array (or FITS file)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   xsize      - Draw window xsize (pixels)
;   ysize      - Draw window ysize (pixels)
;   wave=      - wavelength array
;   INFLG=     - Specifies the type of input files (or arrays)
;   /LLS       - Use Lyman limit line list
;   /QSO       - Use Quasar line list
;   /GAL       - Use galaxy line list
;   /QAL       - Use quasar absorption line list
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   outr_spectratool_ 'spec.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  XGETX_PLT
;  XGETY_PLT
;  XGETXPIX_PLT
;  XGETYPIX_PLT
;
; REVISION HISTORY:
;   22-Apr-2008 Written by JXP
;-
;------------------------------------------------------------------------------
pro outr_spectra_lines_modify, pdv, lines

  ;; 
  sprs = strsplit(pdv, '.', /extract, count=nsplit)
  if strmatch(sprs[nsplit-1], 'None') then begin
      lines.j = -1
      return
  endif

  ;; Emission or Abs?
  if strmatch(sprs[1],'Emission') then gd = where(lines.A LT 0) $
  else gd = where(lines.A GT 0) 

  ;;
  mt = where(strmatch(lines[gd].ion, sprs[nsplit-1]))
  lines[gd[mt]].j = lines[gd[mt]].j * (-1)

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro outr_spectra_lines, pd, lines

  ;; Initialize

  ;; Create the string array
  pd = ['1\Spectral Lines', $
        '0\None', $
        '1\Emission']

  ;; Read in emission lines
  emiss_fil = getenv('XIDL_DIR')+'Outreach/Spectroscopy/Lines/emission.lst'
  readcol, emiss_fil, wv, elm, format='D,A'
  e_lines = replicate({newabslinstrct}, n_elements(wv))
  e_lines.wrest = wv
  e_lines.ion = elm
  e_lines.A = -1.

  ;; Add on emission lines
  elem = x_uniqstr(e_lines.ion, count=nelem, /sort)
  for ww=0L,nelem-1 do begin
      if ww NE nelem-1 then tmps = '0\'+elem[ww] else  tmps = '2\'+elem[ww] 
      pd = [pd, tmps]
  endfor

  ;; Read in galaxy abs lines
  gala_fil = getenv('XIDL_DIR')+'Outreach/Spectroscopy/Lines/galx_abs.lst'
  readcol, gala_fil, wv, elm, format='D,A'
  ga_lines = replicate({newabslinstrct}, n_elements(wv))
  ga_lines.wrest = wv
  ga_lines.ion = elm
  ga_lines.A = 1.

  ;; Add on absorption lines
  pd = [pd, '1\Absorption']
  elem = x_uniqstr(ga_lines.ion, count=nelem, /sort)
  for ww=0L,nelem-1 do begin
      if ww NE nelem-1 then tmps = '0\'+elem[ww] else  tmps = '2\'+elem[ww] 
      pd = [pd, tmps]
  endfor

  ;; Lines
  lines = [e_lines, ga_lines]
  lines.j = -1.

return
end
