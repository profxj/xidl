;+ 
; NAME:
; esi_echmkarc   
;     Version 1.1
;
; PURPOSE:
;    Process and combine arc files  
;
; CALLING SEQUENCE:
;   
;  esi_echmkarc, esi, slit, /CLOBBER
;
; INPUTS:
;   esi     -  ESI structure
;   slit    -  Slit width  (0.5, 0.75, 1.0)
;
; RETURNS:
;
; OUTPUTS:
;  One processed, combined Arc Image (e.g. Arcs/Arc_ECH##.fits)
;
; OPTIONAL KEYWORDS:
;   /CLOBBER - Overwrite exisiting arc image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_echmkarc, esi, 0.5
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   12-Aug-2002 Written by JXP
;   01-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echmkarc, esi, slit, CLOBBER=clobber, FLATFIL=flatfil

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echmkarc, esi, slit, FLATFIL=, /CLOBBER [v1.1]'
      return
  endif 
  
;  Optional Keywords
  
  c_s = esi_slitnm(slit)

; Grab all ECH Arc files

  arcs = where(esi.mode EQ 2 AND esi.flg_anly NE 0 AND $
              esi.slit EQ slit AND strtrim(esi.type,2) EQ 'ARC', narc)
  if narc EQ 0 then begin
      print, 'esi_echmkarc: No Arcs found! Returning' 
      return
  endif

; Check for prior image
  outfil = 'Arcs/ArcECH_'+c_s+'.fits'
  a = findfile(outfil, count=na)
  if na NE 0 AND not keyword_set( CLOBBER ) then begin
      print, 'esi_echmkarc: Arc ', outfil, ' exists.  Returning'
      return
  endif

; Open Flat

  if not keyword_set( FLATFIL ) then $
    flatfil = 'Flats/FlatECH'+c_s+'N.fits'
  if x_chkfil(flatfil+'*') EQ 0 then begin
      print, 'esi_echmkarc: Flat does not exist!', flatfil, '  Returning...'
      return
  endif
  fhead = xheadfits(flatfil)
  print, 'esi_echmkarc: Checking the flat: ', flatfil
  scatt = sxpar(fhead,'SCATTER')
  norm = sxpar(fhead,'NORM')
  if scatt NE 1 OR norm NE 1 then begin
      print, 'esi_echmkarc: Flat not processed! Returning...'
      return
  endif

; Add all arc lamps

  ;; Create the TMP Arc images
  flg = bytarr(8)
  for i=1L,7 do begin
      if i EQ 3 or i EQ 5 then continue

      gd = where(esi[arcs].arclamp EQ i, ngd)
      if ngd EQ 0 then continue else flg[i] = 1B
      
      gdarc = arcs[gd]

      ;; Combine
      esi_echcombarc, esi, gdarc, FLATFIL=flatfil
  endfor

  ;; Add em up
  fin_arc = fltarr(2048L,4096L)
  fin_var = fltarr(2048L,4096L)
  for i=1L,7 do begin
      if flg[i] NE 1B then continue
      ;; Read image
      flg_lmp = strtrim(i,2)
      tmp_fil = 'Arcs/ATMP_'+c_s+'_'+flg_lmp+'.fits'
      img = xmrdfits(tmp_fil, 0, head, /silent)
      var = xmrdfits(tmp_fil, 1, /silent)
      ;; Add em up
      fin_arc = fin_arc + img
      fin_var = fin_var + var
      ;; Delete TMP file
      spawn, '\rm '+tmp_fil
  endfor

  ;; Output
  mwrfits, fin_arc, outfil, head, /create, /silent
  mwrfits, fin_var, outfil, /silent
      
  ;; Cards
  objstd = where(esi.mode EQ 2 AND esi.flg_anly NE 0 AND $
                 esi.slit EQ slit AND $
                 (esi.type EQ 'STD' OR esi.type EQ 'OBJ'), nobj)
  if nobj NE 0 then esi[objstd].arc_fil = outfil

  print, 'esi_echmkarc: All Done! '
  return
end

