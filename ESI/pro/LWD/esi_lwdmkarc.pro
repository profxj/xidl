;+ 
; NAME:
; esi_lwdmkarc   
;     Version 1.0
;
; PURPOSE:
;    Process arc file
;
; CALLING SEQUENCE:
;   
;  esi_lwdmkarc, esi, slit
;
; INPUTS:
;   esi     -  ESI structure
;
; RETURNS:
;
; OUTPUTS:
;  One normalized flat per slit width
;
; OPTIONAL KEYWORDS:
;   DFLAT      - Use Dome flats where possible
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_lwdmkarc, esi
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   27-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_lwdmkarc, esi, slit, CLOBBER=clobber

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'esi_lwdmkarc, esi, slit, /CLOBBER [v1.0]'
      return
  endif 
  
;  Optional Keywords
  
  c_s = esi_slitnm(slit)

; Check for prior image
  outfil = 'Arcs/ArcLWD_'+c_s+'.fits'
  a = findfile(outfil, count=na)
  if na NE 0 AND not keyword_set( CLOBBER ) then begin
      print, 'esi_lwdmkarc: Arc ', outfil, ' exists.  Returning'
      return
  endif

; Grab all ECH Arc files

  ;; Take HgXe only; CuAr is too much
  arcs = where(esi.mode EQ 1 AND esi.flg_anly NE 0 AND $
              esi.slit EQ slit AND strtrim(esi.type,2) EQ 'ARC', narc)
	
  gdarc = where(esi[arcs].arclamp EQ 6, ngd)
  if NGD NE 0 then begin  ; Take arc
      ;; Reset gdarc
      gdarc = arcs[gdarc]
      ;; BIAS (FORCE)
      esi_subbias, esi, gdarc, /FORCE
      ;; Median Combine
      if narc GT 1 then begin
          xcombine, 'OV/ov_'+esi[gdarc].img_root, fin_arc, head, $
            FCOMB=2, SCALE=esi[gdarc].exp, GAIN=esi[gdarc[0]].gain, $
            RN=esi[gdarc[0]].readno
      endif else fin_arc = mrdfits('OV/ov_'+esi[gdarc].img_root, /silent)
      ;; DELOV
      esi_delov, esi, gdarc
  endif else begin    ; Combine Hg, Xe
      gd2 = where(esi[arcs].arclamp EQ 2, ngd2)
      gd4 = where(esi[arcs].arclamp EQ 4, ngd4)
      if ngd2 EQ 0 OR ngd4 EQ 0 then begin
          print, 'esi_echmkarc: Wrong arcs taken!'
          stop
          return
      endif
      gdarc = [arcs[gd2[0]], arcs[gd4[0]]]
      print, 'esi_echmkarc: Combining HgNe+Xe exposures: ', $
        esi[gdarc[0]].img_root, ' ', esi[gdarc[1]].img_root
      
      ;; BIAS (FORCE)
      esi_subbias, esi, gdarc, /FORCE

      ;; Open Arcs
      img2_fil = 'OV/ov_'+esi[gdarc[0]].img_root
      img2_arc = mrdfits(img2_fil, 0, head, /silent)
      sz_arc = size(img2_arc, /dimensions)
      img4_fil = 'OV/ov_'+esi[gdarc[1]].img_root
      img4_arc = mrdfits(img4_fil, /silent)

      ;; Create VAR image
      var2 = (img2_arc*esi[gdarc[0]].gain + esi[gdarc[0]].readno^2) > 0.
      var4 = (img4_arc*esi[gdarc[0]].gain + esi[gdarc[0]].readno^2) > 0.

      ;; Sum Arcs
      fin_arc = img2_arc + img4_arc
;      fin_var = var2 + var4


      ;; DELOV
      esi_delov, esi, gdarc
  endelse

  ;; Output
  mwrfits, fin_arc, outfil, head, /create, /silent
;      mwrfits, fin_var, outfil, /silent
  print, 'esi_echmkarc: Arc created ', outfil

  print, 'esi_echmkarc: All Done! '
  return
end
