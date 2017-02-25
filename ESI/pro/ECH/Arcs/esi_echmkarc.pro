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
;   29-Apr-2007 Removed option to flat field arc (JFH)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echmkarc, esi, slit, CLOBBER=clobber, BIASFIL=biasfil $
                  , SEDG_FIL = SEDG_FIL, FLATFIL=flatfil

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echmkarc, esi, slit, /CLOBBER [v1.1]'
      return
  endif 
  
;  Optional Keywords
  
;; Binning

  aarcs = where(esi.mode EQ 2 AND esi.flg_anly NE 0 AND $
              esi.slit EQ slit AND strtrim(esi.type,2) EQ 'ARC', narc)
  if narc EQ 0 then begin
      print, 'esi_echmkarc: No Arcs found! Returning' 
      return
  endif
  bintot = esi[aarcs].cbin + 10*esi[aarcs].rbin
  bintyp = bintot[uniq(bintot, sort(bintot))]
  nbin = n_elements(bintyp)

  for qq=0L,nbin-1 do begin
      rbin = bintyp[qq]/10
      cbin = bintyp[qq] - rbin*10
  
; Grab all ECH Arc files
      arcs = where(esi.mode EQ 2 AND esi.flg_anly NE 0 AND $
                   esi.rbin EQ rbin AND esi.cbin EQ cbin AND $
                   esi.slit EQ slit AND strtrim(esi.type,2) EQ 'ARC', narc)
      if narc EQ 0 then begin
          print, 'esi_echmkarc: No Arcs found for this binning! Continuing' 
          continue
      endif
      
; Check for prior image
      outfil = esi_getfil('arc_fil', SLIT=slit, $
                          cbin=cbin, rbin=rbin, /name)
      a = findfile(outfil+'*', count=na)
      if na NE 0 AND not keyword_set( CLOBBER ) then begin
          print, 'esi_echmkarc: Arc ', outfil, ' exists.  Returning'
          return
      endif
      if not keyword_set( SEDG_FIL ) then begin
          sedg_fil = esi_getfil('sedg_fil', SLIT = slit $
                                , cbin = cbin, rbin = rbin, /name)
          slit_edg = xmrdfits(sedg_fil, 0)
          tset_slits = xmrdfits(sedg_fil, 1)
      endif else begin
          slit_edg = xmrdfits(strtrim(SEDG_FIL, 2), 0, /silent)
          tset_slits = xmrdfits(SEDG_FIL, 1)
      endelse

; Open Flat
;      if not keyword_set( FLATFIL ) then $
;        flatfil = esi_getfil('finflat_fil', SLIT=slit, cbin=cbin, rbin=rbin, /name)
;      if x_chkfil(flatfil+'*') EQ 0 then begin
;          print, 'esi_echmkarc: Flat does not exist!', flatfil, '  Returning...'
;          return
;      endif
;      fhead = xheadfits(flatfil)
;      print, 'esi_echmkarc: Checking the flat: ', flatfil
;      scatt = sxpar(fhead,'SCATTER')
;      norm = sxpar(fhead,'NORM')
;      if scatt NE 1 OR norm NE 1 then begin
;          print, 'esi_echmkarc: Flat not processed! Returning...'
;          return
;      endif

; Add all arc lamps
      
      ;; Create the TMP Arc images
      flg = bytarr(8)
      for i=1L,7 do begin
          if i EQ 3 or i EQ 5 then continue
          gd = where(esi[arcs].arclamp EQ i, ngd)
          if ngd EQ 0 then continue else flg[i] = 1B
          
          gdarc = arcs[gd]
          ;; Combine
          esi_echcombarc, esi, gdarc, FLATFIL=flatfil, CLOBBER = CLOBBER, BIASFIL=biasfil
      endfor

      ;; Add em up
      fin_arc = fltarr(2048L/cbin,4096L/rbin)
      fin_var = fltarr(2048L/cbin,4096L/rbin)
      for i=1L,7 do begin
          if flg[i] NE 1B then continue
          ;; Read image
          flg_lmp = strtrim(i,2)
          c_s = esi_slitnm( slit )
          tmp_fil = 'Arcs/ATMP_'+c_s+'_'+flg_lmp+'.fits'
          img = xmrdfits(tmp_fil, 0, head, /silent)
          var = xmrdfits(tmp_fil, 1, /silent)
          ;; Add em up
          fin_arc = fin_arc + img
          fin_var = fin_var + var
          ;; Delete TMP file
          spawn, '\rm '+tmp_fil
      endfor

      ;; This step added by JFH 05/08
      ;plate_scale = reverse([0.168, 0.163, 0.158, 0.153, 0.149, 0.144, 0.137 $
      ;                       , 0.134, 0.127, 0.120])
      ;fwhm = slit/plate_scale
      ;plate_med = median(plate_scale)
      ;pkwdth = 1.3*slit/plate_med
      ;TOLER = pkwdth/3.0D
      ;mask = (fin_arc GT -20.0 AND fin_arc LT 1d5)
      ;wset = long_wavepix(mask*fin_arc, tset_slits, FWHM = FWHM $
      ;                    , pkwdth = pkwdth, toler = toler)
      ;; Output
      print, 'esi_echmkarc:  Creating ', outfil
      mwrfits, fin_arc, outfil, head, /create, /silent
      mwrfits, fin_var, outfil, /silent
;      mwrfits, wset, outfil, /silent
      ;; Cards
      ;objstd = where(esi.mode EQ 2 AND esi.flg_anly NE 0 AND $
      ;               esi.slit EQ slit AND $
      ;               (esi.type EQ 'STD' OR esi.type EQ 'OBJ'), nobj)
      ;if nobj NE 0 then esi[objstd].arc_fil = outfil
      
      print, 'esi_echmkarc: All Done! '
  endfor
  return
end

