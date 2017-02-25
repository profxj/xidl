;+ 
; NAME:
; esi_echmkflat   
;     Version 1.1
;
; PURPOSE:
;    Identifies flats, processes them.  Creates one flat file
;      per slit width.  Takes Dome flats as the default
;
; CALLING SEQUENCE:
;   
;  esi_echmkflat, esi, /IFLAT
;
; INPUTS:
;   esi     -  ESI structure
;   [slit]  -  Slit size (e.g. 0.5, 0.75, 1.0)
;
; RETURNS:
;
; OUTPUTS:
;  One normalized flat per slit width
;
; OPTIONAL KEYWORDS:
;   /IFLAT   - Use internal flats
;   /REDOOV  - Overwrite OV files if they exist for the flats
;   /SVOV    - Save the OV files created during this step
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;  Only set for 1x1 binning
;
; EXAMPLES:
;   esi_echmkflat, esi, 1.0
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   04-Aug-2002 Written by JXP
;   01-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echmkflat, esi, slit, IFLAT = iflat , BIASFIL=biasfil $
                   , TWIFLAT = TWIFLAT, REDOOV = redoov, SVOV = svov $
                   , CLOBBER = clobber, NOHOT = NOHOT, HOT_THRESH = HOT_THRESH

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'esi_echmkflat, esi, [slit], /IFLAT, /REDOOV, /CLOBBER [v1.1]'
      return
  endif 
  
  ;;  Optional Keywords
  IF keyword_set(IFLAT) THEN ftype = 'IFLT' $
  ELSE IF KEYWORD_SET(TWIFLAT) THEN ftype = 'TWI' $
  ELSE ftype = 'DFLT'
  
  ;; Grab all ECH FLATS files
  aflt = where(esi.mode EQ 2 AND esi.flg_anly NE 0 AND $
              strtrim(esi.type,2) EQ ftype , nflt)
  if nflt EQ 0 then begin
      print, 'esi_echmkflat: No Flats of type ', ftype, ' found!' 
      if strmatch(ftype,'IFLT') then ftype = 'DFLT' else ftype = 'IFLT'
      print, 'esi_echmkflat: Trying to look for ', ftype
      aflt = where(esi.mode EQ 2 AND esi.flg_anly NE 0 AND $
                   strtrim(esi.type,2) EQ ftype , nflt)
      if nflt EQ 0 then begin
          print, 'esi_echmkflat: No Flats of type ', ftype, ' found!' 
          print, 'esi_echmkflat: You may need to set /TWIFLAT or /IFLAT'
          return
      endif else begin
          if strmatch(ftype,'IFLT') then IFLAT = 1
      endelse
  endif

; Loop on binning
  bintot = esi[aflt].cbin + 10*esi[aflt].rbin
  bintyp = bintot[uniq(bintot, sort(bintot))]
  nbin = n_elements(bintyp)

  for qq=0L,nbin-1 do begin
      rbin = bintyp[qq]/10
      cbin = bintyp[qq] - rbin*10
  
      flt = where(esi.mode EQ 2 AND esi.flg_anly NE 0 AND $
                   esi.rbin EQ rbin AND esi.cbin EQ cbin AND $
                   strtrim(esi.type,2) EQ ftype , nflt)
      if nflt EQ 0 then begin
          print, 'esi_echmkflat: No Flats of type ', ftype, ' found! Punting..' 
          return
      endif

      ;; Loop on Slit width
      if not keyword_set( SLIT ) then begin
          slit = uniq(esi[flt].slit, sort(esi[flt].slit))
          slit = esi[flt[slit]].slit
      endif
      nslit = n_elements(slit)
      
      
      for ii=0,nslit-1 do begin
          ;; Slit
          tslit = slit[ii]
          
          ;; Outfil
          IF ftype EQ 'IFLT' then mode = 0 $
          ELSE IF ftype EQ 'DFLT' THEN mode = 1 $
          ELSE IF ftype EQ 'TWI' THEN mode = 2
          outfil = esi_getfil('echflat_fil', mode, SLIT=tslit, $
                              cbin=cbin, rbin=rbin, /name)

          if x_chkfil(outfil+'*',/silent) NE 0 $
            and not keyword_set( CLOBBER ) then begin
              print, 'esi_echmkflat:  File exists ', outfil, '  Continuing..'
              continue
          endif
          
          ;; Grab the flats
          gdflt = where(esi.mode EQ 2 AND esi.flg_anly NE 0 AND $
                        esi.rbin EQ rbin AND esi.cbin EQ cbin AND $
                        esi.slit EQ tslit AND esi.type EQ ftype , ngdflt)
          if ngdflt EQ 0 then begin
              print, 'esi_echmkflat:  No flats of type ', ftype, $
                ' found for slit width', tslit
              print, 'esi_echmkflat:  Continuing...'
              continue
          endif

          ;; Bias Subtract
          if not keyword_set( REDOOV ) then $
            bias = where(esi[gdflt].flg_ov EQ 0, nbias) $
          else begin
              nbias = n_elements(gdflt)
              bias = lindgen(nbias)
          endelse
          if nbias NE 0 then esi_subbias, esi, gdflt[bias] $
            , HOT_THRESH = HOT_THRESH, NOHOT = NOHOT, BIASFIL=biasfil

          ;; Median Combine
          statcol = round([1080,1180L]/float(cbin))
          statrow = round([1900,2200L]/float(rbin))
          STATSEC='['+strtrim(statcol[0],2)+':'+strtrim(statcol[1],2)+','+$
            strtrim(statrow[0],2)+':'+strtrim(statrow[1],2)+']'
          if ngdflt GT 1 then begin
              xcombine, 'OV/ov_'+esi[gdflt].img_root, img_flat, head, $
                FCOMB=2, SCALE='MED', GAIN=esi[gdflt[0]].gain, $
                RN=esi[gdflt[0]].readno, STATSEC=statsec
          endif else img_flat = xmrdfits('OV/ov_'+esi[gdflt].img_root,0,$
                                         head, /silent)
          ;; Output
          mwrfits, img_flat, outfil, head, /create, /silent
          spawn, 'gzip -f '+outfil
          print, 'esi_echmkflat: Flat created ', outfil+'.gz.'
          
          ;; Set Obj+Std
          ;objstd = where(esi.mode EQ 2 AND esi.flg_anly NE 0 AND $
          ;               esi.rbin EQ rbin AND esi.cbin EQ cbin AND $
          ;               esi.slit EQ tslit AND $
          ;               (esi.type EQ 'STD' OR esi.type EQ 'OBJ'), nobj)
          ;IF NOT KEYWORD_SET(TWIFLAT) AND nobj NE 0 THEN $
          ;  esi[objstd].flat_fil = outfil
          
          ;; DEL OV
          if not keyword_set(SVOV) then esi_delov, esi, gdflt
          
      endfor
  endfor

  print, 'esi_echmkflat: All done!'
  return
end
