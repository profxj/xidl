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

pro esi_echmkflat, esi, slit, IFLAT=iflat, REDOOV=redoov, SVOV=svov

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'esi_echmkflat, esi, [slit], /IFLAT, /REDOOV [v1.1]'
      return
  endif 
  
;  Optional Keywords
  if keyword_set( IFLAT ) then ftype = 'IFLT' else ftype = 'DFLT'
  
; Grab all ECH FLATS files

  flt = where(esi.mode EQ 2 AND esi.flg_anly NE 0 AND $
              strtrim(esi.type,2) EQ ftype , nflt)
  if nflt EQ 0 then begin
      print, 'esi_echmkflat: No Flats of type ', ftype, ' found!' 
      return
  endif

; Loop on Slit width

  if not keyword_set( SLIT ) then begin
      slit = uniq(esi[flt].slit, sort(esi[flt].slit))
      slit = esi[flt[slit]].slit
  endif
  nslit = n_elements(slit)


  for ii=0,nslit-1 do begin
      ;; Slit
      tslit = slit[ii]
      c_s = esi_slitnm( tslit )

      ;; Grab the flats
      gdflt = where(esi.mode EQ 2 AND esi.flg_anly NE 0 AND $
                    esi.slit EQ tslit AND esi.type EQ ftype , ngdflt)
      if ngdflt EQ 0 then begin
          print, 'esi_echmkflat:  No flats of type ', ftype, ' found for slit width',$
            tslit
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
      if nbias NE 0 then esi_subbias, esi, gdflt[bias]

      ;; Median Combine
      statcol = round([1080,1180L]/float(esi[gdflt[0]].cbin))
      statrow = round([1900,2200L]/float(esi[gdflt[0]].rbin))
      STATSEC='['+strtrim(statcol[0],2)+':'+strtrim(statcol[1],2)+','+$
        strtrim(statrow[0],2)+':'+strtrim(statrow[1],2)+']'
      xcombine, 'OV/ov_'+esi[gdflt].img_root, img_flat, head, $
        FCOMB=2, SCALE='MED', GAIN=esi[gdflt[0]].gain, $
        RN=esi[gdflt[0]].readno, STATSEC=statsec
    
      ;; Output
      if ftype EQ 'IFLT' then outfil = 'Flats/FlatECH'+c_s+'_I.fits' $
      else outfil = 'Flats/FlatECH'+c_s+'_D.fits'
      mwrfits, img_flat, outfil, head, /create, /silent
      spawn, 'gzip -f '+outfil
      print, 'esi_echmkflat: Flat created ', outfil+'.gz.'

      ;; Set Obj+Std
      objstd = where(esi.mode EQ 2 AND esi.flg_anly NE 0 AND $
                     esi.slit EQ tslit AND $
                     (esi.type EQ 'STD' OR esi.type EQ 'OBJ'), nobj)
      if nobj NE 0 then esi[objstd].flat_fil = outfil

      ;; DEL OV
      if not keyword_set(SVOV) then esi_delov, esi, gdflt

  endfor

  print, 'esi_echmkflat: All done!'
  return
end
