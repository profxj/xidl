;+ 
; NAME:
; wfccd_procobj   
;  Version 1.0
;
; PURPOSE:
;    Process the object frames -- Outputs into Final
;       OV subtracts, Divides by normalized flat, creates VAR, adds in
;       Arc image
;
; CALLING SEQUENCE:
;   
;   wfccd_procobj, wfccd, mask_id
;
; INPUTS:
;   wfccd -- wfccd_strct defining the images of interest
;
; RETURNS:
;
; OUTPUTS:
;   finimg - fits file in the dir Final named: final_###.fits
;   VAR    - Variance in the obj (in electrons)
;   AIMG   - Writes the arc image in too
;
; OPTIONAL KEYWORDS:
;   SVOV - save ov files
;   NOFITS - No FITS output
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_procobj, wfccd, mask_id
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   21-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_procobj, wfccd, mask_id, SVOV=svov, CLOBBER=clobber

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'wfccd_procobj, struct, mask_id, /SVOV, /CLOBBER [v1.0]'
      return
  endif 
  
;  Optional Keywords
  
; Subsets
  obj = where(wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0 AND $
              wfccd.mask_id EQ mask_id, nobj)

;;;;;;;;;;;;
;  LOOP

  for qq=0L,nobj-1 do begin
      ;  Filenames
      finfil = 'Final/f_'+strtrim(wfccd[obj[qq]].img_root,2)
      a = findfile(finfil, count=count)
      if count GT 0 AND not keyword_set( CLOBBER ) then continue
      wfccd[obj[qq]].img_final = finfil

      ; OV Subtract
      if  wfccd[obj[qq]].flg_ov EQ 0 then $
        wfccd_over, wfccd, obj[qq], IMG=objimg, /NOFITS
      if not keyword_set( OBJIMG ) then $
        objimg = xmrdfits(wfccd[obj[qq]].img_ov, /silent)

      ; Create VAR image
      var = (objimg*wfccd[obj[qq]].gain + wfccd[obj[qq]].readno^2) > 0.

      ; Flatten
      if not keyword_set( SILENT ) then $
        print, 'wfccd_procobj: Flattening...'
      flat = xmrdfits(wfccd[obj[qq]].flat_fil, /silent)
      nwobj = x_normspec(objimg, flat, var, nrmvar)

      ; Read in Arc
      ii = strpos(wfccd[obj[qq]].arc_fil, '_')
      arcfil = 'Arcs/ArcI'+strmid(wfccd[obj[qq]].arc_fil,ii)
      arc = xmrdfits(strtrim(arcfil,2), /silent)

      ; Eliminate bad Arc pixels
      badpix = where(arc LE 0., nbad)
      if nbad NE 0 then begin
          nwobj[badpix] = 0.
          nrmvar[badpix] = 0.
      endif

      ; Write
      if not keyword_set( SILENT ) then $
        print, 'wfccd_procobj: Writing ', finfil
      mwrfits, nwobj, finfil, /create, /silent
      mwrfits, nrmvar, finfil, /silent
      mwrfits, arc, finfil, /silent
      ;; COMPRESS
      spawn, 'gzip -f '+finfil
  endfor

;;;;;;;;;;;;;;

  ; Delete OV files
  if not keyword_set( SVOV ) then wfccd_delov, wfccd, obj

  if not keyword_set( SILENT) then print, 'wfccd_procobj: All done!'

end
