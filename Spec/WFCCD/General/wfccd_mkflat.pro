;+ 
; NAME:
; wfccd_mkflat   
;   Version 1.0
;
; PURPOSE:
;    Creates Bias frame given structure
;
; CALLING SEQUENCE:
;   
;   wfccd_mkflat, struct, mask_id, SVOV=svov
;
; INPUTS:
;   struct -- wfccd_strct defining the images of interest
;
; RETURNS:
;
; OUTPUTS:
;   flat - fits file in the dir Flats named 'Flat_##.fits'
;                 where ## is the mask_id value
;   VAR  - Variance in the flat (in electrons)
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
;   wfccd_mkflat, nght1_strct, mask_id
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_mkflat, struct, mask_id, SVOV=svov, VAR=var, FLAT=comb, NOFITS=nofits

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'wfccd_mkflat, struct, mask_id, SVOV= [v1.0]'
      return
  endif 
  
;  Optional Keywords
  
;  Find the Flat frames

  flts = where(struct.type EQ 'FLT' AND struct.flg_anly NE 0 AND $
               struct.mask_id EQ mask_id, nflts)

  if nflts LE 0 then begin
      print, 'wfccd_mkflat: No flat files available!'
      stop
      return
  endif

; Set gain and RN
  gain = struct[flts[0]].gain
  rn = struct[flts[0]].readno

;  Overscan

  ovflt = where(struct[flts].flg_ov EQ 0, nov)
  if nov NE 0 then wfccd_over, struct, flts[ovflt], ORDR=11

; Create directory if necessary and set out file

  a = findfile('Flats/..', count=count)
  if count EQ 0 then file_mkdir, 'Flats'

  if mask_id LT 10 then outfil = 'Flats/Flat_0'+$
    string(mask_id,format='(i1)')+'.fits' $
  else outfil = 'Flats/Flat_'+string(mask_id,format='(i2)')+'.fits'

; Status
  if not keyword_set( SILENT ) then print, 'wfccd_mkflat: Combining images: '

; Combine Images
  if nflts GT 2 then begin
      xcombine, struct[flts].img_ov, comb, head  ; Median
      var = gain*(comb>0)/float(nflts) + RN^2
  endif else begin 
      head = headfits(struct[flts[0]].img_ov, /silent)
      if nflts EQ 2 then comb = x_addtwoflats(struct[flts[0]].img_ov, $
                                              struct[flts[1]].img_ov, $
                                             GAIN=gain, RN=rn, VAR=var $
                                             ) $
      else begin
          comb = xmrdfits(struct[flts[0]].img_ov, /silent)
          var = (gain*(comb>0) + RN^2)
      endelse
  endelse

; Zero out bad columns
  if struct[flts[0]].ccd EQ 'WFTek5' then begin
      comb[1012:1013, 644:2046] = 0.
      var[1012:1013, 644:2046] = -1.
      comb[1011:1015, 623:716] = 0.
      var[1011:1015, 623:716] = -1.
  endif

; Delete the ov files
  if not keyword_set( SVOV ) then wfccd_delov, struct, flts

; Update structure
  struct[flts].flg_final = 1

  if not keyword_set( NOFITS ) then begin
; Update structure
      struct[flts].img_final = outfil
      struct[flts].flat_fil = outfil
      ;; Output
      mwrfits, comb, outfil, head, /create
      mwrfits, var, outfil, head
      ;; Compress
      spawn, 'gzip -f '+outfil
  endif

  print, 'wfccd_mkflat:  All done!'

end
