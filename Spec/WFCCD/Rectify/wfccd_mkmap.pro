;+ 
; NAME:
; wfccd_mkmap   
;   Version 1.0
;
; PURPOSE:
;    Creates a Map between distorted and undistorted frames
;
; CALLING SEQUENCE:
;   
;   wfccd_mkmap, struct, mask_id, FLAT=, VAR=, /NOCHK, MAP=
;
; INPUTS:
;   struct -- wfccd_strct defining the images of interest
;   mask_id -- Mask ID value (scalar)
;   [FLAT] - Flat (in memory)
;   [VAR]  - Variance of the Flat
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   NOFITS - No FITS output of the map
;   YSTRT  - Starting column for the tracing
;   OUTFIL - Name of map output file
;   NOTRCFITS - No trace structure written to fits
;
; OPTIONAL OUTPUTS:
;   TRCSTRCT - Trace structure
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_mkmap, wfccd_str, mask_id
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_mkmap, struct, mask_id, FLAT=img, NOCHK=nochk, NOFITS=nofits, $
                 VAR=var, OUTFIL=outfil, TRCSTRCT=trcstrct, $
                 NOTRCFITS=notrcfits, YSTRT=ystrt

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'wfccd_mkmap, struct, mask_id, SVOV=, FLAT=, VAR=, OUTFIL=, '
      print, '          TRCSTRCT=, YSTRT= [v1.0]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( YSTRT ) then ystrt = 400L
  
  ; INPUT FLAT, VAR AS NECESSARY
  if not keyword_set(TRCSTRCT) then begin
      flts = where(struct.type EQ 'FLT' AND struct.flg_anly NE 0 AND $
                   struct.mask_id EQ mask_id, nflts)
      if nflts LE 0 then begin
          print, 'wfccd_mkmap: No flat files available!'
          return
      endif
      if not keyword_set( FLAT ) then $
        img = xmrdfits(struct[flts[0]].img_final,/silent)
      if not keyword_set(VAR) then $
        var = xmrdfits(struct[flts[0]].img_final, 1,/silent)
  endif else begin
      sz = size(trcstrct, /type)
      if sz EQ 7 then trcstrct = xmrdfits(trcstrct, /silent)
  endelse

; Create directory if necessary and set out file

  a = findfile('Maps/..', count=count)
  if count EQ 0 then file_mkdir, 'Maps'

  if not keyword_set( OUTFIL ) then begin
      if mask_id LT 10 then outfil = 'Maps/Map_0'+$
        string(mask_id,format='(i1)')+'.fits' $
      else outfil = 'Maps/Map_'+string(mask_id,format='(i2)')+'.fits'
      ; Trace output file
      if mask_id LT 10 then trcfil = 'Maps/Trc_0'+$
        string(mask_id,format='(i1)')+'.fits' $
      else trcfil = 'Maps/Trc_'+string(mask_id,format='(i2)')+'.fits'
  endif
      
; Status
  if not keyword_set( SILENT ) then print, 'wfccd_mkmap: Tracing... '

; Tracing
  if not keyword_set( TRCSTRCT ) then $
    x_traceflat, img, trcstrct, VAR=var, SAW=saw, YSTRT=ystrt

; Release Memory
  if keyword_set(IMG) then delvarx, img
  if keyword_set(VAR) then delvarx, var

; Adjust the trace
  if not keyword_set( NOCHK ) then x_ydtortgui, saw, trcstrct

; Delete saw
  if keyword_set( SAW ) then delvarx, saw

; Output the trace
  if not keyword_set( NOTRCFITS ) then mwrfits, trcstrct, trcfil, /create

; Create the map
  map = x_fit2dtrace(trcstrct, NX=4, NY=4)
  
; Create the Inverse map

  if not keyword_set( SILENT ) then $
    print, 'wfccd_mkmap: Creating inverse map...'
  imap = -1.d*x_rectify(map, map, /DBL)

; Output

  if not keyword_set( NOFITS ) then begin
      ; Update structure
      allimg = where(struct.mask_id EQ mask_id AND struct.flg_anly NE 0)
      struct[allimg].map_fil = outfil
      ;; Output to FITS
      mwrfits, map, outfil, /create
      mwrfits, imap, outfil
  endif

  print, 'wfccd_mkmap: All done!'
  print, '----------------------------------------'

end
