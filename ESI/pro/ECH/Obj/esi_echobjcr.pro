;+ 
; NAME:
; esi_echobjcr   
;   Version 1.0
;
; PURPOSE:
;    Flags CRs given 2 or more object images
;
; CALLING SEQUENCE:
;   
;   esi_echobjcr, esi, obj_id
;
; INPUTS:
;   struct -- esi_echstrct defining the images of interest
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
;   esi_echobjcr, esi, 0L
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   06-Sep-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echobjcr, esi, obj_id, RTIO=rtio, LTHRESH=lthresh

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echobjcr, esi, obj_id, [v1.0]'
      return
  endif 
  
;  Optional Keywords

  if not keyword_set(RTIO) then rtio = 4.
  if not keyword_set(LTHRESH) then lthresh = 100.
  
;  Find the frames

  indx = where(esi.flg_anly NE 0 AND esi.mode EQ 2 AND $
               esi.obj_id EQ obj_id AND strtrim(esi.type,2) EQ 'OBJ', nindx)

  if nindx LT 2 then begin
      print, 'esi_echobjcr: Less than 2 images!'
      return
  endif

;  Read in first image
  if not keyword_set( SILENT ) then print, 'esi_echobjcr: Reading images..'

  tmp  = xmrdfits(esi[indx[0]].img_final, /silent)
  sz   = size(tmp, /dimensions)

  flux = fltarr(sz[0],sz[1],nindx)
  var = fltarr(sz[0],sz[1],nindx)

  flux[*,*,0] = temporary(tmp)
  var[*,*,0]  = xmrdfits(esi[indx[0]].img_final, 1, /silent)

; Read in the rest
  for q=1L,nindx-1 do begin
      flux[*,*,q]  = xmrdfits(esi[indx[q]].img_final, 0, /silent)* $
        (esi[indx[0]].exp/esi[indx[q]].exp)  ; Scale by exposure
      var[*,*,q]  = xmrdfits(esi[indx[q]].img_final, 1, /silent)
  endfor

  if nindx EQ 2 then begin
      tmp0 = flux[*,*,0]
      tmp1 = flux[*,*,1] 
      ;; First image
      a = where(tmp1 NE 0)
      rto = fltarr(sz[0],sz[1])
      rto[a] = tmp0[a]/tmp1[a]
      ;; CR
      crmap = var[*,*,0]
      cr = where(rto GT rtio AND tmp0 GT lthresh, ncr)
      if ncr NE 0 then crmap[cr] = -1.
      var[*,*,0] = crmap
      ;; Second Image
      a = where(tmp0 NE 0)
      rto = fltarr(sz[0],sz[1])
      rto[a] = tmp1[a]/tmp0[a]
      ;; CR
      crmap = var[*,*,1]
      cr = where(rto GT rtio AND tmp1 GT lthresh, ncr)
      if ncr NE 0 then crmap[cr] = -1.
      var[*,*,1] = crmap
  endif else begin
      ; Median
      medimg = djs_median(flux,3)
      ; Loop on Images
      print, 'esi_echobjcr: Looping!'
      for q=0L,nindx-1 do begin
          tmp = flux[*,*,q]
          rto = tmp/medimg
          ;; CR
          cr = where(rto GT rtio AND tmp GT lthresh, ncr)
          crmap = var[*,*,q]
          if ncr NE 0 then crmap[cr] = -1.
          ;; Set VAR
          var[*,*,q] = crmap
      endfor
  endelse

  ; OUTPUT

  for q=0L,nindx-1 do begin
      mwrfits, flux[*,*,q]*(esi[indx[q]].exp/esi[indx[0]].exp), $
        esi[indx[q]].img_final, /create                 ; Put exposure back
      mwrfits, var[*,*,q], esi[indx[q]].img_final
      spawn, 'gzip -f '+esi[indx[q]].img_final
  endfor

  if not keyword_set( SILENT ) then print, 'esi_echobjcr: All done!'
  return
end
