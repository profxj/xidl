;+ 
; NAME:
; wfccd_objcr   
;   Version 1.0
;
; PURPOSE:
;    Flags CRs given 2 or more object images
;
; CALLING SEQUENCE:
;   
;   wfccd_objcr, struct, mask_id, [expsr]
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
;   wfccd_objcr, nght1_strct, mask_id
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   14-May-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_objcr, wfccd, mask_id, expsr

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'wfccd_objcr, wfccd, mask_id, [expsr] [v1.0]'
      return
  endif 
  
;  Optional Keywords

  if not keyword_set(RTIO) then rtio = 4.
  if not keyword_set(LTHRESH) then lthresh = 200.
  
;  Find the frames

  allexp = where(wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0 AND $
                 wfccd.mask_id EQ mask_id, nflts)
  if not keyword_set(EXPSR) then exp = allexp else exp = allexp[expsr]

  nexp = n_elements(exp)
  if nexp LT 2 then begin
      print, 'wfccd_objcr: Less than 2 images!'
      return
  endif

;  Read in first image
  if not keyword_set( SILENT ) then print, 'wfccd_objcr: Reading images..'

  tmp  = xmrdfits(wfccd[exp[0]].img_final, /silent)
  tmp2 = xmrdfits(wfccd[exp[0]].img_final, 2, /silent)
  sz   = size(tmp, /dimensions)
  sz2  = size(tmp, /type)

  flux = fltarr(sz[0],sz[1],nexp)
 var = fltarr(sz[0],sz[1],nexp)
  wave = make_array(sz[0],sz[1],nexp,type=sz2)

  flux[*,*,0] = temporary(tmp)
  var[*,*,0]  = xmrdfits(wfccd[exp[0]].img_final, 1, /silent)
  wave[*,*,0]  = temporary(tmp2)

; Read in the rest
  for q=1L,nexp-1 do begin
      flux[*,*,q]  = xmrdfits(wfccd[exp[q]].img_final, 0, /silent)
      var[*,*,q]  = xmrdfits(wfccd[exp[q]].img_final, 1, /silent)
      wave[*,*,q]  = xmrdfits(wfccd[exp[q]].img_final, 2, /silent)
  endfor

  if nexp EQ 2 then begin
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
      print, 'wfccd_objcr: Looping!'
      ; Loop on Images
      for q=0L,nexp-1 do begin
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

  for q=0L,nexp-1 do begin
      mwrfits, flux[*,*,q], wfccd[exp[q]].img_final, /create
      mwrfits, var[*,*,q], wfccd[exp[q]].img_final
      mwrfits, wave[*,*,q], wfccd[exp[q]].img_final
      ;; COMPRESS
      spawn, 'gzip -f '+wfccd[exp[q]].img_final
  endfor

  if not keyword_set( SILENT ) then print, 'wfccd_objcr: All done!'
  return
end
