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
;   20-Sep-2004 Major rewrite for 2 images
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echobjcr, esi, obj_id, RTIO=rtio, LTHRESH=lthresh, CHK=chk

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echobjcr, esi, obj_id, [v1.0]'
      return
  endif 
  
;  Optional Keywords

  if not keyword_set(RTIO) then rtio = 9.
  if not keyword_set(LTHRESH) then lthresh = 80.
  
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
  svexp = fltarr(nindx)

  flux[*,*,0] = temporary(tmp) / esi[indx[0]].exp
  var[*,*,0]  = xmrdfits(esi[indx[0]].img_final, 1, /silent) / $
    (esi[indx[0]].exp^2)
  svexp[0] = esi[indx[0]].exp

  ;; Read in the rest
  for q=1L,nindx-1 do begin
      flux[*,*,q]  = xmrdfits(esi[indx[q]].img_final, 0, /silent) / $
        esi[indx[q]].exp
      var[*,*,q]  = xmrdfits(esi[indx[q]].img_final, 1, /silent) / $
        (esi[indx[q]].exp^2)
      svexp[q] = esi[indx[q]].exp
  endfor

  if nindx EQ 2 then begin
      diff = (flux[*,*,0] - flux[*,*,1])
      hivar = (var[*,*,0]>0.) > (var[*,*,1]>0.)
      ;; First image
      snr = diff / sqrt(hivar)
      cr = where(SNR GT rtio AND flux[*,*,0] GT lthresh/svexp[0], na)
      ;; CR
      crmap = var[*,*,0] * (svexp[0]^2)
      if na NE 0 then crmap[cr] = -1.
      print, 'esi_objcr:  Identified ', strtrim(na,2), ' cosmic rays in image'+$
        esi[indx[0]].img_final
      var[*,*,0] = crmap
      ;; Second Image
      crtot = cr
      cr = where(-1.*SNR GT rtio AND flux[*,*,1] GT lthresh/svexp[1], na)
      crtot = [crtot,cr]
      ;; CR
      crmap = var[*,*,1] * (svexp[1]^2)
      if na NE 0 then crmap[cr] = -1.
      print, 'esi_objcr:  Identified ', strtrim(na,2), ' cosmic rays in image'+$
        esi[indx[1]].img_final
      var[*,*,1] = crmap
      if keyword_set( CHK ) then begin
          snr[crtot] = -99.
          xatv, snr, /block, min=-99., max=20.
          print, 'esi_objcr:  Continue as you wish.. Or set NOWRT=1 '
          stop
      endif
      
  endif else begin
      ;; Median several images
      medimg = djs_median(flux,3)
      print, 'esi_objcr: Looping!'
      ;; Loop on Images
      crtot = [0L]
      for q=0L,nindx-1 do begin
          tmp = flux[*,*,q]
          rto = tmp/medimg
;              if keyword_set( CHK ) then xatv, rto, /block
          ;; CR
          cr = where(rto GT rtio AND tmp GT lthresh/svexp[q], ncr)
          crtot = [crtot, cr]
          crmap = var[*,*,q] * (svexp[q]^2)
          if ncr NE 0 then crmap[cr] = -1.
          print, 'esi_objcr:  Identified ', ncr, ' cosmic rays..'
          ;; Set VAR
          var[*,*,q] = crmap
      endfor
      if keyword_set( CHK ) then begin
          print, 'esi_objcr: Check for issues!'
          mdd = median(medimg)
          medimg[crtot] = -99.
          xatv, medimg, /block, min=-3*mdd, max=3*mdd
          print, 'esi_objcr:  Continue as you wish.. Or set NOWRT=1 '
          stop
      endif
  endelse
      
  ;; OUTPUT
  if not keyword_set( NOWRT ) then begin
      if not keyword_set( SILENT ) then print, 'esi_objcr: Writing images...'
      for q=0L,nindx-1 do begin
          mwrfits, flux[*,*,q]*esi[indx[q]].exp, $
            esi[indx[q]].img_final, /create
          mwrfits, var[*,*,q], esi[indx[q]].img_final 
          ;; COMPRESS
          spawn, 'gzip -f '+esi[indx[q]].img_final
      endfor
  endif

  if not keyword_set( SILENT ) then print, 'esi_echobjcr: All done!'
  return
end
