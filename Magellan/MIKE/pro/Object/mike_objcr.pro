;+ 
; NAME:
; mike_objcr   
;   Version 1.1
;
; PURPOSE:
;    Flags CRs given 2 or more object images.  
;
; CALLING SEQUENCE:
;   
;  mike_objcr, mike, setup, obj_id, side, [iexp]
;
; INPUTS:
;   mike    -  MIKE structure
;   setup   -  Setup ID
;   obj_id  -  Object ID  (e.g. 1L)  (or STD index if /STD set)
;   side    -  Blue (1) and/or Red (2) side  (Default:  [1L,2L])
;   [iexp]  -  Exposure frames (e.g. [0L, 1L])
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
;   /NOFITS - No FITS output
;   RTIO -  Ratio CR must exceed to be flagged (default: 9)
;   /NOCHK -- Do not show the image with CR's
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Algorithm for 3 or more images is not well tested.
;   It would be nice to add a GROW option.
;
; EXAMPLES:
;   mike_objcr, mike, setup, obj_id, side, 
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Aug-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_objcr, mike, setup, obj_id, side, iexp, NOCHK=nochk, RTIO=rtio

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'mike_objcr, mike, setup, obj_id, side, [exp], /NOCHK, RTIO= [v1.1]'
      return
  endif 
  
;  Optional Keywords

  if not keyword_set(RTIO) then rtio = 9.
  if not keyword_set(LTHRESH) then lthresh = 80.
  if not keyword_set(SIDE) then side = [1L,2L]
  if not keyword_set(GROW) then grow = 0. else stop
  
  for ii=0L,n_elements(side)-1 do begin
      NOWRT = 0
      qq = side[ii]
      if qq EQ 1 then print, 'mike_objcr: Blue side' $
      else print, 'mike_objcr: Red side'

      ;;  Find the frames
      exp = where(mike.type EQ 'OBJ' AND mike.flg_anly NE 0 AND $
                   mike.setup EQ setup AND mike.side EQ qq AND $
                   mike.obj_id EQ obj_id, nindx)
      if keyword_set( IEXP ) then exp = exp[iexp]
 
      nindx = n_elements(exp)

      if nindx LT 2 then begin
          print, 'mike_objcr: Less than 2 images!'
          return
      endif

;  Read in first image
      if not keyword_set( SILENT ) then print, 'mike_objcr: Reading images..'

      tmp  = xmrdfits(mike[exp[0]].img_final, /silent)
      sz   = size(tmp, /dimensions)

      flux = fltarr(sz[0],sz[1],nindx)
      ivar = fltarr(sz[0],sz[1],nindx)
      svexp = fltarr(nindx)

      flux[*,*,0] = temporary(tmp) / mike[exp[0]].exp
      ivar[*,*,0]  = xmrdfits(mike[exp[0]].img_final, 1, /silent) * $
            (mike[exp[0]].exp^2)
      svexp[0] = mike[exp[0]].exp

      ;; Read in the rest
      for q=1L,nindx-1 do begin
          flux[*,*,q]  = xmrdfits(mike[exp[q]].img_final, 0, /silent) / $
            mike[exp[q]].exp
          ivar[*,*,q]  = xmrdfits(mike[exp[q]].img_final, 1, /silent) * $
            (mike[exp[q]].exp^2)
          svexp[q] = mike[exp[q]].exp
      endfor

      if nindx EQ 2 then begin
          diff = (flux[*,*,0] - flux[*,*,1])
          hivar = (1./(ivar[*,*,0]>0.)) > (1./(ivar[*,*,1]>0.))
          ;; First image
          snr = diff / sqrt(hivar)
          cr = where(SNR GT rtio AND flux[*,*,0] GT lthresh/svexp[0], na)
          ;; CR
          crmap = ivar[*,*,0] / (svexp[0]^2)
          if na NE 0 then crmap[cr] = -1.
          print, 'mike_objcr:  Identified ', strtrim(na,2), ' cosmic rays in image'+$
            mike[exp[0]].img_final
          ivar[*,*,0] = crmap
          ;; Second Image
          crtot = cr
          cr = where(-1.*SNR GT rtio AND flux[*,*,1] GT lthresh/svexp[1], na)
          crtot = [crtot,cr]
          ;; CR
          crmap = ivar[*,*,1] / (svexp[1]^2)
          if na NE 0 then crmap[cr] = -1.
          print, 'mike_objcr:  Identified ', strtrim(na,2), ' cosmic rays in image'+$
            mike[exp[1]].img_final
          ivar[*,*,1] = crmap
          if not keyword_set( NOCHK ) then begin
              snr[crtot] = -99.
              xatv, snr, /block, min=-99., max=20.
              print, 'mike_objcr:  Continue as you wish.. Or set NOWRT=1 '
              stop
          endif
              
      endif else begin
          ;; Median several images
          medimg = djs_median(flux,3)
          print, 'mike_objcr: Looping!'
          ;; Loop on Images
          crtot = [0L]
          for q=0L,nindx-1 do begin
              tmp = flux[*,*,q]
              rto = tmp/medimg
;              if keyword_set( CHK ) then xatv, rto, /block
              ;; CR
              cr = where(rto GT rtio AND tmp GT lthresh/svexp[q], ncr)
              crtot = [crtot, cr]
              crmap = ivar[*,*,q] / (svexp[q]^2)
              if ncr NE 0 then crmap[cr] = -1.
              print, 'mike_objcr:  Identified ', ncr, ' cosmic rays..'
              ;; Set VAR
              ivar[*,*,q] = crmap
          endfor
          if not keyword_set( NOCHK ) then begin
              print, 'mike_objcr: Check for issues!'
              mdd = median(medimg)
              medimg[crtot] = -99.
              xatv, medimg, /block, min=-3*mdd, max=3*mdd
              print, 'mike_objcr:  Continue as you wish.. Or set NOWRT=1 '
              stop
          endif
      endelse
      
      ;; OUTPUT
      if not keyword_set( NOWRT ) then begin
          if not keyword_set( SILENT ) then print, 'mike_objcr: Writing images...'
          for q=0L,nindx-1 do begin
              head  = xheadfits(mike[exp[q]].img_final, /silent)
              sxaddpar, head, 'CRREJ', 'T'
              mwrfits, flux[*,*,q]*mike[exp[q]].exp, $
                mike[exp[q]].img_final, head, /create
              mwrfits, ivar[*,*,q], mike[exp[q]].img_final 
              ;; COMPRESS
              spawn, 'gzip -f '+mike[exp[q]].img_final
          endfor
      endif
  endfor

  if not keyword_set( SILENT ) then print, 'mike_objcr: All done!'
  return
end
