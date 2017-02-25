;+ 
; NAME:
; esi_echbadpix   
;     Version 1.1
;
; PURPOSE:
;    Set bad pixels to 0
;
; CALLING SEQUENCE:
;   
;  esi_echbadpix, img
;
; INPUTS:
;   img     - 2D IMG array
;
; RETURNS:
;
; OUTPUTS:
;  Full image with badpixels suppressed
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Only set for 1x1 binning for now
;
; EXAMPLES:
;   esi_echbadpix, img
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   12-Aug-2002 Written by JXP
;   03-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echbadpix, esi, img, NOHOT = NOHOT, HOT_THRESH = HOT_THRESH

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'esi_echbadpix, img [v1.1]'
      return
  endif 

;  Optional Keywords
  cbin = esi[0].cbin
  rbin = esi[0].rbin
  IF NOT KEYWORD_SET(HOT_THRESH) THEN BEGIN
      IF cbin EQ 1 THEN HOT_THRESH = 3.0 $
      ELSE IF cbin EQ 2 THEN HOT_THRESH = 3.0D
  ENDIF

;  Set bad pix to 0
  case esi.cbin of
      1: begin
         case esi.rbin of 
             1: begin
                 ;; HOT SPOT
                 dims = size(img, /dim)
                 nx = dims[0]
                 ny = dims[1]
                 mask = lonarr(nx, ny) + 1L
                 ximg = findgen(nx) # replicate(1.0, ny)
                 yimg = replicate(1.0, nx) # findgen(ny)
                 IF NOT KEYWORD_SET(NOHOT) THEN BEGIN
                     sky_sec = WHERE(ximg LE 60 AND ximg GT 3 AND $
                                     yimg GE 3980 AND yimg LE 4096)
                     djs_iterstat, img[sky_sec], mean = avg_sky
                     region = (ximg LE 205 AND ximg GT 0 AND $
                               yimg GE 3800 AND yimg LE 3850) $
                       OR (ximg LE 205 AND ximg GT 0 AND $
                           yimg GE 3850 AND yimg LE 3890) $
                       OR  (ximg LE 205 AND ximg GT 0 AND $
                            yimg GE 3890 AND yimg LE 3940) $
                       OR  (ximg LE 205 AND ximg GT 0 AND $
                            yimg GE 3940 AND yimg LE 3980)
                     img_bad = WHERE(region AND $
                                     img GE (abs(avg_sky) + $
                                             HOT_THRESH*sqrt(abs(avg_sky))) $
                                     , nhot)
                     IF nhot GT 0 THEN mask[img_bad] = 0
                  ENDIF
                 ;;NEW masking
                 ;; Bad CLMS order 14
                 mask[420:438, 2646:4095] = 0
                 mask[414:444, 2647:2659] = 0
                 ;; Bad CLMS order 11
                 mask[878:890, 3832:3954] = 0
                 mask[891:906, 3814:3865] = 0
                 mask[897, 3814:3975] = 0
                 badpix = WHERE(mask EQ 0, nbad)
                 IF nbad NE 0 THEN  img[badpix] = -1.0d8
                 ;; OLD MASKING
;;                  mask[421:425, 2647:4095] = 0
;;                  mask[432:437, 2647:4095] = 0
;;                  mask[421:437, 2647:2671] = 0
;;                  mask[878:890, 3832:3954] = 0
;;                  mask[891:906, 3814:3865] = 0
;;                  mask[897, 3814:3975] = 0
;;                  mask[418:439, 2647:2659] = 0
;;                  badpix = WHERE(mask EQ 0, nbad)
;;                  IF nbad NE 0 THEN  img[badpix] = -1.0d8
;;               OLDEST MASKING
                 ;; HOT SPOT
;                 img[35:150,3850:3930] = 0.
;                img[151:209,3808:3880] = 0.
;                img[151:190,3881:3960] = 0.
;                ;; Bad CLMS
;                img[421:437,2647:4095] = 0.
;                img[875:905,3827:4095] = 0.
;                img[891:910,3814:3827] = 0.
             end
             else: stop
         endcase
     end
      2: begin
         case esi.rbin of 
             1: begin
                 ;; HOT SPOT
                 dims = size(img, /dim)
                 nx = dims[0]
                 ny = dims[1]
                 mask = lonarr(nx, ny) + 1L
                 ximg = findgen(nx) # replicate(1.0, ny)
                 yimg = replicate(1.0, nx) # findgen(ny)
                 IF NOT KEYWORD_SET(NOHOT) THEN BEGIN
                     sky_sec = WHERE(ximg LE 110 AND $
                                     yimg GE 3530 AND yimg LE 3760)
                     djs_iterstat, img[sky_sec], mean = avg_sky
                     img_bad = WHERE(ximg LE 110 AND $
                                     yimg GE 3760 AND yimg LE 4000 AND $
                                     img GE (avg_sky + $
                                             HOT_THRESH*sqrt(avg_sky)) $
                                     , nhot)
                     IF nhot GT 0 THEN mask[img_bad] = 0
                 ENDIF
                 ;; Bad CLMS 
                 mask[212:213, 2647: *] = 0.
                 mask[217:219, 2647: *] = 0.
                 mask[211:220, 2647:2671] = 0.0
                 badpix = WHERE(mask EQ 0, nbad)
                 IF nbad NE 0 THEN  img[badpix] = -1.0d8
;;                Old hot spot masking
;                 mask[17:75, 3850:3930] = 0.
;                 mask[75:105, 3808:3880] = 0.
;                 mask[50:95, 3881:3960] = 0.
;                 mask[30:56, 3875:3916] = 0.
;;               Old bad clms masking 
;;               mask[211:221, 2647:*] = 0.
             end
             else: stop
         endcase
     end
     4: begin
         case esi.rbin of 
             4: begin
                 mskfil = 'Bias/MSK_'+strtrim(cbin,2)+'x'+strtrim(rbin,2)+'.fits'
                 mskimg = xmrdfits(mskfil, 0, /silent)
                 msk = where(mskimg EQ 0B)
                 img[msk] = 0.
             end
             else: stop
         endcase
     end
     else: stop
 endcase
                 

  return
end
