;+
; NAME:
;   gmos_trnimg1to3
;
; PURPOSE:
;
;   Transforms a single transformed image into 3-ccd images which are
;   displayed in a gmos_mosaic mode. 
;
; CALLING SEQUENCE:
;  
;
; INPUTS:
;  transimg -- Transformed image. 
;
; OPTIONAL INPUTS:
; /NOFLEX  -- Do not apply flexure correction [necessary if your setup
;             has not been calibrated.  Contact JH or JXP for help if
;             this is the case.]
;  HAND_FWHM -- Set the FWHM of the object profile to this value (in
;               pixels)
; /NOHELIO -- Do not correct to heliocentric velocities
; /NOZAP   -- Do not flag CRs
;
; OUTPUTS:
;   outimg = mosaiced 3-ccd image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   
; PROCEDURES CALLED:
;   
; REVISION HISTORY:
;   28-Jul-2008  Written by Joe Hennawi and Mike Gladders
;-  
;-----------------------------------------------------------------------------

FUNCTION gmos_trnimg1to3, transimg, tset_slits

  
  nx = tset_slits[0].dims[0]
  ny = tset_slits[0].dims[1]
  idim = size(tset_slits[0].COEFF, /dim)
  nslit = idim[1]
  yarr = findgen(ny) ## replicate(1.0, nx)
  xarr = findgen(nx) # replicate(1.0, ny)
  traceset2xy, tset_slits[0], yy1, xx1
  traceset2xy, tset_slits[1], yy2, xx2
  ;; slitmask in tranformed frame
  ;;tslitmask = long_slits2mask(tset_slits)
  utransimg = 0*transimg

  IF      (ny-2L*9L)  MOD  512 EQ 0 THEN specbin = 4 $
  ELSE IF (ny-2L*18L) MOD  512 EQ 0 THEN specbin = 2 $
  ELSE IF (ny-2L*36L) MOD  512 EQ 0 THEN specbin = 1 $
  ELSE message, 'Transformation not supported for your binning'
  
  spatbin = long(4L*1152/nx)
  
  ;; Binning is spec x spatial (i.e. X x Y in Gemini convention)
  IF specbin EQ 2 AND spatbin EQ 2 THEN BTAG = '2x2' $
  ELSE IF specbin EQ 1 AND spatbin EQ 1 THEN BTAG = '1x1' $
  ELSE IF specbin EQ 2 AND spatbin EQ 1 THEN BTAG = '2x1' $
  ELSE message, 'Your binning is not supported'
  
  ;; This is like nearest grid-point interpolation. We create a new
  ;; image in untransformed frame with same pixel values as
  ;; transformed image with one-to-one mapping between pixels
  FOR islit = 0L, nslit-1L DO BEGIN
     xx1_islit = replicate(1.0, nx) # xx1[*, islit]
     xx2_islit = replicate(1.0, nx) # xx2[*, islit]
     not_edge = xarr GT round(xx1_islit) AND xarr LT round(xx2_islit)
     ;; Transform the non-edge pixels on the slit
     ipix = WHERE(not_edge)
     gmos_trnxy1to3, xarr[ipix], yarr[ipix], x3 = x3, y3 = y3 $
                     , BTAG = BTAG 
     utransimg[round(x3), round(y3)] = transimg[ipix]
     ;; Transform the edge pixels. We treat the edge pixels
     ;; differently since otherwise we would be rounding off twice and
     ;; slightly screwing up the edges. 
     gmos_trnxy1to3, xx1[*, islit], yy1[*, islit], x3 = xx1_3, y3 = yy1_3 $
                     , BTAG = BTAG 
     utransimg[round(xx1_3), round(yy1_3)] = $
        transimg[round(xx1[*, islit]), round(yy1[*, islit])]
     gmos_trnxy1to3, xx2[*, islit], yy2[*, islit], x3 = xx2_3, y3 = yy2_3 $
                     , BTAG = BTAG 
     utransimg[round(xx2_3), round(yy2_3)] = $
        transimg[round(xx2[*, islit]), round(yy2[*, islit])]
  ENDFOR
  ;; This was the old way. Notice that the edges of the trace get
  ;; smoothed since they are rounded both in the slitmask image and in
  ;; the transformation. 
  ;;FOR islit = 0L, nslit-1L DO BEGIN
  ;;   ipix = WHERE(tslitmask EQ (islit+1L))
  ;;   gmos_trnxy1to3, xarr[ipix], yarr[ipix], x3 = x3, y3 = y3 $
  ;;                   , BTAG=BTAG 
  ;;   utransimg[round(x3), round(y3)] = transimg[ipix]
  ;;ENDFOR
  
  RETURN, utransimg
END
