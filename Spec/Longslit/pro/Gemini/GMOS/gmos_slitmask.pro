;+
; NAME:
;   gmos_slitmask
;
; PURPOSE:
;   Determine the positions of the slits on the image, and return tracesets
;
; CALLING SEQUENCE:
;   gmos_slitmask, domefile,outfile 
;
; INPUTS:
;   domefile   - Image for finding the slits, which will be a domeflat
;   outfile    - Output file name with slit mask images and structures
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   tset_slits - 2-element array of trace sets, where the first defines
;                the starting slit positions, and the second one defines
;                the ending slit positions
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   long_proc
;   long_slitmask
;   mwrfits
;   splog
;   trace_crude()
;   xy2traceset
;
; REVISION HISTORY:
;   21-Jul-2008  Written by J. F. Hennawi, Berkeley
;-
;------------------------------------------------------------------------------

;domefile = 'N20080407S0018.fits.gz'
;outfile = 'slitmask.fits'
PRO GMOS_SLITMASK, domefile, maskfile, indxfile, slitfile, CHK = CHK

  ;; put in sensible values????
  IF NOT KEYWORD_SET(MINSLIT) THEN MINSLIT = 0.0
  IF NOT KEYWORD_SET(MAXSLIT) THEN MAXSLIT = 1000.0

  ;;IF NOT KEYWORD_SET(MINSLIT) THEN MINSLIT = 5.0
  ;;IF NOT KEYWORD_SET(MAXSLIT) THEN MAXSLIT = 8.0

  IF NOT KEYWORD_SET(NCOEFF1) THEN NCOEFF = 2 $
  ELSE NCOEFF = NCOEFF1
  ;;IF NOT KEYWORD_SET(SLIT_WIDTH) THEN SLIT_WIDTH = 1.0 ;; slit width in arcsec
  
  ;; Compute DTAX offset
  hdr_mask = xheadfits(maskfile)
  hdr_dome = xheadfits(domefile)
  ;; pixels are 13.5 microns in size
  dtax_dome = float(sxpar(hdr_dome, 'DTAX'))
  dtax_mask = float(sxpar(hdr_mask, 'DTAX'))
  DTAX = round((dtax_dome - dtax_mask)/13.5)
  grating = strcompress(sxpar(hdr_dome, 'GRATING'), /rem)
  long_proc, domefile, flat, /TRANSFORM, bin = bin
  IF strmatch(grating, '*R150*') THEN BEGIN
     readcol, indxfile, slitind, xslit, yslit, skyflag, objno $
              , format = 'L,F,F,L,L'
     plate_scale = bin[0]*0.072
     slitpix = 1.0d/plate_scale ;; 1.0" microslits for old setup
     xslit_sta = xslit - slitpix/2.0d
     xslit_end = xslit + slitpix/2.0d
  ENDIF ELSE BEGIN
     ;; read in slit positions from index file
     readcol, indxfile, slitind, xslit_sta, xslit_end, yslit, skyflag, objno $
              , format = 'L,F,F,L,L,L'
     ;; Mike's slit positions use an iraf convention, subtract off a 1.0
     xslit_sta = xslit_sta-1.0d
     xlist_end = xslit_end-1.0d
  ENDELSE
  ;; GMOS uses 
  pix_dtax = dtax/bin[0]
  xstart = (xslit_sta + pix_dtax) ;;- slitpix/2.0
  xend   = (xslit_end + pix_dtax) ;;+ slitpix/2.0
  nslit = n_elements(slitind)
  ;; Shift slits to align start/end positions with the flat image
  xshift = gmos_alignmask(flat,xstart,xend)
  xstart = xstart + xshift
  xend   = xend + xshift
  ;; should peakthresh be 0.02?
  y1 = 1040L*(bin[1]/2)
  y2 = 2060L*(bin[1]/2)
  ;; look for slits at pre-identified positions
  splog, 'Finding slits.....'
  long_slitmask, domefile, tset_slits = tset_slits0 $
                 , ksize = 1L, /TRANSFORM $
                 , y1 = y1, y2 = y2, ncoeff = ncoeff $
                 , xstart = xstart, xend = xend
  nx = tset_slits0[0].dims[0]
  ny = tset_slits0[0].dims[1]
  dim = size(tset_slits0[0].coeff, /dim)
  ncoeff = dim[0]
  FOR iter = 1L, 2L DO BEGIN
     ;; tweak the slits for two iterations
     splog, 'Tweaking slit edges'
     tset_slits = long_slittweak(flat, tset_slits0, /GMOS, chk = chk)
     splog, 'Checking for bad slits'
     ;; Adjust the broken slits so that they correspond to the 
     ;; what we found in the mask image (with the average slit slope)
     xdiff_sta = abs(xstart - tset_slits[0].coeff[0, *])
     xdiff_end   = abs(xend   - tset_slits[1].coeff[0, *])
     width = tset_slits[1].coeff[0, *]-tset_slits[0].coeff[0, *]
     good_slits = WHERE(width GE MINSLIT AND $
                        width LE MAXSLIT AND $
                        xdiff_sta LE 2.0/bin[0] AND $
                        xdiff_end LE 2.0/bin[0], nkeep $
                        , complement = bad_slits, ncomp = nbad)
     IF nbad GT 0 THEN BEGIN
        splog, 'Iteration #', iter
        splog, 'Found nbad=', strcompress(string(nbad), /rem) $
               , ' bad slits out of nslit=', strcompress(string(nslit), /rem)
        splog, 'Bad slits: '
        splog, '-----------'
        splog, ' slitind     width'
        forprint, string(slitind[bad_slits]), string(width[bad_slits]) $
                  , textout = 2
        splog, 'Bad slits: '
        splog, 'Correcting with mask image positions'
        tset_slits[0].COEFF[0, bad_slits] = xstart[bad_slits]
        tset_slits[1].COEFF[0, bad_slits] = xend[bad_slits]
        FOR kk = 1L, ncoeff-1L DO BEGIN 
           coeff_left  = tset_slits[0].COEFF[kk, good_slits]
           coeff_right = tset_slits[1].COEFF[kk, good_slits]
           djs_iterstat, [coeff_left, coeff_right], median = med_coeff
           tset_slits[0].COEFF[kk, bad_slits] = med_coeff
           tset_slits[1].COEFF[kk, bad_slits] = med_coeff
        ENDFOR
     ENDIF
     tset_slits0 = tset_slits
  ENDFOR
  ;; assign the ystruct values from the index catalog file
  ystruct = replicate(create_struct('SLIT', 0L $
                                    , 'YSLIT_MIN', 0.0, 'YSLIT_MAX', 0.0 $
                                    , 'MASK_XSTA', 0.0, 'MASK_XEND', 0.0 $
                                    , 'MASK_Y', 0.0 $
                                    , 'SKYFLAG', 0L, 'OBJNO', 0L), nslit)
  ystruct.SLIT = long(slitind)
  ystruct.SKYFLAG = long(skyflag)
  ystruct.OBJNO = long(objno)
  ystruct.MASK_XSTA = XSTART
  ystruct.MASK_XEND = XEND
  ystruct.MASK_Y = YSLIT
  
  tslitmask = long_slits2mask(tset_slits)
;; transform slitmask to untransformed frame
  slitmask = gmos_trnimg1to3(tslitmask, tset_slits)
;; deal with hot pixels in slitmask
  ;;badpix = WHERE(slitmask LT 0.0 OR slitmask GT (nslit))
  ;;slitmask[badpix] = 0.0

  mwrfits, slitmask, slitfile, /create
  mwrfits, tset_slits, slitfile
  mwrfits, ystruct, slitfile
  
;atvplot, left1, rows1, psym = 4, color = 2
;atvplot, left2, rows2, psym = 4, color = 2
;atvplot, left3, rows3, psym = 4, color = 2
;atvplot, right1, rows1, psym = 4, color = 1
;atvplot, right2, rows2, psym = 4, color = 1
;atvplot, right3, rows3, psym = 4, color = 1

;trace1 = (left1 + right1)/2.0d
;ind = where(tset_slits1[0].COEFF[1, *] EQ 0.0, COMP = inorm)
;atvplot, left1[*, inorm], rows1[*, inorm], psym = 3, color = 2
;atvplot, right1[*, inorm], rows1[*, inorm], psym = 3, color = 1
;atvplot, left1[*, ind], rows1[*, ind], psym = 3, color = 5
;atvplot, right1[*, ind], rows1, psym = 3, color = 4



END
