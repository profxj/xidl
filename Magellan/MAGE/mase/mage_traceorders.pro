;+
; NAME:
;   mage_traceorders
;
; PURPOSE:
;   Determine the positions of the slits on the image, and return tracesets
;
; CALLING SEQUENCE:
;   mage_traceorders
;
; INPUTS:
;   filename   - Image for finding the slits, which would typically
;                be a flat-field image, an arc image, or a sum of those
;   outfile    - Output file name with slit mask positions
;
; OPTIONAL INPUTS:
;   mislit     - Minimum slit width. Default is to return all slits.
;   biasfile   - Bias file to apply to raw images
;   y1         - Starting row number for smashing image to initially
;                identify slits; default to 0.40*NY
;   y2         - Ending row number for smashing image to initially
;                identify slits; default to 0.60*NY
;   nmed       - Width for median-filtering the image first in the wavelength
;                direction (to remove cosmics)
;   ksize      - Half-kernel size for sharpness filter; default to 5 pix
;   peakthresh - Flux threshhold for finding slits; the flux must be
;                at least this fraction of the brightest slit; default to 0.02
;   radius     - Keyword for TRACE_CRUDE; default to same value as KSIZE
;   nave       - Keyword for TRACE_CRUDE; default to 3
;   maxshifte  - Keyword for TRACE_CRUDE; default to 0.1
;   maxshift0  - Keyword for TRACE_CRUDE; default to 1.0
;   func       - Keyword for XY2TRACESET; default to 'legendre'
;   ncoeff     - Keyword for XY2TRACESET; default to 4
;   verbose    - Verbose if set
;
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
;   ;
; REVISION HISTORY:
;   
;-
;------------------------------------------------------------------------------
FUNCTION mage_traceorders, image $
                           , minslit = minslit $
                           , y1 = y1, y2 = y2, nmed = nmed, ksize = ksize $
                           , peakthresh = peakthresh, radius = radius $
                           , nave = nave, maxshifte = maxshifte $
                           , maxshift0 = maxshift0, func = func $
                           , ncoeff = ncoeff, nslitmax = nslitmax $
                           , verbose = verbose, ORDERMASK = ORDERMASK $
                           , CRUDE = CRUDE, ARCHIVE_FILE=ARCHIVE_FILE, CHK = CHK

;mage_proc, '/Users/jhennawi/mage/nite1/mage0071.fits', image, ivar, hdr = hdr
t0 = systime(1)

dims = size(image, /dimens)
nx = dims[0]
ny = dims[1]

;----------
; Set defaults

IF NOT KEYWORD_SET(CRUDE) THEN begin 
   mage_dir=strcompress(getenv('MAGE_DIR'),/rem)
   archive_file = mage_dir + '/Calib/mage_archive_orders.fits'
   print, 'USING ARCHIVE ORDERS'
ENDIF

if (n_elements(y1) EQ 0) then y1 = 1250.0
if (n_elements(y2) EQ 0) then y2 = 1275.0
if (y1 GT y2 OR y1 LT 0 OR y2 GT ny-1) then $
  message, 'Invalid values for Y1,Y2'
if (NOT keyword_set(ksize)) then ksize = 3
if (ksize LT 1 OR ksize GE nx/2-1) then $
  message, 'Invalid kernel size KSIZE'
if (NOT keyword_set(peakthresh)) then peakthresh = 0.01
if (NOT keyword_set(nmed)) then nmed = 5

;  Default values for TRACE_CRUDE
if (NOT keyword_set(radius)) then radius = ksize
if (NOT keyword_set(nave)) then nave = 5
if (NOT keyword_set(maxshifte)) then maxshifte = 5.0
if (NOT keyword_set(maxshift0)) then maxshift0 = 5.0
if (NOT keyword_set(maxerr)) then maxerr = 1.0

;  Default values for parametrizing the traces (XY2TRACESET)
if (NOT keyword_set(func)) then func = 'legendre'
if (NOT keyword_set(ncoeff)) then ncoeff = 4

;   Convolve the input image with the sharpness filter
kern = [-(findgen(ksize)+1), reverse(findgen(ksize)+1)] / (2.*ksize)
image_tmp1 = djs_maskinterp(image, (image LT 10.0), iaxis = 0)
image_tmp2 = djs_maskinterp(image_tmp1, (image_tmp1 LT 10.0), iaxis = 1)
cimg = convol(image_tmp2, kern, /center, /edge_truncate)


norders = 15L
;xleft  = [20, 100, 200, 300, 375, 450, 525, 590, 650, 705, 760, 808, 855, 902 $
;          , 954]
;xright = [100, 200, 300, 375, 450, 525, 590, 650, 705, 760, 808, 855, 902, 962 $
;          , 1019]

;What are these?? Points in between orders..
xleft  = [20, 105, 205, 305, 380, 455, 530, 595, 655, 711, 763, 816, 863, 910 $
          , 957]
xright = [100, 205, 305, 380, 455, 530, 595, 655, 711, 760, 816, 863, 911, 955 $
         ; , 1019]
          , 1002]
         ; , nx-1L]  ; -- KHRR changed this!!!!

;xleft  = [20, 105, 205, 305, 380, 455, 530, 595, 660, 720, 770, 822, 870, 915 $
;          , 964]
;xright = [100, 205, 305, 380, 455, 530, 595, 660, 720, 770, 822, 870, 915, 963 $
;          , nx-1L]

   
;   Smash the center of the sharpness-filtered image to a 1-dimensional vector
fsum = djs_median(cimg[*, y1:y2], 2)

xstart = dblarr(norders)
xend   = dblarr(norders)

FOR iorder = 0L, norders-1L DO BEGIN
   ;; Find candidate starting positions
   fsum_iord = djs_median(fsum[xleft[iorder]:xright[iorder]] $
                          , width = 4, boundary = 'reflect')
   xstart1 = long_find_nminima(-fsum_iord, nfind = 1, minsep = 40L, width = 3)
   ystart = interpolate(fsum_iord, double(xstart1))
   xstart[iorder] = double(xstart1) + double(xleft[iorder])
   xend1  = long_find_nminima(fsum_iord, nfind = 1, minsep = 40L, width = 3)
   
   if iorder EQ 14 then begin 
      if xend1 LT xstart1 then xend1 = xstart1+35.0 
   endif 

   xend[iorder]   = double(xend1) + double(xleft[iorder])
   yend   = interpolate(fsum_iord, double(xend1))
   ;if iorder EQ 14 then stop
ENDFOR

;; Make sure that xstart is always before xend
FOR iorder = 0L, norders-1L DO BEGIN
   IF xend[iorder] LT xstart[iorder] THEN BEGIN
      xs1 = xstart[iorder]
      xe1 = xend[iorder]
      xend[iorder] = xs1
      xstart[iorder] = xe1 
  ENDIF
ENDFOR


   ;; Find candidate ending positions
;  peakcut = peakthresh*interpolate(fsum, xmax)
;   peakcut = peakthresh*interpolate(-fsum_iord, xmax)
;   xstart = long_find_nminima(smooth(-fsum, 3) < peakcut, nfind = 1 $
;   , minsep = 100L, width = 5)
;   nstart = n_elements(xstart)
;   xend = long_find_nminima(smooth(fsum, 3) < peakcut, nfind = 20 $
;                         , minsep = 100L, width = 5)
;xend = double(xend)
;yend = interpolate(fsum, xend)
;nend = n_elements(xend)

;peakval = 1
;if (nstart GT 0) then peakval = peakval > max(fsum[xstart])
;if (nend GT 0) then peakval = peakval > max((-fsum[xend]))
;if (nstart GT 0) then begin
;    igood = where(fsum[xstart] GT peakthresh*peakval, ngood)
;    if (ngood GT 0) then xstart = xstart[igood] $
;    else xstart = [0]
;endif else xstart = [0]
;if (nend GT 0) then begin
;    igood = where(fsum[xend] LT -peakthresh*peakval, ngood)
;    if (ngood GT 0) then xend = xend[igood] $
;    else xend = [nx-1]
;endif else xend = [nx-1]
; Make sure we start with a "starting" position, and end with an "ending"
;if (xstart[0] GT xend[0]) then $
;  xstart = [0, xstart]
;if (xend[n_elements(xend)-1] LT xstart[n_elements(xstart)-1]) then $
;  xend = [xend, nx-1]
; Find pairs of starting+ending slit positions, discarding pairs
; of start or end positions that should not exist.
; For book-keeping, set ALLQ=+1 for the start of a slit, =-1 for the end.
;allx = float([xstart, xend])
;allq = [replicate(1, n_elements(xstart)) $
;        , replicate(-1, n_elements(xend))]
;isort = sort(allx)
;allx = allx[isort]
;allq = allq[isort]

;   ; Keep the last starting position if there are several in a row,
;   ; and the first ending position if there are several in a row.
;igood = where((allq EQ 1 AND allq NE shift(allq, -1)) $
;              OR (allq EQ -1 AND allq NE shift(allq, 1)), ngood)

;allx = allx[igood]
;xstart = allx[lindgen(ngood/2)*2]
;xend = allx[lindgen(ngood/2)*2+1]
;nslit = n_elements(xstart)
;stop
;IF KEYWORD_SET(NSLITMAX) THEN BEGIN
;    xstart = xstart[0:NSLITMAX-1L]
;    xend   = xend[0:NSLITMAX-1L]
;    nslit = NSLITMAX
;ENDIF

; Trim to slits with width > 5 pixels
ycen = (y1+y2)/2.0
minval = -djsig(cimg)     ; Clip negative values on image when tracing


IF KEYWORD_SET(ARCHIVE_FILE) THEN BEGIN 
;; Use archived slit traces as crutch for tracing
    tset_guess = mrdfits(archive_file[0], 1)  
    traceset2xy, tset_guess[0], rows, left_edge
    traceset2xy, tset_guess[1], rows, right_edge
    left_crutch = 0*left_edge
    right_crutch = 0*right_edge
    FOR islit = 0L, norders-1L DO BEGIN
        left_crutch[*, islit] = left_edge[*, islit] - $
          (interpol(left_edge[*, islit], rows[*, islit], ycen))  + xstart[islit]
        right_crutch[*, islit] = right_edge[*, islit] - $
          (interpol(right_edge[*, islit], rows[*, islit], ycen))  + xend[islit]
    ENDFOR
ENDIF ELSE BEGIN
    ;; use trace crude to determine the crutch for slit tracing
    radius = 3
;   Trace the starting+end positions on the CCD.
    xpos1 = trace_crude(cimg > minval, xstart = xstart $
                        , ystart = replicate(ycen, n_elements(xstart)) $
                        , radius = radius, maxerr = maxerr  $
                        , nave = nave, maxshifte = maxshifte $
                        , maxshift0 = maxshift0, xerr = xerr1)
    xpos2 = trace_crude(-cimg > minval, xstart = xend $
                        , ystart = replicate(ycen, n_elements(xstart)) $
                        , radius = radius, maxerr = maxerr $
                        , nave = nave, maxshifte = maxshifte $
                        , maxshift0 = maxshift0, xerr = xerr2)
;  Fit these traces by a polynomial function
    xy2traceset, rebin(findgen(ny), ny, norders), xpos1, sset1 $
                 , invvar = 1./xerr1^2, func = func, ncoeff = ncoeff
    xy2traceset, rebin(findgen(ny), ny, norders), xpos2, sset2 $
                 , invvar = 1./xerr2^2, func = func, ncoeff = ncoeff
    tset_guess = [sset1, sset2]
    tset_guess = struct_addtags(tset_guess $
                                , replicate(create_struct('DIMS', dims) $
                                            , size(tset_guess, /dimens)))
    traceset2xy, tset_guess[0], rows, left_crutch
    traceset2xy, tset_guess[1], rows, right_crutch
 ENDELSE
;; Use a tighter kernel for edge tracing
;ksize = 3
;kern = [-(findgen(ksize)+1), reverse(findgen(ksize)+1)] / (2.*ksize)
;cimg = convol(image_tmp2, kern, /center, /edge_truncate)
niter = 12L
xfit1 = left_crutch
FOR i = 1L, niter DO BEGIN
    IF i LT niter/3 THEN rad = 3.0*ksize $
    ELSE IF (i GE niter/3) AND (i LT 2*niter/3) THEN rad = ksize $
    ELSE rad = ksize/2.0d
    xpos1 = trace_fweight(cimg > 0.0, xfit1, rows, radius = rad)
    xy2traceset, rows, xpos1, left_set1, ncoeff = 5, yfit = xfit1 $
                 , maxdev = 1.0D, /silent
ENDFOR
xfit2 = xfit1
FOR i = 1L, niter DO BEGIN
    xpos2 = trace_gweight(cimg > 0.0, xfit2, rows $
                          , sigma = KSIZE/3.0D/2.3548D)
    xy2traceset, rows, xpos2, left_set, ncoeff = 5, yfit = xfit2 $
                 , maxdev = 1.0D, /silent
ENDFOR
xfit1 = right_crutch
FOR i = 1L, niter DO BEGIN
    IF i LT niter/3 THEN rad = 3.0*ksize $
    ELSE IF (i GE niter/3) AND (i LT 2*niter/3) THEN rad = ksize $
    ELSE rad = ksize/2.0d
    xpos1 = trace_fweight(-cimg > 0.0, xfit1, rows, radius = rad)
    xy2traceset, rows, xpos1, right_set1, ncoeff = 5, yfit = xfit1 $
                 , maxdev = 1.0D, /silent
ENDFOR
xfit2 = xfit1
FOR i = 1L, niter DO BEGIN
    xpos2 = trace_gweight(-cimg > 0.0, xfit2, rows $
                          , sigma = KSIZE/3.0D/2.3548D)
    xy2traceset, rows, xpos2, right_set, ncoeff = 5, yfit = xfit2 $
                 , maxdev = 1.0D, /silent
ENDFOR
tset_slits = [left_set, right_set]
tset_slits = struct_addtags(tset_slits $
                            , replicate(create_struct('DIMS', dims) $
                                        , size(tset_slits, /dimens)))
;  Generate a mask of the order numbers
ordermask = long_slits2mask(tset_slits)
ordermask[WHERE(ordermask GT 0)] = -ordermask[WHERE(ordermask GT 0)] + 21L

splog, 'Elapsed time = ', systime(1)-t0, ' sec'

tmp = {ordr_str, $
       order: 0, $
       flg_anly: 0, $
       xcen: 0.0, $
       ycen: 1024.0, $
       ymin: 0L, $
       ymax: 2048L, $
       lcen: 0L, $
       lhedg: fltarr(2048), $
       rhedg: fltarr(2048), $
       arc_m: fltarr(2048), $
       profile0: fltarr(251), $
       profile1: fltarr(251), $
       lhc: fltarr(6), $
       rhc: fltarr(6) $
      }

ordr_str = replicate(tmp, 15)

traceset2xy, left_set, xl, yl
traceset2xy, right_set, xr, yr

for i=0, 14 do begin
   ordr_str[i].order = 20-i
   ordr_str[i].lhedg = yl[*,i]
   ordr_str[i].rhedg = yr[*,i] 
   ordr_str[i].lhc[0:4] = left_set.coeff[*,i]
   ordr_str[i].rhc[0:4] = right_set.coeff[*,i]
endfor

mwrfits, ordr_str, "OStr_mage.fits", /create

traceset2xy, tset_slits[0], rows, left_edge
traceset2xy, tset_slits[1], rows, right_edge
slit_edg = fltarr(ny, norders, 2)
slit_edg[*, *, 0] = left_edge
slit_edg[*, *, 1] = right_edge

;; CHK (Plot slit edges)
if keyword_set( CHK ) then begin
   tmp = image
   for i = 0L, norders-1L DO BEGIN
      for j = 0, 1 do begin
         rnd_trc2 = 0 > round(slit_edg[*, i, j]) < nx-1
         trc_msk = rnd_trc2 + lindgen(ny)*nx
         tmp[trc_msk] = -10000
         if keyword_set(IFU) then begin
            if j EQ 0 then begin
               rnd_trc2 = round(slit_edg[*, i, 0]+ $
                                (slit_edg[*, i, 1]-slit_edg[*, i, 0])*0.53)
               trc_msk = rnd_trc2 + lindgen(ny)*nx
               tmp[trc_msk] = -10000
               rnd_trc2 = round(slit_edg[*, i, 0]+ $
                                (slit_edg[*, i, 1]-slit_edg[*, i, 0])*0.57)
               trc_msk = rnd_trc2 + lindgen(ny)*nx
               tmp[trc_msk] = -10000
            endif
         endif
      endfor
   endfor
   xatv, tmp, /block, min = -10.0, max = 4.0d4, /histeq
endif


return, tset_slits
end
;------------------------------------------------------------------------------
