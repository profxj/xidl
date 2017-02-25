;+
; NAME:
;   long_slitmask
;
; PURPOSE:
;   Determine the positions of the slits on the image, and return tracesets
;
; CALLING SEQUENCE:
;   long_slitmask, filename, outfile, $
;    [minslit=,biasfile=, y1=, y2=, nmed=, ksize=, peakthresh=, $
;    radius=, nave=, maxshifte=, maxshift0=, func=, ncoeff=, tset_slits= ]
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
;   ncoeff     - Keyword for XY2TRACESET; default to 3
;   verbose    - Verbose if set
;   /SNGL      - Forces the code to assume one (very long) long slit
;   EDIT_SEDGE_FIL 
;              - Name of a text file describing how to trim specific
;                slits.  The file should have 5 columns:
;                1. x value near the middle of the 'problem' slit
;                2. '1' for the left edge of the slit, or '2' for the
;                right edge
;                3. number of pixels to trim from that particular edge
;                4. pixel number (in y direction) at which to start
;                the adjustment
;                5. pixel number (in y direction) at which to
;                end the adjustment
;   ADD_SLITS - Name  of a text file describing new slits to force
;               into the mask.  The file should have 2 columns:
;               1. x value at the start of the new slit
;               2. x value at the end of the new slit
;   COMBINE_SLITS - Used to combine two adjacent slits into one.  This
;                   is a vector, which should be set to the lower slit
;                   ID of the two you want to combine.  (I.e., if you
;                   want to combine slits 9 and 10, set COMBINE_SLITS=[9])
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
;   Need to deal with cases where slits appear to overlap,
;   or if a slit position goes crazy.  At a minimum, we should
;   toss such cases???
;   Allow several input files, like a dome and an arc, and then
;   align those images and add them???
;
;   I think peakthresh should be much lower like 0.02, especially
;   if slitwidths are widely varying  SMB, 3-25-05
;
; PROCEDURES CALLED:
;   long_proc
;   long_slits2mask()
;   mwrfits
;   splog
;   trace_crude()
;   xy2traceset
;
; REVISION HISTORY:
;   10-Mar-2005  Written by D. Schlegel, LBL
;   13-Oct-2005  Edited by J. Simon, Caltech to assume that solo starting 
;                and ending slit positions are real, and need a complementary
;                 ending or starting position generated for them
;-
;------------------------------------------------------------------------------
; Add the vector XPOS_ADD (and XERR_ADD) to the array XPOS (XERR)
; as the i-th row of that array.
pro lucifer_slitmask_add1, i, xpos, xerr, xpos_add, xerr_add

   ; NSLIT is the number of rows before adding the one additional row.
   nslit = n_elements(xpos) / n_elements(xpos_add)

   if (i EQ 0) then begin
      xpos = [[xpos_add], [xpos]]
      xerr = [[xerr_add], [xerr]]
   endif else if (i EQ nslit) then begin
      xpos = [[xpos], [xpos_add]]
      xerr = [[xerr], [xerr_add]]
   endif else begin
      xpos = [[xpos[*,0:i-1]], [xpos_add], [xpos[*,i:nslit-1]]]
      xerr = [[xerr[*,0:i-1]], [xerr_add], [xerr[*,i:nslit-1]]]
   endelse

   return
end
;------------------------------------------------------------------------------
;pro lucifer_slitmask, filename, outfile, biasfile = biasfile $
;                   , xstart = xstart1, xend = xend1 $
;                   , minslit = minslit, y1 = y1_in, y2 = y2_in, nmed = nmed $
;                   , ksize = ksize, peakthresh = peakthresh $
;                   , radius = radius, nave = nave, maxshifte = maxshifte $
;                   , maxshift0 = maxshift0, func = func, ncoeff = ncoeff $
;                   , verbose = verbose, tset_slits = tset_slits $
;                   , nfind = nfind1, GMOSLONG = GMOSLONG $
;                   , minsep = minsep1, CCDONLY = CCDONLY, CIMG = CIMG $
;                   , TRANSFORM = TRANSFORM, TRIM_SEDGE=trim_sedge $
;                   ;, REMOVE_OVERLAP=remove_overlap, $
;                   , EDIT_SEDGE_FIL=edit_sedge_fil, ADD_SLITS=add_slits $;
;                   , SPLIT_SLITS = split_slits, COMBINE_SLITS=combine_slits

t0 = systime(1)

;----------
; Read raw image
outfile = '/Users/joe/luci_data/slits-luci.20100312.0098.fits'
filename = '/Users/joe/luci_data/calib/luci.20100312.0098.fits'
obj = mrdfits(filename, 0, hdr, /fscale)
image = reverse(transpose(obj), 2)
gain = float(sxpar(hdr, 'ELECGAIN'))
rdnoise = float(sxpar(hdr, 'ENOISE'))
image = gain*image
ivar = 1.0/(abs(image - sqrt(2.0)*rdnoise) +rdnoise^2)
peakthresh = 0.1D
minslit = 20
;; GMOS-NX
;long_proc, filename, image, ivar $
;           , hdr = hdr, biasfile = biasfile, CCDONLY = CCDONLY $
;           , verbose = verbose, TRANSFORM = TRANSFORM
    
dims = size(image, /dimens)
nx = dims[0]
ny = dims[1]

telescope  = strcompress(sxpar(hdr[*, 0], 'TELESCOP'), /rem)
instrument = strcompress(sxpar(hdr[*, 0], 'INSTRUME'), /rem)
detector   =  strcompress(sxpar(hdr[*, 0], 'DETECTOR'), /rem)
fpa        = strcompress(sxpar(hdr[*, 0], 'FPA'), /rem)

;; Set default parameters
if (stregex(sxpar(hdr[*,0],'INSTRUME'),'.*kast.*',/boolean,/fold_case) eq 1) or $
   (stregex(sxpar(hdr[*,0],'VERSION'),'kast*', /bool, /fold_case) EQ 1) then begin
;      tset_slits = kast_slitset(nx, ny) $  ;; Too much Kast red fringing..
    peakthresh = 0.1
    nfind1 = 1
endif
 
IF strmatch(instrument, 'MMT Blue*') THEN $
  tset_slits = mmt_slitset(nx, ny, CHANNEL = 'BLUE') $
ELSE IF strmatch(instrument, 'mmtredchan') THEN $
   tset_slits = mmt_slitset(nx, ny, CHANNEL = 'NEWRED') $
ELSE IF strmatch(instrument, '*TWIN*') THEN $
   tset_slits = caha_slitset(nx, ny, hdr) $
ELSE IF strmatch(telescope, '*Gemini-North*') AND KEYWORD_SET(GMOSLONG) $
THEN tset_slits = gmos_slitset(nx, ny) $
ELSE IF sxpar(hdr[*,0], 'XGMOS')  then begin 
   cbin = strtrim(strsplit(sxpar(hdr[*,0],'CCDSUM'),/extract),2)
   bin = reverse(long(cbin))
   tset_slits = xgmos_slitset(nx, ny, bin[0]) 
ENDIF ELSE BEGIN
;; Putting back this kludge..
;  OR strmatch(fpa, 'DBSP_BLUE') OR $
;   strmatch(fpa, 'DBSP_RED') THEN $

    
   ;;----------
   ;; Set defaults

    y1_default = floor(0.40*ny)
    y2_default = ceil(0.50*ny)

   if (NOT keyword_set(ksize)) then ksize = 3
   if (ksize LT 1 OR ksize GE nx/2-1) then $
    message, 'Invalid kernel size KSIZE'
   if (NOT keyword_set(peakthresh)) then peakthresh = 0.02
   if (NOT keyword_set(nmed)) then nmed = 5
   IF NOT KEYWORD_SET(nfind1) THEN nfind = 100L $
   ELSE nfind = nfind1
   ; Default values for TRACE_CRUDE
   if (NOT keyword_set(radius)) then radius = ksize
   if (NOT keyword_set(nave)) then nave = 3
   if (NOT keyword_set(maxshifte)) then maxshifte = 0.1
   if (NOT keyword_set(maxshift0)) then maxshift0 = 1.0
   if keyword_set(minsep1) then minsep = minsep1 $
   ELSE minsep = 10L
   
   ; Default values for parametrizing the traces (XY2TRACESET)
   if (NOT keyword_set(func)) then func = 'legendre'
   if (NOT keyword_set(ncoeff)) then ncoeff = 3

   ; Median-filter the image in the wavelength direction, just to
   ; remove cosmics
   if (nmed GT 1) then $
    for i=0L, nx-1L do $
     image[i,*] = djs_median(reform(image[i,*],ny), $
     width=nmed, boundary='reflect')

   ;; Convolve the input image with the sharpness filter
   kern = [-(findgen(ksize)+1), reverse(findgen(ksize)+1)] / (2.*ksize)
   mask = double(ivar GT 0.0)
   cimg1 = convol(image*mask, kern, /center, /edge_truncate)
   cmask = convol(mask, kern, /center, /edge_truncate)
   cimg = cimg1*(cmask EQ 0)
   ;; Unless y1 and y2 are set, crudely determine y values for which there is 
   ;; light in the flat
   IF n_elements(y1_in) GT 0 AND n_elements(y2_in) GT 0 THEN BEGIN
       y1 = y1_in
       y2 = y2_in
       if (y1 GT y2 OR y1 LT 0 OR y2 GT ny-1) then $
         message, 'Invalid values for Y1,Y2'
   ENDIF ELSE BEGIN
       smash_x = djs_avsigclip(image, 2)
       djs_iterstat, smash_x, mean = mean, median = median, sigma = sigma $
                     , sigrej = 2.0
       slitinds = WHERE(smash_x GT mean, nind)
       smash_y = djs_avsigclip(image[slitinds, *], 1)
       smax = max(smash_y, y_max)
       y1 = y_max-0.05*ny > 0L
       y2 = y_max+0.05*ny < (ny-1L)
;       lightinds = WHERE(smash_y GT 0.6*max(smash_y))
;       y1 = min(lightinds) >  0L
;       y2 = max(lightinds) <  (ny-1L)
       splog, 'Flat field image is illuminated between y1=' + $
              strcompress(string(y1), /rem) + ' and y2=' + $
              strcompress(string(y2), /rem)
   ENDELSE 
   
   IF KEYWORD_SET(XSTART1) AND KEYWORD_SET(XEND1) THEN BEGIN 
      xstart = xstart1
      xend = xend1
      nslit = n_elements(xstart)
      splog, 'Using nslit=', strcompress(string(nslit), /rem)$
             , ' input slit positions for slit tracing'
   ENDIF ELSE BEGIN
      ;; Smash center of the sharpness-filtered image to a 1-dimensional vector
      ymid = 0.5 * (y1 + y2)
      fsum = total(cimg[*, y1:y2], 2)
      ;; Find candidate starting positions
      xmax = long_find_nminima(-fsum, nfind = 1, minsep = minsep, width = 3)
      peakcut = peakthresh*interpolate(-fsum, xmax)
      xstart = long_find_nminima(-fsum < peakcut, nfind = nfind $
                                 , minsep = minsep $
                                 , width = 3)
      
      nstart = n_elements(xstart)
      ;; Find candidate ending positions
      xmax = long_find_nminima(fsum, nfind = 1, minsep = minsep, width = 3)
      peakcut = peakthresh*interpolate(fsum, xmax)
      xend = long_find_nminima(fsum < peakcut, nfind = nfind, minsep = minsep $
                               , width = 3)
      nend = n_elements(xend)
      peakval = 1
      ;stop
      if (nstart GT 0) then peakval = peakval > max(fsum[xstart])
      if (nend GT 0) then peakval = peakval > max((-fsum[xend]))
      if (nstart GT 0) then begin
         igood = where(fsum[xstart] GT peakthresh*peakval, ngood)
         if (ngood GT 0) then xstart = xstart[igood] $
         else xstart = [0]
      endif else xstart = [0]
      if (nend GT 0) then begin
         igood = where(fsum[xend] LT -peakthresh*peakval, ngood)
         if (ngood GT 0) then xend = xend[igood] $
         else xend = [nx-1]
      endif else xend = [nx-1]
      
      xstart = xstart[sort(xstart)]
      xend   = xend[sort(xend)]
      ;; Make sure we start with a "starting" position, and end with an "ending"
      if (xstart[0] GT xend[0]) then $
         xstart = [0, xstart]
      if (xend[n_elements(xend)-1] LT xstart[n_elements(xstart)-1]) then $
         xend = [xend, nx-1]
      ;; Find pairs of starting+ending slit positions, discarding pairs
      ;; of start or end positions that should not exist.
      ;; For book-keeping, set ALLQ=+1 for the start of a slit, =-1 for the end.
      allx = float([xstart, xend])
      allq = [replicate(1, n_elements(xstart)), replicate(-1, n_elements(xend))]
      isort = sort(allx)
      allx = allx[isort]
      allq = allq[isort]
      ;;JOSH'S TEST TO MAKE SURE NO SLITS ARE MISSED:
      ;;IF THERE ARE TWO STARTS FOLLOWED BY TWO ENDS (WHICH SEEMS TO BE WHAT
      ;;HAPPENS SOMETIMES WHEN TWO SLITS GET COMBINED INTO ONE), REVERSE THE SECOND
      ;;START AND THE FIRST END
      jgood = where( (allq EQ 1 and allq EQ shift(allq, -1) and allq NE $
                      shift(allq, -2) and allq NE shift(allq, -3)), njgood)
      if njgood gt 0 then begin
         allq[jgood+2] = 1
         allq[jgood+1] = -1
         allx[jgood+1] = allx[jgood+2]
         allx[jgood+2] = allx[jgood+2]+2
         
         isort = sort(allx)
         allx = allx[isort]
         allq = allq[isort]
      endif

      
;   ; Keep the last starting position if there are several in a row,
;   ; and the first ending position if there are several in a row.
;   igood = where((allq EQ 1 AND allq NE shift(allq,-1)) $
;    OR (allq EQ -1 AND allq NE shift(allq,1)), ngood)
      ;; Keep the first starting position if there are several in a row,
      ;; and the last ending position if there are several in a row.
      allind = lindgen(n_elements(allq))
      igood = where(((allq EQ 1 AND allq NE shift(allq, 1)) OR allind EQ 0) $
                    OR ((allq EQ -1 AND allq NE shift(allq, -1)) OR $
                        allind EQ n_elements(allq)-1L), $
                    ngood, complement = bad)
;IF THE PREVIOUS LINE FOUND ANY INSTANCES OF CONSECUTIVE STARTING POSITIONS,
;GENERATE CORRESPONDING END POSITIONS SO THAT EACH START HAS AN END
;(AND LIKEWISE FOR CONSECUTIVE ENDING POSITIONS)
;THIS AMOUNTS TO MAKING THE DEFAULT ASSUMPTION THAT EVERY STARTING/ENDING 
;POSITION IDENTIFIED EARLIER IN THIS ROUTINE IS REAL, BUT THAT SOME MAY
;BE MISSING THE MATCHING END OR START FOR THAT SLIT
      if ngood ne N_elements(allq) then begin 
         newq = allq 
      newx = allx 
      for i = 0, N_elements(bad) - 1 do begin 
         print, i 
         if (newq[bad[i]-1] eq 1 AND newq[bad[i]] eq -1 AND $
             newq[bad[i]+1] eq -1) then begin 
            newq = [newq[0:bad[i]], 1, allq[(bad[i]+1-i):*]] 
            newx = [newx[0:bad[i]], -999, allx[(bad[i]+1-i):*]] 
            newx[bad[i]] = allx[bad[i]-i] - 3 
            newx[bad[i]+1] = allx[bad[i]-i] + 2
            bad[i:(N_elements(bad)-1)] = bad[i:(N_elements(bad)-1)] + 1 
            print, 'Missing end position -- adding slit end'
         endif 
         if (newq[bad[i]-1] eq 1 AND newq[bad[i]] eq 1 AND $
             newq[bad[i]+1] eq -1) then begin 
            newq = [newq[0:(bad[i]-1)], -1, allq[(bad[i]-i):*]] 
            newx = [newx[0:(bad[i]-1)], -999, allx[(bad[i]-i):*]] 
            newx[bad[i]] = allx[bad[i]-i] - 3 
            newx[bad[i]+1] = allx[bad[i]-i] + 2
            bad[i:(N_elements(bad)-1)] = bad[i:(N_elements(bad)-1)] + 1 
            print, 'Missing start position -- adding slit start'
         endif 
      endfor 
      allq = newq 
      allx = newx 
;<<<<<<< long_slitmask.pro
      ;stop
;=======
;>>>>>>> 1.29
   endif
      allind = lindgen(n_elements(allq))
      igood0 = (allq EQ 1 AND allq NE shift(allq, 1) AND $
                allind NE 0 AND allind NE n_elements(allq)-1L) 
      igood1 = (allq EQ -1 AND allq NE shift(allq, -1) AND $
                allind NE 0 AND allind NE n_elements(allq)-1L) 
      igood2 = allind EQ 0 AND allq EQ 1
      igood3 = allind EQ n_elements(allq)-1L AND allq EQ -1
      igood = where(igood0 OR igood1 OR igood2 OR igood3, ngood)
      
      allx = allx[igood]
      xstart = allx[lindgen(ngood/2)*2]
      xend = allx[lindgen(ngood/2)*2+1]
      nslit = n_elements(xstart)
      
      ;; Trim to slits with width > 5 pixels
      narrow = where(xend - xstart gt 5)
      xstart = xstart[narrow]
      xend = xend[narrow]
      nslit = N_elements(xstart)
   ENDELSE
     ; stop
   ycen = (y1+y2)/2.0D
   radius=2
   ; Trace the starting+end positions on the CCD.
   minval = -djsig(cimg) ; Clip negative values on image when tracing
   xpos1 = trace_crude(cimg>minval, xstart=xstart>0.5, $
    ystart=ycen, radius=radius, $
    nave=nave, maxshifte=maxshifte, maxshift0=maxshift0, xerr=xerr1)
   xpos2 = trace_crude(-cimg>minval, xstart=xend<(nx-0.5), $
    ystart=ycen, radius=radius, $
    nave=nave, maxshifte=maxshifte, maxshift0=maxshift0, xerr=xerr2)
;; Fit these initial traces by a polynomial function
   xy2traceset, rebin(findgen(ny), ny, nslit), xpos1, sset1 $
                , func = func, ncoeff = ncoeff, yfit = yfit1, /silent
   xy2traceset, rebin(findgen(ny), ny, nslit), xpos2, sset2 $
                , func = func, ncoeff = ncoeff, yfit = yfit2, /silent
   slitcen = (yfit1 + yfit2)/2.0d
   flat_spec = extract_boxcar(image, slitcen, rebin(findgen(ny), ny, nslit) $
                              , radius = 2.5)
   ;;ycen = fltarr(nslit)
   yvec = lindgen(ny)
   xposivar = fltarr(ny, nslit)
   FOR kk = 0L, nslit-1L DO BEGIN
      flat_now = djs_median(flat_spec[*, kk], width = 100, boundary = 'reflect')
      max_now = max(flat_now, kmax)
      ;;ycen[kk] = float(kmax)
      keep = where(flat_now GT 0.10*max_now, nkeep)
      xposivar[keep, kk] = 1.0d
   ENDFOR
   ;; re-trace using better starting position
   xpos1 = trace_crude(cimg > minval, xstart = xstart > 0.5 $
                       , ystart = ycen, radius = radius $
                       , nave = nave, maxshifte = maxshifte $
                       , maxshift0 = maxshift0, xerr = xerr1)
   xpos2 = trace_crude(-cimg > minval, xstart = xend < (nx-0.5) $
                       , ystart = ycen, radius = radius $
                       , nave = nave, maxshifte = maxshifte $
                       , maxshift0 = maxshift0, xerr = xerr2)
   ;; Special-case for LRISBLUE, whose images are the join of two CCDs
   ;; with a physical gap between them.
   ;; Also a split the LRIS Red upgrade CCD. Currently doing this
   ;; because of CTE problems but may not be necessary once that is
   ;; fixed. 
   if strcmp(instrument, 'LRISBLUE') OR $
      (strcmp(instrument, 'LRIS') AND $
       strmid(sxpar(hdr[*, 0], 'DATE'), 10) GT '2009-07-01') THEN BEGIN
      xmid = nx/2. + 0.5
      i = (where(xstart LT xmid AND xend GT xmid, ct))[0]
      if (ct EQ 1) then begin
         splog, 'Splitting the slit at the LRISBLUE chip boundary'
         xpos1add = fltarr(ny)+xmid+0.5 ; added slit ending position
         xpos2add = fltarr(ny)+xmid-0.5 ; added slit starting position
         lucifer_slitmask_add1, i+1, xpos1, xerr1, xpos1add, 0*xpos1add+1.
         lucifer_slitmask_add1, i, xpos2, xposivar, xpos2add, 0*xpos2add+1.
         nslit = nslit + 1L
      endif
  endif
   ;stop
   ;; Trim slits
   if keyword_set(TRIM_SEDGE) and not keyword_set(EDIT_SEDGE_FIL) then begin
      print, 'long_slitmask:  Trimming each edge by ', trim_sedge, 'pixels'
      xpos1 = xpos1 + TRIM_SEDGE
      xpos2 = xpos2 - TRIM_SEDGE
   endif

   if keyword_set(COMBINE_SLITS) then begin
      n_remove = n_elements(COMBINE_SLITS)
      slitlist = lindgen(nslit)
      nslit = nslit - n_remove
      newxpos1 = fltarr(ny,nslit)
      newxpos2 = fltarr(ny,nslit)
      newxposivar = fltarr(ny,nslit)
     
      for nn=0L,n_remove-1L do begin
         if nn eq 0L then begin
            slitind = where(slitlist lt COMBINE_SLITS[nn],ct)
            newxpos1[*,slitind] = xpos1[*,slitind]
            newxpos2[*,slitind[0:ct-2]] = xpos2[*,slitind[0:ct-2]]
            newxposivar[*,slitind] = xposivar[*,slitind]
            newxpos2[*,slitind[ct-1]] = xpos2[*,max(slitind)+1]
            ;stop
         endif else begin
            slitind = where(slitlist gt COMBINE_SLITS[nn-1] and $
                            slitlist lt COMBINE_SLITS[nn],ct)
            newxpos1[*,slitind-nn] = xpos1[*,slitind]
            newxpos2[*,slitind[0:ct-2]-nn] = xpos2[*,slitind[0:ct-2]]
            newxposivar[*,slitind-nn] = xposivar[*,slitind]
            newxpos2[*,slitind[ct-1]-nn] = xpos2[*,max(slitind)+1]
            ;stop
         endelse 
            
        
      endfor 
      slitind = where(slitlist gt COMBINE_SLITS[n_remove-1],ct)
      if ct gt 0 then begin
         newxpos1[*,slitind-n_remove] = xpos1[*,slitind]
         newxpos2[*,slitind-n_remove] = xpos2[*,slitind]
         newxposivar[*,slitind-n_remove] = xposivar[*,slitind]
         ;stop
      endif 
      xpos1 = newxpos1
      xpos2 = newxpos2
      xposivar = newxposivar
   endif 
         

   if keyword_set(SPLIT_SLITS) then begin
      readcol, SPLIT_SLITS, spl_add_end, spl_add_start, format='F,F'
      medianxpos1 = djs_median(xpos1,1)
      medianxpos2 = djs_median(xpos2,1)
      spl_add_nslit = n_elements(spl_add_end)
      print, 'long_slitmask: Splitting ', spl_add_nslit, $
             ' slits from file ', SPLIT_SLITS
      nslit = nslit + spl_add_nslit
      newxpos1 = fltarr(ny,nslit)
      newxpos2 = fltarr(ny,nslit)
      newxposivar = fltarr(ny,nslit)
      ;stop
      for mm=0L,spl_add_nslit-1 do begin
         if mm eq 0L then begin
            slitind = where(medianxpos1 lt spl_add_end[mm],ct)
         endif else begin
            slitind = where(medianxpos1 lt spl_add_end[mm] and $
                            medianxpos1 ge spl_add_end[mm-1],ct)
         endelse 
         ;stop
         newxpos1[*,slitind+mm] = xpos1[*,slitind]
         newxpos2[*,slitind[0:ct-2]+mm] = xpos2[*,slitind[0:ct-2]]
         newxposivar[*,slitind+mm] = xposivar[*,slitind]

         selslit = max(slitind)
         newxpos2[*,selslit+mm] = xpos1[*,selslit] $
                                    + (spl_add_end - medianxpos1[selslit])
         newxpos1[*,selslit+1+mm] = xpos1[*,selslit] $
                                    + (spl_add_start - medianxpos1[selslit])
         newxpos2[*,selslit+1+mm] = xpos2[*,selslit]
         newxposivar[*,selslit+1+mm] = xposivar[*,selslit] 
         ;stop
      endfor 
      slitind = where(medianxpos1 ge spl_add_end[spl_add_nslit-1],ct)
      ;stop
      if ct gt 0 then begin
         newxpos1[*,slitind+spl_add_nslit] = xpos1[*,slitind]
         newxpos2[*,slitind+spl_add_nslit] = xpos2[*,slitind]
         newxposivar[*,slitind+spl_add_nslit] = xposivar[*,slitind]
      endif 
      xpos1 = newxpos1
      xpos2 = newxpos2
      xposivar = newxposivar
      
   endif

   if keyword_set(ADD_SLITS) then begin
      readcol, ADD_SLITS, add_start, add_end, format='F,F'
      medianxpos1 = djs_median(xpos1,1)
      medianxpos2 = djs_median(xpos2,1)
      add_nslit = n_elements(add_start)
      print, 'long_slitmask: Adding ', add_nslit, ' slits from file', ADD_SLITS
      nslit = nslit + add_nslit
      newxpos1 = fltarr(ny,nslit)
      newxpos2 = fltarr(ny,nslit)
      newxposivar = fltarr(ny,nslit)
      ;slit_reloc = lonarr(nslit)
      for mm=0L,add_nslit-1 do begin
         if mm eq 0L then begin
            slitind = where(medianxpos2 lt add_start[mm])
            ;stop
         endif else begin
            slitind = where(medianxpos2 lt add_start[mm] and medianxpos2 $
                            ge add_start[mm-1])
            ;stop
         endelse 
         newxpos1[*,slitind+mm] = xpos1[*,slitind]
         newxpos2[*,slitind+mm] = xpos2[*,slitind]
         newxposivar[*,slitind+mm] = xposivar[*,slitind]

         selslit = max(slitind)       
         newxpos1[*,selslit+1+mm] = xpos2[*,selslit]+3.0
         newxpos2[*,selslit+1+mm] = xpos1[*,selslit+1]-3.0
         newxposivar[*,selslit+1+mm] = xposivar[*,selslit]
      endfor 
      ;; Fill in the end
      slitind = where(medianxpos2 ge add_start[add_nslit-1],ct)
      if ct gt 0 then begin
         newxpos1[*,slitind+add_nslit] = xpos1[*,slitind]
         newxpos2[*,slitind+add_nslit] = xpos2[*,slitind]
         newxposivar[*,slitind+add_nslit] = xposivar[*,slitind]
      endif 
      xpos1 = newxpos1
      xpos2 = newxpos2
      xposivar = newxposivar
   endif 
   ;stop

  ;; Trim specific slits by hand
  if keyword_set(EDIT_SEDGE_FIL) then begin
     if keyword_set(TRIM_SEDGE) then begin
        print, 'long_slitmask:  Trimming each edge by ', trim_sedge, 'pixels'
        xpos1 = xpos1 + TRIM_SEDGE
        xpos2 = xpos2 - TRIM_SEDGE
     endif 
     mtchtol = 8.0
     print, 'long_slitmask: Reading in file to trim specific slits: ', edit_sedge_fil
     print, '               with tolerance ', mtchtol, '   pixels'
     readcol, edit_sedge_fil, hand_xmid, hand_flg, hand_npix, hand_lower, $
              hand_upper, format='F,L,F,F,F'
     for kk=0L,nslit-1 do begin
        tmp_xpos2 = xpos2[*,kk]
        tmp_xpos1 = xpos1[*,kk]
        ;slitsec = where(findgen(ny) ge hand_lower[0] and findgen(ny) le hand_upper[0])
        tmp_xposmid = (tmp_xpos2+tmp_xpos1)/2.
        tmp_xposavg = djs_mean(tmp_xposmid)
        mtchbadslit = where(hand_xmid gt tmp_xposavg-mtchtol and $
                            hand_xmid lt tmp_xposavg+mtchtol, ctmtch)
        
        CASE ctmtch OF 
           1L: BEGIN
              slitsec = where(findgen(ny) ge hand_lower[mtchbadslit[0]] and findgen(ny) $
                                 le hand_upper[mtchbadslit[0]])
              if hand_flg[mtchbadslit] eq 1L then begin
                 if keyword_set(TRIM_SEDGE) then tmp_xpos1 = tmp_xpos1 - TRIM_SEDGE
                 xpos1[slitsec,kk] = tmp_xpos1[slitsec]+hand_npix[mtchbadslit[0]]
                 print, 'long_slitmask: editing left edge of slit at xpos = ', tmp_xposavg
              endif else begin
                 if hand_flg[mtchbadslit] eq 2L then begin
                    if keyword_set(TRIM_SEDGE) then tmp_xpos2 = tmp_xpos2 + TRIM_SEDGE
                    xpos2[slitsec,kk] = tmp_xpos2[slitsec]-hand_npix[mtchbadslit[0]]
                    print, 'long_slitmask: editing right edge of slit at xpos = ', tmp_xposavg
                 endif else begin
                    print, 'long_slitmask: Bad slit side specification (should be 1 or 2)'
                    stop
                 endelse 
              endelse 
              ;stop
           END
           0L: BEGIN
              print, 'You dont seem to want to change the slit at xpos = ', tmp_xposavg
              ;stop
           END 
           ELSE: print, 'long_slitmask: Ambiguity in slit selection-- you may want to lower mtchtol'
        ENDCASE
     endfor 
  endif 
;stop
  ;; Deal w/ overlapping slits
  ;if keyword_set(REMOVE_OVERLAP) then begin
     ;stop
  ;   for kk=0L,nslit-2 do begin
  ;      tmp_xpos2 = xpos2[*,kk]
  ;      tmp_xpos1 = xpos1[*,kk+1]
  ;      overlap = where(tmp_xpos2 gt tmp_xpos1,n_overlap,complement=nonoverlap)
  ;      if n_overlap gt 0 and n_overlap lt 0.7*ny then begin
  ;         ;stop
  ;         print, 'long_slitmask: Removing overlap between slits ', kk+1, ' and ', kk+2
  ;         tmp_xpos2[overlap] = mean(tmp_xpos2[nonoverlap])
  ;         tmp_xpos1[overlap] = mean(tmp_xpos1[nonoverlap])
  ;         xpos2[*,kk] = tmp_xpos2
  ;         xpos1[*,kk+1] = tmp_xpos1
  ;      endif 
  ;   endfor 
  ;endif
  ;; Fit these traces by a polynomial function
  ;stop
  xy2traceset, rebin(findgen(ny), ny, nslit), xpos1, sset1 $
               , func = func, ncoeff = ncoeff, invvar = xposivar, yfit = yfit1
  xy2traceset, rebin(findgen(ny), ny, nslit), xpos2, sset2 $
               , func = func, ncoeff = ncoeff, invvar = xposivar, yfit = yfit2
  tset_slits = [sset1, sset2]
  tset_slits = struct_addtags(tset_slits $
                              , replicate(create_struct('DIMS', dims) $
                                          , size(tset_slits, /dimens)))
  ;stop
  ; Trim to slits with width > MINSLIT pixels
   IF KEYWORD_SET(MINSLIT) THEN BEGIN
       old_slits = tset_slits
       width = old_slits[1].coeff[0, *]-old_slits[0].coeff[0, *]
       keep_slits = WHERE(width GE MINSLIT, nkeep)
       ;stop
       IF nkeep EQ 0 THEN begin
           print, 'ERROR: No good slits on this mask'
           print, 'If you continue, the code will make one slit '
           print, 'covering the entire chip.  This is often the case'
           print, 'for MMT and KAST data.'
           print, 'In this case, you should include a standard star in'
           print, 'the reductions.'
           nkeep = 1
           tset_proto = $
             { func    :    old_slits[0].FUNC, $
               xmin    :    old_slits[0].XMIN, $
               xmax    :    old_slits[0].XMAX, $
               dims    :    old_slits[0].DIMS, $
               coeff   :    dblarr(3, nkeep),  $
               xcorr_coeff : dblarr(3, nkeep)  $
             }
           tset_slits = replicate(tset_proto, 2)
           tset_slits[0].coeff[0, 0] = 5.0D
           tset_slits[1].coeff[0, 0] = (float(nx-1)-5.0) ;; JXP (handles trim)
       endif else begin
           tset_proto = $
             { func    :    old_slits[0].FUNC, $
               xmin    :    old_slits[0].XMIN, $
               xmax    :    old_slits[0].XMAX, $
               coeff   :    dblarr(3, nkeep),  $
               dims    :    old_slits[0].DIMS, $
               xcorr_coeff : dblarr(3, nslit)  $
             }
           tset_slits = replicate(tset_proto, 2)
           tset_slits[0].COEFF[*, *] = old_slits[0].COEFF[*, keep_slits]
           tset_slits[1].COEFF[*, *] = old_slits[1].COEFF[*, keep_slits]
           tset_slits[0].XCORR_COEFF = old_slits[0].COEFF
           tset_slits[1].XCORR_COEFF = old_slits[1].COEFF
       endelse
       nslit = nkeep
   ENDIF
ENDELSE
   ;; Construct this image to write to HDU #0 of the output file
   slitmask = long_slits2mask(tset_slits, nslit = nslit)
   
   ;; Write output file
   IF KEYWORD_SET(OUTFILE) THEN BEGIN
       splog, 'Writing output file'
       mwrfits, slitmask, outfile, hdr[*, 0], /create
       mwrfits, tset_slits, outfile
   ENDIF
   splog, 'Number of slits = ', nslit
   splog, 'Elapsed time = ', systime(1)-t0, ' sec'
;   return

   traceset2xy, tset_slits[0], rows, left_edge
   traceset2xy, tset_slits[1], rows, right_edge
   
end
;------------------------------------------------------------------------------
