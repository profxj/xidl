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
pro sinfoni_slitmask, flatfile, outfile, darkfile = darkfile $
                      , xstart = xstart1, xend = xend1 $
                      , minslit = minslit, y1 = y1_in, y2 = y2_in, nmed = nmed $
                      , ksize = ksize, peakthresh = peakthresh $
                      , radius = radius, nave = nave, maxshifte = maxshifte $
                      , maxshift0 = maxshift0, func = func, ncoeff = ncoeff $
                      , verbose = verbose, tset_slits = tset_slits $
                      , nfind = nfind1, GMOSLONG = GMOSLONG $
                      , minsep = minsep1, CCDONLY = CCDONLY, CIMG = CIMG $
                      , TRANSFORM = TRANSFORM, TRIM_SEDGE = trim_sedge $
                                ;, REMOVE_OVERLAP=remove_overlap, $
                      , EDIT_SEDGE_FIL = edit_sedge_fil, ADD_SLITS = add_slits $
                      , SPLIT_SLITS = split_slits $
                      , COMBINE_SLITS = combine_slits $
                      , YRMASK = yrmask $
                      , NO_HEFIXES = no_hefixes, slitmask = slitmask ;;$
  ;;, SLIT_EDG_TWK = SLIT_EDG_TWK
  
  
   t0 = systime(1)

   ;;----------
   ;; Read raw image
   sinfoni_proc, flatfile, image, ivar $
                 , hdr = hdr, verbose = verbose, darkfile = darkfile

   optic = strcompress(esopar(hdr, 'HIERARCH ESO INS OPTI1 NAME'), /rem)
   grating = strcompress(esopar(hdr, 'HIERARCH ESO INS GRAT1 NAME'), /rem)

   ;; This has been tuned to work, but probably we need a bad pixel mask
   ;; to make it more robust, and not identify spurious slits or miss
   ;; slits
   setup = optic + '-' + grating
   CASE setup OF
      '0.025-K': peakthresh = 0.2
      '0.25-K': BEGIN
         peakthresh = 0.05
         combine_slits = [4, 18, 20] ;; This kludge is necessary to fix some problem slits
      END
      ELSE: peakthresh = 0.05 ;; this is kludgy for now
   ENDCASE
   slit_edg_twk = 1.0
   ksize = 2.0
   long_slitmask_work, image, ivar, hdr, outfile, biasfile = biasfile $
                       , xstart = xstart1, xend = xend1 $
                       , minslit = minslit, y1 = y1_in, y2 = y2_in $
                       , nmed = nmed $
                       , ksize = ksize, peakthresh = peakthresh $
                       , radius = radius, nave = nave $
                       , maxshifte = maxshifte $
                       , maxshift0 = maxshift0, func = func, ncoeff = ncoeff $
                       , verbose = verbose, tset_slits = tset_slits $
                       , nfind = nfind1, GMOSLONG = GMOSLONG $
                       , minsep = minsep1, CCDONLY = CCDONLY, CIMG = CIMG $
                       , TRANSFORM = TRANSFORM, TRIM_SEDGE = trim_sedge $
                                ;, REMOVE_OVERLAP=remove_overlap, $
                       , EDIT_SEDGE_FIL = edit_sedge_fil $
                       , ADD_SLITS = add_slits $
                       , SPLIT_SLITS = split_slits $
                       , COMBINE_SLITS = combine_slits $
                       , YRMASK = yrmask $
                       , NO_HEFIXES = no_hefixes, slitmask = slitmask $
                       , SLIT_EDG_TWK = SLIT_EDG_TWK, nslit = nslit

   IF nslit NE 32 THEN message, 'Error in slitmask construction. Need to massage inputs to long_slitmask_work'
   ;; Useful tip: For debugging the slitmask construction, turn off the
   ;; slit_edg_twk to make it faster, then turn it back on once things
   ;; are fixed

   
   splog, 'Elapsed time = ', systime(1)-t0, ' sec'
   return
end
;------------------------------------------------------------------------------
