;+
; NAME:
;   gmos_reduce
;
; PURPOSE:
;
;   Main program for the Low-redux pipeline.  This set of algorithms
;   runs mainly as a black box.
;
; CALLING SEQUENCE:
;  long_reduce, planfile, /clobber, /NOZAP, /NOFLEX, /NOHELIO 
;
; INPUTS:
;  planfile  -- File created by long_plan which guides the reduction
;               process
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
;  (1) Various calibration files
;  (2) One multi-extension FITS file in Science per exposure containing
;  the extracted data and processed images
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
;   11-Mar-2005  Written by JH + SB
;-  
;-----------------------------------------------------------------------------
PRO gmos_reduce, planfile, clobber = clobber, verbose = verbose $
                 , NOFLEX = noflex1, NOZAP = NOZAP1, NOHELIO = NOHELIO1 $
                 , PROF_NSIGMA = PROF_NSIGMA1, SKYTRACE = SKYTRACE1 $
                 , NOSHIFT = NOSHIFT1, STD = STD1 $
                 , PSF_CR_FUDGE = PSF_CR_FUDGE1, NSIGMA_CR = NSIGMA_CR1 $
                 , ONLYSCI = ONLYSCI, HAND_X = HAND_X, HAND_Y = HAND_Y $
                 , HAND_FWHM = HAND_FWHM, LINELIST = linelist $
                 , REID_FILE= REID_FILE, TRIM_SEDGE=trim_sedge $
                 , REMOVE_OVERLAP=remove_overlap $
                 , USAGE = USAGE, FILESTD = FILESTD1, CHK = CHK $
                 , TRCCHK = TRCCHK, BIN_RATIO = bin_ratio $
                 , ISLIT = ISLIT1, EDIT_SEDGE_FIL=EDIT_SEDGE_FIL $
                 , ADD_SLITS = ADD_SLITS, SPLIT_SLITS = SPLIT_SLITS $
                 , COMBINE_SLITS = COMBINE_SLITS $
                 , ARC_INTER = arc_inter, PROCARC = procarc $
                 , TWEAK_ARC = tweak_arc, PIXFLAT_ARCHIVE = pixflat_archive $
                 , _EXTRA = extra 
  if  KEYWORD_SET(USAGE) THEN BEGIN
      print,'Syntax - ' + $
        'long_reduce, planfile, /CLOBBER, /VERBOSE, /NOFLEX, /NOHELIO [v1.0]' 
      return
  endif 

if (NOT keyword_set(planfile)) then planfile = findfile('plan*.par')

   ;----------
   ; If multiple plan files exist, then call this script recursively
   ; for each such plan file.

   if planfile[0] EQ '' then begin
      print, 'could not find plan file'
      print, 'try running long_plan'
      return
   endif

   if (n_elements(planfile) GT 1) then begin
       FOR i = 0L, n_elements(planfile)-1L do BEGIN
           splog, '------------------------'
           splog, 'Reducing ' + planfile[i]
           splog, '------------------------'
           gmos_reduce, planfile[i], clobber = clobber, verbose = verbose $
                       , NOFLEX = noflex1, NOZAP = NOZAP1, NOHELIO = NOHELIO1 $
                       , PROF_NSIGMA = PROF_NSIGMA1, SKYTRACE = SKYTRACE1 $
                       , NOSHIFT = NOSHIFT1, STD = STD1 $
                       , PSF_CR_FUDGE = PSF_CR_FUDGE1, NSIGMA_CR = NSIGMA_CR1 $
                       , ONLYSCI = ONLYSCI, HAND_X = HAND_X, HAND_Y = HAND_Y $
                       , HAND_FWHM = HAND_FWHM, LINELIST = linelist $
                       , REID_FILE = REID_FILE $
                       , USAGE = USAGE, FILESTD = FILESTD1, CHK = CHK $
                       , TRCCHK = TRCCHK, BIN_RATIO = bin_ratio $
                       , ISLIT = ISLIT1 $
                       , ARC_INTER = arc_inter, PROCARC = procarc $
                       , TWEAK_ARC = tweak_arc, PIXFLAT_ARCHIVE = pixflat_archive $
                       , _EXTRA = extra 
       ENDFOR
       RETURN
   endif

   ;----------
   ; Read the plan file
   planstr = yanny_readone(planfile, hdr=planhdr, /anonymous)
   if (NOT keyword_set(planstr)) then begin
       splog, 'Empty plan file ', planfile
       return
   endif

   ;; Truncate the gz stuff
   nfil = n_elements(planstr)
   for qq=0L,nfil-1 do begin
       slen = strlen(planstr[qq].filename)
       if strmid(planstr[qq].filename,slen-3) EQ '.gz' then $
         planstr[qq].filename = strmid(planstr[qq].filename,0,slen-3)
   endfor
   logfile = yanny_par(planhdr, 'logfile')
   plotfile = yanny_par(planhdr, 'plotfile')
   indir = yanny_par(planhdr, 'indir')
   tempdir = yanny_par(planhdr, 'tempdir')
   scidir  = yanny_par(planhdr, 'scidir')
   minslit1 = long(yanny_par(planhdr, 'minslit'))
   slitthresh1 =  float(yanny_par(planhdr, 'slitthresh'))
   reduxthresh = float(yanny_par(planhdr, 'reduxthresh'))
   sig_thresh = float(yanny_par(planhdr, 'sig_thresh'))
   nolocal = long(yanny_par(planhdr, 'nolocal'))
   slity1_1 = float(yanny_par(planhdr, 'slity1'))
   slity2_1 = float(yanny_par(planhdr, 'slity2'))
   maxflex1 = float(yanny_par(planhdr, 'maxflex'))
   maxgood = long(yanny_par(planhdr, 'maxgood'))
;   ksize_1  = yanny_par(planhdr, 'slitksize')
   maxobj = long(yanny_par(planhdr, 'maxobj'))
   box_rad1 = float(yanny_par(planhdr, 'box_rad'))
   IF n_elements(SKYTRACE1) GT 0 THEN SKYTRACE = SKYTRACE1 $
   ELSE skytrace = long(yanny_par(planhdr, 'skytrace'))
   IF n_elements(NOFLEX1) GT 0 THEN NOFLEX = NOFLEX1 $
   ELSE noflex = long(yanny_par(planhdr, 'noflex'))
   IF n_elements(NOHELIO1) GT 0 THEN NOHELIO = NOHELIO1 $
   ELSE nohelio = long(yanny_par(planhdr, 'nohelio'))
   IF n_elements(NOSHIFT1) GT 0 THEN NOSHIFT = NOSHIFT1 $
   ELSE noshift = long(yanny_par(planhdr, 'noshift'))
   IF n_elements(PROF_NSIGMA1) GT 0 THEN PROF_NSIGMA = PROF_NSIGMA1 $
   ELSE prof_nsigma = float(yanny_par(planhdr, 'prof_nsigma'))
   IF n_elements(STD1) GT 0 THEN STD = STD1 $
   ELSE STD = long(yanny_par(planhdr, 'std'))
   IF n_elements(PSF_CR_FUDGE1) GT 0 THEN PSF_CR_FUDGE = PSF_CR_FUDGE1 $
   ELSE PSF_CR_FUDGE = float(yanny_par(planhdr, 'psf_cr_fudge'))
   IF n_elements(NOZAP1) GT 0 THEN NOZAP = NOZAP1 $
   ELSE NOZAP = float(yanny_par(planhdr, 'nozap'))
   IF n_elements(NSIGMA_CR1) GT 0 THEN NSIGMA_CR = NSIGMA_CR1 $
   ELSE NSIGMA_CR = float(yanny_par(planhdr, 'nsigma_cr'))
   IF n_elements(ISLIT1) GT 0 THEN ISLIT = ISLIT1 $
   ELSE ISLIT = long(yanny_par(planhdr, 'islit'))

   IF KEYWORD_SET(slity1_1) THEN slity1 = long(slity1_1)
   IF KEYWORD_SET(slity2_1) THEN slity2 = long(slity2_1)
;   IF KEYWORD_SET(ksize_1) THEN ksize = long(ksize_1)
   IF KEYWORD_SET(box_rad1) THEN box_rad = long(box_rad1)
   IF KEYWORD_SET(minslit1) THEN minslit = long(minslit1)
   IF KEYWORD_SET(maxflex1) THEN maxflex = long(maxflex1)
   
   ;----------
   ; Create science dir
   IF keyword_set(scidir) THEN spawn, '\mkdir -p '+scidir

   ;----------
   ; Open log file
   if (keyword_set(logfile)) then begin
       splog, filename = logfile
       splog, 'Log file ' + logfile + ' opened ' + systime()
   endif
   splog, 'IDL version: ' + string(!version,format='(99(a," "))')
   spawn, 'uname -a', uname
   splog, 'UNAME: ' + uname[0]
   splog, 'idlutils version ' + idlutils_version()
   splog, 'Longslit version ' + longslit_version()
   
   plotfile = 0
   if (keyword_set(plotfile)) then begin
       thisfile = findfile(plotfile, count = ct)
       IF (ct EQ 0 OR KEYWORD_SET(CLOBBER)) THEN BEGIN
           splog, 'Plot file ' + plotfile
           dfpsplot, plotfile, /color
       ENDIF ELSE BEGIN
           cpbackup, plotfile
           splog, 'Plot file already exists. Creating backup'
           splog, 'Plot file ' + plotfile
           dfpsplot, plotfile, /color
       ENDELSE
   ENDIF
      
   ;----------
   ; Loop over each INSTRUMENT (e.g., CCD)

   ccd_list = planstr.instrument
   ccd_list = ccd_list[uniq(ccd_list, sort(ccd_list))]
   nccd = n_elements(ccd_list)

   for iccd=0L, nccd-1L do begin
      indx = where(planstr.instrument EQ ccd_list[iccd])

      ;;-----------------------------
      ;; Make a superbias if possible
      ;;-----------------------------

      ii = where(planstr[indx].flavor EQ 'bias', nbias)
      if (nbias GT 0) then begin
          ibias = indx[ii]
          superbiasfile = 'superbias-' + planstr[ibias[0]].filename
          thisfile = findfile(superbiasfile, count = ct)
          if (ct EQ 0 OR keyword_set(clobber)) then begin
              splog, 'Generating superbias for INSTRUMENT=', ccd_list[iccd]
              long_superbias $
                , djs_filepath(planstr[ibias].filename $
                               , root_dir = indir) $
                , superbiasfile, verbose = verbose
          endif else begin
              splog, 'Do not overwrite existing superbias ', thisfile
          endelse
      endif else begin
          superbiasfile = ''
          splog, 'No input biases for superbias for INSTRUMENT=', $
                 ccd_list[iccd]
      endelse
  
      ; Loop over each GRATING+MASK+WAVE for this INSTRUMENT

      mask_list = planstr[indx].grating 
; + planstr[indx].wave + planstr[indx].maskname $
      mask_list = mask_list[uniq(mask_list, sort(mask_list))]
      nmask = n_elements(mask_list)
      
      for imask = 0L, nmask-1L do begin
          jndx = indx[where(planstr[indx].grating  EQ mask_list[imask]) ]
          
          qboth = planstr[jndx].flavor EQ 'bothflat'
          qtwi = (planstr[jndx].flavor EQ 'twiflat') OR qboth
          qpix = (planstr[jndx].flavor EQ 'domeflat' OR $
                   planstr[jndx].flavor EQ 'iflat') OR qboth
;          iflat = where(qtwi OR qdome)
          itwi  = where(qtwi, ntwi)
          ipix  = where(qpix, npix)
          iboth = WHERE(qtwi OR qpix, nboth)
         ;;---------------------------
         ;; Create slitmask file 
         ;;---------------------------
          slitfile = 'slits-GMOS.fits'
          thisfile = findfile(slitfile, count = ct)
          if (ct EQ 0 OR keyword_set(clobber)) then begin
              splog, 'Generating slits for INSTRUMENT=', ccd_list[iccd], $
                     ' GRATING+MASK+WAVE=', mask_list[imask]
              long_slitmask, djs_filepath(planstr[jndx[0]].filename, $
                                          root_dir = indir) $
                             , slitfile, verbose = verbose
          endif else begin
              splog, 'Do not overwrite existing slitmask file ', $
                     thisfile
          endelse
         ;---------------------------
         ; Make a wavelength solution
         ;---------------------------
         ii = where(planstr[jndx].flavor EQ 'arc', narc)

         if (narc GT 0) then begin
             iarc = jndx[ii]
             wavefile = 'wave-' + planstr[iarc[0]].filename
             if strpos(wavefile, '.fits') LT 0 then wavefile = wavefile+'.fits'
             thisfile = findfile(wavefile, count = ct)
             if (ct EQ 0 OR keyword_set(clobber)) then begin
                 splog, 'Generating wavelengths for INSTRUMENT=' $
                        , ccd_list[iccd], $
                        ' GRATING+MASK+WAVE=', mask_list[imask]
                 fitfile = 'database/fcgs' + $
                   repstr(planstr[iarc[0]].filename, '.fits', '_001')
                 thisfitfile = findfile(fitfile, count = fct)
                 idfile = 'database/idgs' + $
                   repstr(planstr[iarc[0]].filename, '.fits', '_001')
                 thisidfile = findfile(fitfile, count = idct)
                 gsfile = 'gs' + planstr[iarc[0]].filename
                 thisgsfile = findfile(fitfile, count = gct)
                 IF fct NE 0 AND idct NE 0 AND gct NE 0 THEN BEGIN
                     splog $
                       , 'Found IRAF files for GMOS-N wavelength solution: ' $
                       , thisfitfile, thisidfile, thisgsfile
                     wavehdr = headfits(planstr[iarc[0]].filename)
                     gmos_wave, fitfile, idfile, gsfile, slitfile $ 
                                , wavefile, hdr = wavehdr $
                                , FWHMSET = FWHMSET 
                 ENDIF ELSE message, 'Could not find IRAF database file: ' $
                   + thisidfile
             ENDIF ELSE BEGIN
                 splog, 'Do not overwrite existing wavelength file ', thisfile
             ENDELSE
         ENDIF ELSE BEGIN
             wavefile = ''
             splog, 'No input arcs for wavelengths for INSTRUMENT=', $
                    ccd_list[iccd], ' GRATING+MASK+WAVE=', mask_list[imask]
         ENDELSE
     
         ;;------------------------------------
         ;; Make  pixel and illumination flats 
         ;;------------------------------------
         IF (npix GT 0) THEN BEGIN
             pixflatfile = 'pixflat-' + planstr[jndx[ipix[0]]].filename
             thispixflatfile = findfile(pixflatfile, count = pixct)
             if (pixct EQ 0 OR keyword_set(clobber)) then begin
                 splog, 'Generating pixel flat for INSTRUMENT=' $
                        , ccd_list[iccd], ' GRATING+MASK+WAVE='$
                        , mask_list[imask]
                 gmos_superflat $
                   , djs_filepath(planstr[jndx[ipix]].filename $
                                  , root_dir = indir), pixflatfile $
                   , biasfile = superbiasfile, indir = indir $
                   , tempdir = tempdir $
                   , /write_flats
             ENDIF ELSE BEGIN
                 splog, 'Do not overwrite existing pixel flats ' $
                        , thispixflatfile
             endelse
         endif else begin
             pixflatfile = ''
             splog, 'No input pixel flats for flats for INSTRUMENT=', $
                    ccd_list[iccd], ' GRATING+MASK=', mask_list[imask]
         endelse

         ;;------------------------------
         ;; Make illumination flats if possible
         ;;------------------------------
         
         if (ntwi GT 0) THEN BEGIN 
             illumflatfile = 'illumflat-' + planstr[jndx[itwi[0]]].filename
             thisillumflatfile = findfile(illumflatfile, count = illumct)
             if (illumct EQ 0 OR keyword_set(clobber)) then begin
;                illumination function for gemini
                 splog, 'Generating illumnination flat for INSTRUMENT=' $
                        , ccd_list[iccd], ' GRATING+MASK+WAVE='$
                        , mask_list[imask]
                 gmos_illumflat, djs_filepath(planstr[jndx[itwi]].filename $
                                         , root_dir = indir) $
                            , illumflatfile, pixflatfile $
                            , biasfile = superbiasfile $
                            , wavefile = wavefile, INDIR = INDIR $
                            , TEMPDIR = TEMPDIR $
                            , /write_flats $
                            , verbose = verbose
             endif else begin
                 splog, 'Do not overwrite existing flat ' $
                        , thisillumflatfile
             endelse
         endif else begin
             illumflatfile = ''
             splog, 'No input illumination flats for flats for INSTRUMENT=', $
                    ccd_list[iccd], ' GRATING+MASK=', mask_list[imask]
         endelse

         ;-----------------------------------
         ; Finally, reduce each science image
         ;-----------------------------------

         ii = where(planstr[jndx].flavor EQ 'science', nsci)
         for isci = 0L, nsci-1 do begin
             j = jndx[ii[isci]]
             scifile = djs_filepath('sci-' + planstr[j].filename $
                                    , root_dir = scidir)
             ipos = strpos(scifile, '.fits')
             if ipos LT 0 then scifile = scifile+'.fits'
             IF KEYWORD_SET(PROFILE) THEN $
               profile_filename = 'profile-' + planstr[j].filename $
             ELSE PROFILE_FILENAME = 0
             thisfile = findfile(scifile+'*', count = ct)
             if (ct EQ 0 OR keyword_set(clobber)) THEN BEGIN
                 splog, 'Reducing science frame ', prelog = planstr[j].filename
                 IF KEYWORD_SET(FILESTD1) THEN filestd = filestd1 $
                 ELSE BEGIN
                     istd = where(planstr[jndx].flavor EQ 'std', nstd)
                     IF nstd GT 0 THEN BEGIN
                         stdfile = 'std-' + planstr[istd[0]].filename
                         thisstdfile = findfile('*' + stdfile + '*' $
                                                , count = ct1)
                         IF ct1 GT 0 THEN filestd = thisstdfile[0] 
                     ENDIF 
                  ENDELSE
                 ;; overwrite it filestd was input
                 long_reduce_2dwork, djs_filepath(planstr[j].filename $
                                                , root_dir = indir), scifile $
                                   , slitfile = slitfile, wavefile = wavefile $
                                   , biasfile = superbiasfile $
                                   , pixflatfile = pixflatfile $
                                   , illumflatfile = illumflatfile $
                                   , verbose = verbose $
                                   , maxobj = maxobj, box_rad = box_rad $
                                   , reduxthresh = reduxthresh $
                                   , sig_thresh = sig_thresh $
                                   , noflex = noflex, NOZAP = NOZAP $
                                   , maxflex = maxflex , MAXGOOD=maxgood $
                                   , NOHELIO = NOHELIO $
                                   ;, profile_filename = profile_filename $
                                   , HAND_X = HAND_X, HAND_Y = HAND_Y $
                                   , HAND_FWHM = HAND_FWHM $
                                   , PROF_NSIGMA=PROF_NSIGMA $
                                   , FILESTD = FILESTD $
                                   , STD = (planstr[j].flavor EQ 'std') $
                                   , ISLIT = ISLIT, CHK = CHK $
                                   , TRCCHK = TRCCHK, NOLOCAL=nolocal $
                                   , NOSHIFT = NOSHIFT, SKYTRACE = SKYTRACE $
                                   , _EXTRA = extra
                 splog, prelog = ''
             endif else begin
                 splog, 'Do not overwrite existing science frame ', scifile
             endelse
         endfor
       
   endfor                       ; End loop over GRATING+MASK
   endfor ; End loop over INSTRUMENT

   if (keyword_set(plotfile)) then begin
       dfpsclose
   endif
   
   splog, /close
   x_psclose
   
   return
end
;------------------------------------------------------------------------------
