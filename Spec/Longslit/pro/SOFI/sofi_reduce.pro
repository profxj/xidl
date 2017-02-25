;------------------------------------------------------------------------------
PRO sofi_reduce, planfile, clobber_calib = clobber_calib $
                 , nuclear = nuclear, verbose = verbose $
                 , CHK = CHK, WVCHK = WVCHK, CALIB = CALIB
  
  if (NOT keyword_set(planfile)) then planfile = findfile('plan*.par')
  
  ;;----------
  ;; If multiple plan files exist, then call this script recursively
  ;; for each such plan file.

  if planfile[0] EQ '' then begin
     print, 'ERROR: Could not find plan file'
     print, 'Try running gnirs_plan'
     return
  endif
  
  if (n_elements(planfile) GT 1) then begin
     for i = 0L, n_elements(planfile)-1L do $
        sofi_reduce, planfile[i], clobber_calib = clobber_calib $
                     , nuclear = nuclear, verbose = verbose $
                     , CHK = CHK, WVCHK = WVCHK, CALIB = CALIB
     
     return
  endif

  IF KEYWORD_SET(NUCLEAR) THEN CLOBBER_CALIB = 1
  ;;----------
  ;; Read the plan file
  planstr = yanny_readone(planfile, hdr = planhdr, /anonymous)
  if (NOT keyword_set(planstr)) then begin
     splog, 'Empty plan file ', planfile
     return
  endif
  logfile = yanny_par(planhdr, 'logfile')
  plotfile = yanny_par(planhdr, 'plotfile')
  indir = yanny_par(planhdr, 'indir')
  tempdir = yanny_par(planhdr, 'tempdir')
  scidir  = yanny_par(planhdr, 'scidir')
  use_tell_wave  = yanny_par(planhdr, 'use_tell_wave')
  

  plotfile = 0
  
  ;;----------
  ;; Create science dir
  IF keyword_set(scidir) THEN spawn, '\mkdir -p '+scidir
  
  ;;----------
  ;; Open log file
  if (keyword_set(logfile)) then begin
     splog, filename = logfile
     splog, 'Log file ' + logfile + ' opened ' + systime()
  endif
  splog, 'IDL version: ' + string(!version, format = '(99(a," "))')
  spawn, 'uname -a', uname
  splog, 'UNAME: ' + uname[0]
  
  splog, 'idlutils version ' + idlutils_version()
  splog, 'Longslit version ' + longslit_version()
  
  if (keyword_set(plotfile)) then begin
     thisfile = findfile(plotfile, count = ct)
     IF (ct EQ 0 OR KEYWORD_SET(NUCLEAR)) THEN BEGIN
        splog, 'Plot file ' + plotfile
        dfpsplot, plotfile, /color
     ENDIF ELSE BEGIN
        cpbackup, plotfile
        splog, 'Plot file already exists. Creating backup'
        splog, 'Plot file ' + plotfile
        dfpsplot, plotfile, /color
     ENDELSE
  ENDIF
  
  ;;----------
  ;; Loop over each group
  
  ig = WHERE(planstr.GROUP GE 0, ngd)
  IF ngd EQ 0 THEN RETURN ELSE planstr = planstr[ig]
  group_list = planstr.GROUP
  group_list = group_list[uniq(group_list, sort(group_list))]
  mode_list = planstr[uniq(group_list, sort(group_list))].MODE
  ngroup = n_elements(group_list)
  for igroup = 0L, ngroup-1L do begin
     indx = where(planstr.GROUP EQ group_list[igroup])
     ;; Generate the slitmask file
     qpix = WHERE(planstr[indx].flavor EQ 'dflat-lamp', nqpix)
     IF nqpix GT 0 THEN BEGIN
        ithis = indx[qpix[0]]
        slitfile = 'slits-' + planstr[ithis].filename
        dindx = WHERE(planstr.FILEINDX EQ planstr[ithis].DARKINDX, ndark)
        IF ndark GT 0 THEN fdarkfile = djs_filepath(planstr[dindx].FILENAME, root_dir = indir) $
        ELSE fdarkfile = ''
        thisfile = findfile(slitfile + '*', count = ct)
        if (ct EQ 0 OR keyword_set(clobber_calib)) then begin
           splog, 'Generating slits for GROUP=', group_list[igroup]
           luci_slitmask, djs_filepath(planstr[ithis].filename, root_dir = indir) $
                          , slitfile $
                          , darkfile = fdarkfile $
                          , minslit = minslit $
                          , peakthresh = slitthresh $
                          , y1 = slity1, y2 = slity2 $
                          , ksize = ksize, nfind = nfind $
                          , verbose = verbose, /SOFI
        endif else begin
           splog, 'Do not overwrite existing slitmask file ', $
                  thisfile
        endelse
     endif else begin
        slitfile = ''
        splog, 'No input flats for GROUP=', group_list[igroup], ' MODE =', mode_list[igroup]
     endelse
     ;; Generate the superdark for computing the wavelength images
     qpix_dark = WHERE(planstr[indx].flavor EQ 'dark', nqpix_dark)
     IF nqpix_dark GT 0 THEN BEGIN
        ;; First check that the exposure times match
        sciind = WHERE(WHERE(planstr[indx].FLAVOR EQ 'science'), nsci)
        IF nsci GT 0 THEN BEGIN
           exptimes = planstr[WHERE(planstr[indx].FLAVOR EQ 'science')].EXPTIME
           exptimes = exptimes[uniq(exptimes, sort(exptimes))]
           IF n_elements(exptimes) GT 1 THEN $
              message, 'ERROR: sofi_reduce can only ingest darks with a single exposure time matching the science data'
        ENDIF ELSE $
           message, 'ERROR: cannot construct darks without knowing exptime from science frame'
        qpix_dark = WHERE(planstr[indx].flavor EQ 'dark' AND $
                          planstr[indx].EXPTIME EQ exptimes, nqpix_dark)
        ithis_dark = indx[qpix_dark]
        superdarkfile = 'dark-' + planstr[ithis_dark[0]].filename
        thisfile = findfile(superdarkfile + '*', count = ct)
        if (ct EQ 0 OR keyword_set(clobber_calib)) then begin
           splog, 'Generating superdark for GROUP=', group_list[igroup]
           darkimg = niri_superdark(djs_filepath(planstr[ithis_dark].filename, $
                                                 root_dir = indir) $
                                    , outfile = superdarkfile, /SOFI)
        endif else begin
           splog, 'Do not overwrite existing dark file ', $
                  thisfile
        endelse
     endif else begin
        superdarkfile = ''
        splog, 'No input darks for GROUP=', group_list[igroup], ' MODE =', mode_list[igroup]
     endelse
     ;; Generate the pixel flat and illumination flat
     idom = WHERE(planstr[indx].flavor EQ 'dflat-lamp', ndom)
     itwi = WHERE(planstr[indx].flavor EQ 'twiflat-sky', ntwi)
     ;; Dome flats are required to run the flat making code.
     ;; Twilight flats are not.
     IF ndom GT 0 THEN BEGIN 
        pixflatfile = 'pixflat-' + planstr[indx[idom[0]]].filename
        thispixflatfile = findfile(pixflatfile +'*', count = pixct)
        ;; If twilight flats exist use them, otherwise compute
        ;; illum from domeflat
        IF ntwi GT 0 THEN BEGIN 
           illumflatfile = 'illumflat-' + planstr[indx[itwi[0]]].filename
           thisillumflatfile = findfile(illumflatfile + '*', count = illumct)
        ENDIF ELSE BEGIN
           illumflatfile = 'illumflat-' + planstr[indx[idom[0]]].filename
           thisillumflatfile = findfile(illumflatfile + '*', count = illumct)
        ENDELSE
        if (pixct EQ 0 OR illumct EQ 0 OR keyword_set(clobber_calib)) then begin
           splog, 'Generating pixel and illum flat for GROUP=', group_list[igroup] $
                  , ' MODE =', mode_list[igroup]
           IF ntwi GT 0 THEN BEGIN
              iflat = [idom, itwi]
              use_pixel = [lonarr(ndom) + 1L, lonarr(ntwi)]
              use_illum = [lonarr(ndom), lonarr(ntwi) + 1L]
           ENDIF ELSE BEGIN
              iflat = idom
              use_pixel = lonarr(ndom) + 1L
              use_illum = lonarr(ndom) + 1L
           ENDELSE
           infiles = djs_filepath(planstr[indx[iflat]].filename, root_dir = indir)
           nflat = n_elements(infiles)
           darkfiles = strarr(nflat)
           FOR jj = 0L, nflat-1L DO BEGIN
              dindx = WHERE(planstr.FILEINDX EQ planstr[indx[iflat[jj]]].DARKINDX, ndd)
              IF ndd GT 0 THEN darkfiles[jj] = planstr[dindx].FILENAME $
              ELSE message, 'Could not find the dark file for flat file:' $
                            ,  planstr[indx[iflat[jj]]].FILENAME
           ENDFOR
           darkfiles = djs_filepath(darkfiles, root_dir = indir)
           isci = WHERE(planstr[indx].FLAVOR EQ 'science', nsci)
           IF nsci EQ 0 THEN message, $
              'ERROR: Cannot create flats without a science file for wavelength calibration'
           scifile = djs_filepath(planstr[indx[isci[0]]].FILENAME, root_dir = indir)
           sofi_proc, scifile, arcimg, arcivar, hdr = hdr_arc, darkfile = superdarkfile
           tset_slits = mrdfits(slitfile, 1)
           waveimg = luci_waveimg(arcimg, arcivar, tset_slits, hdr_arc, piximg = piximg $
                                  , CALIB = CALIB)
           splog, 'Generating pixel and illum flat for GROUP=', group_list[igroup]
           long_superflat, infiles, pixflatfile, illumflatfile $
                           , waveimg = waveimg, piximg = piximg  $
                           , slitfile = slitfile $
                           , darkfiles = darkfiles $
                           , use_illum = use_illum, use_pixel = use_pixel $
                           , tempdir = tempdir, slitsamp = 5.0, CHK = CHK, /SOFI
        endif else begin
           splog, 'Do not overwrite existing pixel flat and illum flat files ', $
                  thispixflatfile, thisillumflatfile
        endelse
     endif else begin
        pixflatfile = ''
        illumflatfile = ''
        splog, 'No input flats for GROUP=', group_list[igroup], ' MODE =', mode_list[igroup]
     endelse
     sci_inds = WHERE(planstr[indx].FLAVOR EQ 'science', nsci1)
     tel_inds = WHERE(planstr[indx].FLAVOR EQ 'tell', ntel1)
     IF nsci1 GT 0 THEN BEGIN ;; cannot reduce telluric without science 
        ;; Do Telluric first
        IF ntel1 GT 0 THEN BEGIN
           ;; Find uniq telluric targets
           tel_list = planstr[indx[tel_inds]].TARGET
           tel_list = tel_list[uniq(tel_list, sort(tel_list))]
           ntel = n_elements(tel_list)
           FOR itel = 0L, ntel-1L DO BEGIN
              targdir =  scidir + '/' + strcompress(tel_list[itel], /rem)  $
                         + '_' + strcompress(string(group_list[igroup]), /rem) 
              spawn, '\mkdir -p '+targdir
              ithis = WHERE(planstr[indx[tel_inds]].TARGET EQ tel_list[itel])
              seq_list = planstr[indx[tel_inds[ithis]]].SEQ_INDX
              seq_list = seq_list[uniq(seq_list, sort(seq_list))]
              nseq = n_elements(seq_list)
              FOR iseq = 0L, nseq-1L DO BEGIN 
                 ;; indices of A frame and B frame for this target.
                 aind = where(planstr[indx].TARGET EQ tel_list[itel] AND $
                              planstr[indx].FLAVOR EQ 'tell'  AND $
                              planstr[indx].SEQ_INDX EQ seq_list[iseq] AND $
                              planstr[indx].LOCATION EQ 'A', nA)
                 bind = where(planstr[indx].TARGET EQ tel_list[itel] AND $
                              planstr[indx].FLAVOR EQ 'tell'  AND $
                              planstr[indx].SEQ_INDX EQ seq_list[iseq] AND $
                              planstr[indx].LOCATION EQ 'B', nB)
                 afiles = planstr[indx[aind]].FILENAME
                 bfiles = planstr[indx[bind]].FILENAME
                 ;; Tellurics should always be reduced individaully, i.e. no
                 ;; stacking frames
                 IF nA NE 1 OR nB NE 1 THEN message, 'ERROR: Problem with telluric files.'
                 ;; Wavelength solution from science for telluric
                 wavefile_tell = planstr[WHERE(planstr[indx].FILEINDX EQ $
                                               planstr[indx[aind]].WAVEINDX)].FILENAME 
                 scifile1 = 'tel-' + luci_fileno(afiles) + '-' + luci_fileno(bfiles) + '.fits'
                 waveqafile1 = 'wave-' + luci_fileno(wavefile_tell) + '.ps'
                 scifile = djs_filepath(scifile1, root_dir = targdir)
                 waveqafile = djs_filepath(waveqafile1, root_dir = targdir)
                 thisfile = findfile(scifile, count = ct)
                 ;; Assign the first telluric to the be trace that we
                 ;; will use as the crutch for object tracing. 
                 IF itel EQ 0 AND iseq EQ 0 THEN FILESTD = scifile
                 splog, 'Using sequence ', seq_list[iseq], ' of telluric ', tel_list[itel] $
                        , ' from file:', filestd, ' as the crutch for object tracing'
                 if (ct EQ 0 OR keyword_set(nuclear)) then begin
                    splog, 'Reducing telluric frames ', prelog = scifile
                    luci_reduce_work, djs_filepath(afiles, root_dir = indir) $
                                      , djs_filepath(bfiles, root_dir = indir) $
                                      , slitfile, scifile, waveqafile $
                                      , /TELLURIC $
                                      , WAVEFILE_TELL = djs_filepath(wavefile_tell $
                                                                     , root_dir = indir) $
                                      , pixflatfile = pixflatfile $
                                      , illumflatfile = illumflatfile $
                                      , darkfile = superdarkfile, /SOFI, CALIB = CALIB
                    splog, prelog = ''
                 ENDIF ELSE BEGIN
                    splog, 'Do not overwrite existing telluric frame ' $
                           , scifile
                 ENDELSE
              ENDFOR ;; End loop over telluric frames
           ENDFOR    ;; End loop over tellurics 
        ENDIF        ;; End if for tellurics
        ;; Now do science frames
        sci_list = planstr[indx[sci_inds]].TARGET
        sci_list = sci_list[uniq(sci_list, sort(sci_list))]
        nsci = n_elements(sci_list)
        FOR isci = 0L, nsci-1L DO BEGIN
           targdir =  scidir + '/' + strcompress(sci_list[isci], /rem)  $
                      + '_' + strcompress(string(group_list[igroup]), /rem) 
           spawn, '\mkdir -p '+targdir
           ithis = WHERE(planstr[indx[sci_inds]].TARGET EQ sci_list[isci])
           seq_list = planstr[indx[sci_inds[ithis]]].SEQ_INDX
           seq_list = seq_list[uniq(seq_list, sort(seq_list))]
           nseq = n_elements(seq_list)
           FOR iseq = 0L, nseq-1L DO BEGIN 
              ;; indices of A frame and B frame for this target.
              aind = where(planstr[indx].TARGET EQ sci_list[isci] AND $
                           planstr[indx].FLAVOR EQ 'science'  AND $
                           planstr[indx].SEQ_INDX EQ seq_list[iseq] AND $
                           planstr[indx].LOCATION EQ 'A', nA)
              bind = where(planstr[indx].TARGET EQ sci_list[isci] AND $
                           planstr[indx].FLAVOR EQ 'science'  AND $
                           planstr[indx].SEQ_INDX EQ seq_list[iseq] AND $
                           planstr[indx].LOCATION EQ 'B', nB)
              afiles = planstr[indx[aind]].FILENAME
              bfiles = planstr[indx[bind]].FILENAME
              IF nA EQ 0 OR nB EQ 0 OR (nA NE nB) THEN message, 'ERROR: Problem with science file assignment to A and B sequence.'
              scifile1 = 'sci-' + luci_fileno(afiles[0]) + '-' + $
                         luci_fileno(bfiles[nB-1L]) + '.fits'
              waveqafile1 = 'wave-' + luci_fileno(afiles[0]) + '-' + $
                            luci_fileno(bfiles[nB-1L]) + '.ps'
              scifile = djs_filepath(scifile1, root_dir = targdir)
              waveqafile = djs_filepath(waveqafile1, root_dir = targdir)
              thisfile = findfile(scifile, count = ct)
              if (ct EQ 0 OR keyword_set(nuclear)) then begin
                 splog, 'Reducing telluric frames ', prelog = scifile
                 luci_reduce_work, djs_filepath(afiles, root_dir = indir) $
                                   , djs_filepath(bfiles, root_dir = indir) $
                                   , slitfile, scifile, waveqafile $
                                   , pixflatfile = pixflatfile $
                                   , illumflatfile = illumflatfile $
                                   , darkfile = superdarkfile $
                                   , FILESTD = FILESTD, CALIB = CALIB, /SOFI
                 splog, prelog = ''
              ENDIF ELSE BEGIN
                 splog, 'Do not overwrite existing science frame ' $
                        , scifile
              ENDELSE
           ENDFOR ;; End loop over sequences
        ENDFOR    ;; End loop over science targets 
     ENDIF        ;; End if for having science targets 
  ENDFOR          ;; End loop over groups
  


if (keyword_set(plotfile)) then begin
    dfpsclose
endif

splog, /close

return
end
;------------------------------------------------------------------------------

