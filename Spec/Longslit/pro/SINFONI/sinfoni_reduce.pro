FUNCTION sinfoni_parse_wavefiles, indx, planstr, waveinds = waveinds

  waveframes = planstr[indx].WAVEFRAMES
  ;; remove brackets
  split1 = strsplit(waveframes, '[*]', /extract)
  indstr = strsplit(split1, ',', /extract)
  waveinds = long(indstr)
  nwave = n_elements(waveinds)
  wavefiles = strarr(nwave)
  FOR j = 0L, nwave-1L DO BEGIN
     jnd = WHERE(planstr.FILEINDX EQ waveinds[j], nmatch)
     IF nmatch NE 1 THEN message, 'Error parsing skyfiles' $
     ELSE wavefiles[j] = planstr[jnd].FILENAME
ENDFOR

RETURN, wavefiles
END

PRO sinfoni_reduce1, objfile, skyfile, wavefiles, slitfile, scifile, waveqafile $
                     , TELLURIC = TELLURIC, CHK = CHK, WVCHK = WVCHK $
                     , pixflatfile = pixflatfile, illumflatfile = illumflatfile $
                     , darkfile = darkfile, FWHM = FWHM1, peakthresh = peakthresh
  
  
  IF KEYWORD_SET(CHK) THEN set_plot, 'X'
  ;;----------
  ;; Set defaults
  if (NOT keyword_set(box_rad)) then box_rad = 8L
  IF KEYWORD_SET(FWHM1) THEN FWHM = FWHM1 ELSE FWHM = 4.5
    
  t0 = systime(1)
    
  ;; Read in slitmask and slit structure 
  slitmask = mrdfits(slitfile, 0)
  tset_slits = mrdfits(slitfile, 1)
  dims = tset_slits[0].DIMS
  nx = dims[0]
  ny = dims[1]
  
  ;; Compute the wavelength image for this source from the science frame
  ;; or skyfile stack
  waveimg = sinfoni_waveimg(wavefiles, tset_slits $
                            , piximg = piximg, CHK = WVCHK, QAFILE = waveqafile $
                            , darkfile = darkfile $
                            , pixflatfile = pixflatfile $
                            , illumflatfile = illumflatfile)
  ;; Perform sky subtraction
  obj_min_sky = sinfoni_skysub(objfile, skyfile, tset_slits, piximg $
                               , ivar = ivar $
                               , sky_resids = sky_resids $
                               , sky_model = sky_model $
                               , hdr_obj = scihdr $
                               , plate_scale = plate_scale $
                               , pixflatfile = pixflatfile $
                               , illumflatfile = illumflatfile, chk = chk)
  ;; Search for object traces in sky-subtracted image
  objstruct = long_objfind(obj_min_sky - sky_resids, tset_slits = tset_slits $
                           , FWHM = FWHM $
                           , NPERSLIT = KEYWORD_SET(TELLURIC) $
                           , OBJTHRESH = KEYWORD_SET(TELLURIC)*0.1 $
                           , peakthresh = peakthresh)

  diff = obj_min_sky - sky_resids
  final_struct = 0
  dimt = size(tset_slits.coeff, /dimen)
  nslit = dimt[1]
  nobj = n_elements(objstruct)
  img_profile = fltarr(nx, ny)
  ;;----------
  ;; Loop over each object and extract
  FOR iobj = 0L, nobj-1L DO BEGIN
     thismask = (slitmask EQ objstruct[iobj].SLITID)
     extract_obj = gnirs_extract(diff, ivar, waveimg $
                                 , thismask, sky_model, objstruct[iobj] $
                                 , plate_scale $
                                 , SN_GAUSS = SN_GAUSS, img_profile = img_profile)
     final_struct = struct_append(final_struct, extract_obj)
  ENDFOR
  ;;----------
  ;; Write output file
  splog, 'Writing FITS file ', scifile
  mwrfits, float(diff), scifile, scihdr, /create
  mwrfits, float(ivar), scifile
  mwrfits, float(waveimg), scifile
  mwrfits, float(img_profile), scifile
  mwrfits, final_struct, scifile
;IF nobj NE 0 THEN niri_plotsci, scifile $
;                                , hard_ps = repstr(scifile, '.fits', '.ps') $
;                                , box = keyword_set(TELLURIC), /SINFONI

splog, 'Elapsed time = ', systime(1)-t0, ' sec'

RETURN
END


;------------------------------------------------------------------------------
PRO sinfoni_reduce, planfile, calib_clobber = calib_clobber, clobber = clobber, verbose = verbose $
                 , CHK = CHK, WVCHK = WVCHK, CALIB = CALIB

if (NOT keyword_set(planfile)) then planfile = findfile('plan*.par')

;----------
; If multiple plan files exist, then call this script recursively
; for each such plan file.

if planfile[0] EQ '' then begin
    print, 'ERROR: Could not find plan file'
    print, 'Try running gnirs_plan'
    return
endif

if (n_elements(planfile) GT 1) then begin
   for i = 0L, n_elements(planfile)-1L do $
      sinfoni_reduce, planfile[i], calib_clobber = calib_clobber, clobber = clobber, verbose = verbose $
                      , CHK = CHK, WVCHK = WVCHK, CALIB = CALIB
   
   return
 endif

;----------
; Read the plan file
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

;----------
; Create science dir
IF keyword_set(scidir) THEN spawn, '\mkdir -p '+scidir

;----------
; Open log file
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
; Loop over each group

ig = WHERE(planstr.GROUP GE 0, ngd)
IF ngd EQ 0 THEN RETURN
group_list = planstr[ig].GROUP
group_list = group_list[uniq(group_list, sort(group_list))]
ngroup = n_elements(group_list)
for igroup = 0L, ngroup-1L do begin
   indx = where(planstr.GROUP EQ group_list[igroup])
   ;; Generate the flatdark file for computing slitmask and flats
   qpix_fdark = WHERE(planstr[indx].flavor EQ 'iflat-dark', nqpix_fdark)
   IF nqpix_fdark GT 0 THEN BEGIN
      ithis_fdark = indx[qpix_fdark]
      fdarkfile = 'flatdark-' + planstr[ithis_fdark[0]].filename
      thisfile = findfile(fdarkfile + '*', count = ct)
      if (ct EQ 0 OR keyword_set(calib_clobber)) then begin
         splog, 'Generating flat dark for GROUP=', group_list[igroup]
         darkimg = niri_superdark(djs_filepath(planstr[ithis_fdark].filename, $
                                               root_dir = indir), outfile = fdarkfile, /SINFONI)
      endif else begin
         splog, 'Do not overwrite existing flatdark file ', $
                thisfile
      endelse
   endif else begin
      slitfile = ''
      splog, 'No input darks for GROUP=', group_list[igroup]
   endelse
   ;; Generate the slitmask file
   qpix = WHERE(planstr[indx].flavor EQ 'iflat-lamp', nqpix)
   IF nqpix GT 0 THEN BEGIN
      ithis = indx[qpix[0]]
      slitfile = 'slits-' + planstr[ithis].filename
      thisfile = findfile(slitfile + '*', count = ct)
      if (ct EQ 0 OR keyword_set(calib_clobber)) then begin
         splog, 'Generating slits for GROUP=', group_list[igroup]
         sinfoni_slitmask, djs_filepath(planstr[ithis].filename, $
                                        root_dir = indir) $
                           , slitfile $
                           , darkfile = fdarkfile $
                           , minslit = minslit $
                           , peakthresh = slitthresh $
                           , y1 = slity1, y2 = slity2 $
                           , ksize = ksize, nfind = nfind $
                           , verbose = verbose
      endif else begin
         splog, 'Do not overwrite existing slitmask file ', $
                thisfile
      endelse
   endif else begin
      slitfile = ''
      splog, 'No input flats for GROUP=', group_list[igroup]
   endelse
   ;; Generate the superdark for computing the wavelengths 
   qpix_dark = WHERE(planstr[indx].flavor EQ 'dark', nqpix_dark)
   IF nqpix_dark GT 0 THEN BEGIN
      ithis_dark = indx[qpix_dark]
      darkfile = 'dark-' + planstr[ithis_dark[0]].filename
      thisfile = findfile(darkfile + '*', count = ct)
      if (ct EQ 0 OR keyword_set(calib_clobber)) then begin
         splog, 'Generating superdark for GROUP=', group_list[igroup]
         darkimg = niri_superdark(djs_filepath(planstr[ithis_dark].filename, $
                                               root_dir = indir), outfile = darkfile, /SINFONI)
      endif else begin
         splog, 'Do not overwrite existing dark file ', $
                thisfile
      endelse
   endif else begin
      slitfile = ''
      splog, 'No input darks for GROUP=', group_list[igroup]
   endelse
   ;; Generate the pixel flat and illumination flat
   qpix = WHERE(planstr[indx].flavor EQ 'iflat-lamp', nqpix)
   IF nqpix GT 0 THEN BEGIN
      ithis = indx[qpix] 
      pixflatfile = 'pixflat-' + planstr[ithis[0]].filename
      illumflatfile = 'illumflat-' + planstr[ithis[0]].filename
      thispix = findfile(pixflatfile + '*', count = ct_pix)
      thisillum = findfile(illumflatfile + '*', count = ct_ill)
      if (ct_pix EQ 0 OR ct_ill EQ 0 OR keyword_set(calib_clobber)) then begin
         CASE planstr[ithis[0]].MODE OF
            ;; Just use all the sky frames available
            '0.025': BEGIN isky_flat = WHERE(planstr[indx].FLAVOR EQ 'sky', nsky_flat)
               IF nsky_flat EQ 0 THEN message, 'Sky frames are needed for wavelenghts to construct flats'
            END
            '0.25': BEGIN
               isky_flat1 =  WHERE(planstr[indx].FLAVOR EQ 'science', nsky_flat1)
               IF nsky_flat1 EQ 0 THEN message, 'Science frames are needed for wavelenghts to construct flats'
               isky_flat = isky_flat1[0] ;; pix the first science frame for the flat wavelength soln.
            END
            ELSE: message, 'Mode not supported'
         ENDCASE
         splog, 'Generating pixel and illum flat for GROUP=', group_list[igroup]
         sinfoni_superflat, djs_filepath(planstr[ithis].filename, $
                                         root_dir = indir) $
                            , djs_filepath(planstr[indx[isky_flat]].filename, $
                                           root_dir = indir) $
                            , fdarkfile, slitfile, pixflatfile, illumflatfile $
                            , darkfile = darkfile $
                            , tempdir = tempdir, verbose = verbose, chk = chk
      endif else begin
         splog, 'Do not overwrite existing pixel flat and illum flat files ', $
                thispix, thisillum
      endelse
   endif else begin
      pixflatfile = ''
      illumflatfile = ''
      splog, 'No input flats for GROUP=', group_list[igroup]
   endelse
   targ_inds = WHERE(planstr[indx].FLAVOR EQ 'science' $
                     OR planstr[indx].FLAVOR EQ 'tell', ntfile)
   IF ntfile NE 0 THEN BEGIN 
      targ_list = planstr[indx[targ_inds]].TARGET 
      targ_list = targ_list[uniq(targ_list, sort(targ_list))]
      ntarg = n_elements(targ_list)
      ;; Loop over targets and reduce all files from eqch sequentially
      FOR itarg = 0L, ntarg-1L DO BEGIN
         targdir =  scidir + '/' + strcompress(targ_list[itarg], /rem)  $
                    + '_' + strcompress(string(group_list[igroup]), /rem)
         spawn, '\mkdir -p '+targdir
         ;; indices of science files for this target
         sci_inds = indx[where(planstr[indx].TARGET EQ targ_list[itarg] $
                               AND (planstr[indx].FLAVOR EQ 'science' $
                                    OR planstr[indx].FLAVOR EQ 'tell'), nsci)]
         FOR isci = 0L, nsci-1L DO BEGIN
            skyindx = where(planstr.FILEINDX EQ planstr[sci_inds[isci]].SKYINDX, nsky)
            IF nsky GT 0 THEN $
               skyfile = planstr[skyindx].FILENAME $
            ELSE message, "Cannot find skyfile for this object"
            IF planstr[sci_inds[isci]].FLAVOR EQ 'tell' THEN TELLURIC = 1 ELSE TELLURIC = 0
            IF TELLURIC THEN sci_pref = 'tel-'  ELSE sci_pref = 'sci-'
            scifile1 = sci_pref + $
                       gnirs_fileprefix(planstr[sci_inds[isci]].filename)  $
                       + '-' + gnirs_fileprefix(skyfile) + '.fits'
            scifile = djs_filepath(scifile1, root_dir = targdir)
            
            wavefiles = sinfoni_parse_wavefiles(sci_inds[isci], planstr, waveinds = waveinds)
            IF NOT TELLURIC THEN BEGIN 
               waveqafile1 = 'wave-' + gnirs_fileprefix(wavefiles[0]) + '.ps'
               waveqafile = djs_filepath(waveqafile1, root_dir = targdir)
            ENDIF ELSE waveqafile = 0
            thisfile = findfile(scifile, count = ct)
            if (ct EQ 0 OR keyword_set(clobber)) then begin
               sinfoni_reduce1, djs_filepath(planstr[sci_inds[isci]].filename, root_dir = indir) $
                                , djs_filepath(skyfile, root_dir = indir) $
                                , djs_filepath(wavefiles, root_dir = indir) $
                                , slitfile, scifile, waveqafile $
                                , TELLURIC = (planstr[sci_inds[isci]].flavor EQ 'tell') $
                                , CHK = CHK, WVCHK = WVCHK $
                                , pixflatfile = pixflatfile, illumflatfile = illumflatfile $
                                , darkfile = darkfile, FWHM = FWHM1, peakthresh = peakthresh
            ENDIF ELSE BEGIN
               splog, 'Do not overwrite existing science frame ' $
                      , scifile
            ENDELSE
         ENDFOR ;; end of loop over science frames
      ENDFOR    ;; end of loop over targets
   ENDIF        ;; endif for target frames existence
ENDFOR          ;; end of loop over groups


if (keyword_set(plotfile)) then begin
    dfpsclose
endif

splog, /close

return
end
;------------------------------------------------------------------------------
