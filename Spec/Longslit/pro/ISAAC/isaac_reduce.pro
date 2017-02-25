PRO isaac_reduce1, filename, skyfile, scifile $
                  , pixflatfile = pixflatfile, illumflatfile = illumflatfile $
                  , darkfile = darkfile $
                  , TELLURIC = TELLURIC, CHK = CHK, WVCHK = WVCHK $
                  , VERBOSE = VERBOSE, TARGDIR = TARGDIR, FILESTD = FILESTD1 $
                   , CALIB = CALIB

;chk = 1
;wvchk = 1
IF KEYWORD_SET(CHK) THEN set_plot, 'X'
;----------
; Set defaults
if (NOT keyword_set(box_rad)) then box_rad = 8L

t0 = systime(1)
;chk = 0
;------------
;----------
; Perform sky-subtraction 
;;tset_slits = mrdfits(slitfile,1)
;; Replace with traced slits ?????
;tset_slits = isaac_slitset(nx, ny)

IF KEYWORD_SET(FILESTD1) THEN FILESTD = FILESTD1
IF KEYWORD_SET(FILESTD) THEN BEGIN
    splog, 'Using telluric as crutch from ' + filestd
    stdstruct = xmrdfits(filestd, 4, /silent)
    stdmax = max(stdstruct.PEAKFLUX, stdind)
    stdtrace = stdstruct[stdind].XPOS
ENDIF ELSE stdtrace = 0


sky_model = isaac_skysub(filename, skyfile, scifile, tset_slits = tset_slits $
                         , pixflatfile = pixflatfile $
                         , illumflatfile = illumflatfile $
                         , sciimg = sciimg, ivar_diff = ivar_diff $
                         , objstruct = objstruct, HDR = SCIHDR $
                         , WAVEIMG = WAVEIMG, SKYIMG = SKYIMG  $
                         , SLITMASK = SLITMASK $
                         , TELLURIC = TELLURIC, CHK = CHK, WVCHK = WVCHK $
                         , VERBOSE = VERBOSE, targdir = targdir $
                         , LOCAL_SKY = LOCAL_SKY, STDTRACE = STDTRACE $
                         , CALIB = CALIB)
; Read in order set structure and create ordermask
plate_scale = 0.147D
final_struct = 0
;----------
; Loop over objects and extract
IF KEYWORD_SET(objstruct) THEN nobj = n_elements(objstruct) $
ELSE nobj = 0L

FOR iobj = 0L, nobj -1L DO BEGIN
   extract = gnirs_extract(sciimg-sky_model, ivar_diff, waveimg $
                           , (slitmask EQ objstruct[iobj].SLITID) $
                           , sky_model, objstruct[iobj], plate_scale $
                           , SN_GAUSS = 1.0d) 
;; ??? Does SN_GAUSS=2.0 improve things ????
   final_struct = struct_append(final_struct, extract)
ENDFOR

;----------
; Write output file
splog, 'Writing FITS file ', scifile
mwrfits, float(sciimg), scifile, scihdr, /create
mwrfits, float(sky_model), scifile
mwrfits, float(ivar_diff), scifile
mwrfits, float(waveimg), scifile
mwrfits, final_struct, scifile

IF nobj NE 0 THEN niri_plotsci, scifile $
                                , hard_ps = repstr(scifile, '.fits', '.ps') $
                                , box = keyword_set(TELLURIC)

splog, 'Elapsed time = ', systime(1)-t0, ' sec'

;;RETURN
END


;------------------------------------------------------------------------------
PRO isaac_reduce, planfile, clobber = clobber, verbose = verbose $
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
       isaac_reduce, planfile[i], clobber = clobber, verbose = verbose $
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
   sci_inds = WHERE(planstr[indx].FLAVOR EQ 'science', nsci1)
   tel_inds = WHERE(planstr[indx].FLAVOR EQ 'tell', ntel1)
   IF nsci1 GT 0 THEN BEGIN ;; cannot reduce telluric without science 
      ;; Do Telluric first
      IF ntel1 GT 0 THEN BEGIN 
         tel_list = planstr[indx[tel_inds]].TARGET
         tel_list = tel_list[uniq(tel_list, sort(tel_list))]
         ntel = n_elements(tel_list)
         FOR itel = 0L, ntel-1L DO BEGIN
            targdir =  scidir + '/' + strcompress(tel_list[itel], /rem)  $
                       + '_' + strcompress(string(group_list[igroup]), /rem) 
            spawn, '\mkdir -p '+targdir
;           indices of science files for this target
            inow = where(planstr[indx].TARGET EQ tel_list[itel], nimg)
            jndx = indx[inow]
            FOR jimg = 0L, nimg-1L DO BEGIN
               skyfile = planstr[WHERE(planstr.FILEINDX EQ $
                                       planstr[jndx[jimg]].SKYINDX)].FILENAME
               nsky = n_elements(skyfile)
               scifile1 = 'tel-' + $
                          gnirs_fileprefix(planstr[jndx[jimg]].filename)  $
                          + '-' + gnirs_fileprefix(skyfile) + '.fits'
               scifile = djs_filepath(scifile1, root_dir = targdir)
               thisfile = findfile(scifile, count = ct)
               if (ct EQ 0 OR keyword_set(clobber)) then begin
                  splog, 'Reducing telluric frames ', prelog = scifile
                  ;; set telluric to file science file which will 
                  ;;  give wavelengths. We choose the image that best 
                  ;; matches the airmass of telluric sequence
                  airmass_tell = planstr[jndx[jimg]].AIRMASS
;;                  between telluric and science frame???
;;                  objind = WHERE(planstr[indx].FLAVOR EQ 'science')
                  airmass = planstr[indx[sci_inds]].AIRMASS
                  min_diff = min(abs(airmass_tell-airmass), kk)
                  telluric = djs_filepath(planstr[indx[sci_inds[kk]]].FILENAME $
                                          , root_dir = indir)
                  IF KEYWORD_SET(USE_TELL_WAVE) THEN telluric = 'use_tell_wave'
                  isaac_reduce1 $
                     , djs_filepath(planstr[jndx[jimg]].filename $
                                    , root_dir = indir) $
                     , djs_filepath(skyfile $
                                    , root_dir = indir) $
                     , scifile $
                     , pixflatfile = pixflatfile $
                     , illumflatfile = illumflatfile $
                     , darkfile = darkfile $
                     , targdir = targdir $
                     , TELLURIC = TELLURIC $ 
                     , CHK = CHK, CALIB = CALIB, WVCHK = WVCHK 
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
;           indices of science files for this target
         inow = where(planstr[indx].TARGET EQ sci_list[isci], nimg)
         jndx = indx[inow]
         FOR jimg = 0L, nimg-1L DO BEGIN
            skyfile = planstr[WHERE(planstr.FILEINDX EQ $
                                    planstr[jndx[jimg]].SKYINDX)].FILENAME
            nsky = n_elements(skyfile)
            scifile1 = 'sci-' + $
                       gnirs_fileprefix(planstr[jndx[jimg]].filename)  $
                       + '-' + gnirs_fileprefix(skyfile) + '.fits'
            scifile = djs_filepath(scifile1, root_dir = targdir)
            thisfile = findfile(scifile, count = ct)
            if (ct EQ 0 OR keyword_set(clobber)) then begin
               splog, 'Reducing science frames ', prelog = scifile
               ;; If Tellurics exist use them as a crutch for tracing
               IF ntel1 GT 0 THEN BEGIN 
                  tel_targdir =  scidir + '/' + strcompress(tel_list[0], /rem)  $
                                 + '_' + strcompress(string(group_list[igroup]), /rem) 
                  tel_filelist = findfile(tel_targdir + '/tel-*.fits', count = ct)
                  IF ct EQ 0 THEN message, 'Problem with tellurics'
                  filestd = tel_filelist[0] ;; just take the first one
               ENDIF
               isaac_reduce1 $
                  , djs_filepath(planstr[jndx[jimg]].filename $
                                 , root_dir = indir) $
                  , djs_filepath(skyfile $
                                    , root_dir = indir) $
                  , scifile $
                  , pixflatfile = pixflatfile $
                  , illumflatfile = illumflatfile $
                  , darkfile = darkfile $
                  , targdir = targdir $
                  , FILESTD = FILESTD $
                  , CHK = CHK, WVCHK = WVCHK, CALIB = CALIB

               splog, prelog = ''
            ENDIF ELSE BEGIN
               splog, 'Do not overwrite existing science frame ' $
                      , scifile
            ENDELSE
         ENDFOR  ;; End loop over science frames
      ENDFOR     ;; End loop over targets
   ENDIF         ;; End if for targets
ENDFOR           ;; End loop over groups


if (keyword_set(plotfile)) then begin
    dfpsclose
endif

splog, /close

return
end
;------------------------------------------------------------------------------
