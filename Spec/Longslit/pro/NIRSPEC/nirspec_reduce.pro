; BUGS:
;   Heliocentric corrections!!!???
;
;------------------------------------------------------------------------------
; The code that follows is the science frame reduction code,
; which should get split out to another routine.
FUNCTION nirspec_parse_skyfiles, indx, planstr, skyinds = skyinds
skyframes = planstr[indx].SKYFRAMES
; remove brackets
split1 = strsplit(skyframes, '[*]', /extract)
indstr = strsplit(split1, ',', /extract)
skyinds = long(indstr)
nsky = n_elements(skyinds)
skyfiles = strarr(nsky)
FOR j = 0L, nsky-1L DO BEGIN
   jnd = WHERE(planstr.KECKINDX EQ skyinds[j], nmatch)
   IF nmatch NE 1 THEN message, 'Error parsing skyfiles' $
   ELSE skyfiles[j] = planstr[jnd].FILENAME
ENDFOR

RETURN, skyfiles
END


PRO nirspec_reduce1, filename, skyfile, scifile, slitfile $
                     , pixflatfile = pixflatfile,illumflatfile=illumflatfile $
                     , TELLURIC = TELLURIC, CHK = CHK, WVCHK = WVCHK $
                     , VERBOSE = VERBOSE, TARGDIR = TARGDIR, MXSHFT=mxshft $
                     , simple_sub=simple_sub, MAXOBJ=maxobj, planstr=planstr
  
IF KEYWORD_SET(CHK) THEN set_plot, 'X'
;----------
; Set defaults
if (NOT keyword_set(box_rad)) then box_rad = 8L
if (NOT keyword_set(MAXOBJ)) then MAXOBJ = 10L

t0 = systime(1)
;chk = 0
;------------
;----------
; Perform sky-subtraction 
tset_slits = mrdfits(slitfile,1)
sky_model = nirspec_skysub(filename, skyfile, tset_slits $
                           , pixflatfile=pixflatfile $
                           , illumflatfile=illumflatfile $
                           , simple_sub=simple_sub $
                           , sciimg = sciimg, ivar_diff = ivar_diff $
                           , objstruct = objstruct, HDR = SCIHDR $
                           , WAVEIMG = WAVEIMG, SKYIMG=SKYIMG  $
                           , SLITMASK = SLITMASK, MXSHFT=mxshft $
                           , TELLURIC = TELLURIC, CHK = CHK, WVCHK = WVCHK $
                           , FILTER=planstr.filter $
                           , VERBOSE = VERBOSE, targdir = targdir)
; Read in order set structure and create ordermask
plate_scale = 0.143D
final_struct = 0
;----------
; Loop over objects and extract
IF KEYWORD_SET(objstruct) THEN nobj = n_elements(objstruct) $
ELSE nobj = 0L

nobj = nobj < MAXOBJ

FOR iobj = 0L, nobj -1L DO BEGIN
   ;; GNIRS is now expecting microns, not Ang
   ;stop
   extract = gnirs_extract(sciimg-sky_model, ivar_diff, waveimg/1e4, slitmask $
                           , sky_model, objstruct[iobj], plate_scale)
   ;; Back to Ang
   extract.wave_opt = extract.wave_opt * 1e4 
   extract.wave_box = extract.wave_box * 1e4 
   ;; Save
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

RETURN
END


;------------------------------------------------------------------------------
;; MXSHIFT -- Keyword to allow for larger shifts for non-archived
;;            spectra
PRO nirspec_reduce, planfile, clobber = clobber, verbose = verbose $
                    , CHK = CHK, WVCHK = WVCHK

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
       nirspec_reduce, planfile[i], clobber = clobber, verbose = verbose
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
mxshft  = yanny_par(planhdr, 'mxshft')
simple_sub  = yanny_par(planhdr, 'simple_sub')
maxobj  = yanny_par(planhdr, 'maxobj')

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

group_list = planstr.GROUP
group_list = group_list[uniq(group_list, sort(group_list))]
ngroup = n_elements(group_list)

for igroup = 0L, ngroup-1L do begin
    indx = where(planstr.GROUP EQ group_list[igroup])
    
    qboth = planstr[indx].flavor EQ 'bothflat'
    qtwi = (planstr[indx].flavor EQ 'twiflat') OR qboth
    qpix = (planstr[indx].flavor EQ 'iflat') OR qboth
    itwi  = where(qtwi, ntwi)
    ipix  = where(qpix, npix)
    iboth = WHERE(qtwi OR qpix, nboth)
    idark = WHERE(planstr[indx].flavor EQ 'dark', ndark)
    targ_inds = WHERE(planstr[indx].FLAVOR EQ 'SCIENCE' $
                      OR planstr[indx].FLAVOR EQ 'TELL', ntfile)
    if (npix GT 0 OR ntwi GT 0) then begin
       ;; Use twiflats if possible, otherwise domeflats.
       ;; (Use only the first such file).
       if (ntwi GT 0) then ithis = indx[itwi[0]] $
       else ithis = indx[ipix[0]] 
       id=indx[idark[0]]
       slitfile = 'slits-' + planstr[ithis].filename
       thisfile = findfile(slitfile, count = ct)
       if (ct EQ 0 OR keyword_set(clobber)) then begin
          splog, 'Generating slits for GROUP=', group_list[igroup]
          nirspec_traceorders,djs_filepath(planstr[ithis].filename, $
                                           root_dir = indir) $
                              ,djs_filepath(planstr[id].filename, $
                                            root_dir = indir) $
                              , slitfile, peakthresh = peakresh $
                              , y1 = slity1, y2 = slity2 $
                              , ksize = ksize, radius=radius $
                              , nave = nave, maxshifte = maxshifte $
                              , maxshift0 = maxshift0, func = func $
                              , ncoeff = ncoeff, CRUDE = CRUDE
       endif else begin
          splog, 'Do not overwrite existing slitmask file ', $
                 thisfile
       endelse
    endif else begin
       slitfile = ''
       splog, 'No input flats for slitmask for GROUP=', group_list[igroup]
    endelse
    
    ;;------------------------------------
    ;; Make  pixel and illumination flats 
    ;;------------------------------------
    IF (npix GT 0) OR (ntwi GT 0) THEN BEGIN 
       IF npix GT 0 THEN BEGIN
          pixflatfile = 'pixflat-' + planstr[indx[ipix[0]]].filename
          thispixflatfile = findfile(pixflatfile, count = pixct)
       ENDIF ELSE pixct = 0
       IF ntwi GT 0 THEN BEGIN
          illumflatfile = 'illumflat-' + planstr[indx[itwi[0]]].filename
          thisillumflatfile = findfile(illumflatfile, count = illumct)
       ENDIF ELSE illumct = 0
       if (pixct EQ 0 OR illumct EQ 0 OR keyword_set(clobber)) then begin
          splog, 'Generating pixel flat for GROUP=', group_list[igroup]
          IF illumct EQ 0 AND pixct EQ 0 THEN BEGIN
             infiles =  djs_filepath(planstr[indx[iboth]].filename $
                                     , root_dir = indir)
             use_illum = qtwi[iboth] 
             use_pixel = qpix[iboth] 
          ENDIF ELSE IF illumct EQ 1 AND pixct EQ 0 THEN BEGIN
             infiles =  djs_filepath(planstr[indx[ipix]].filename $
                                     , root_dir = indir)
             use_illum = lonarr(npix)
             use_pixel = lonarr(npix) + 1L
             splog, 'Do not overwrite existing illumination flat ' $
                    , thisillumflatfile
          ENDIF ELSE IF illumct EQ 0 AND pixct EQ 1 THEN BEGIN
             if ntwi NE 0 then begin 
                infiles =  djs_filepath(planstr[indx[itwi]].filename $
                                        , root_dir = indir) 
                use_illum = lonarr(ntwi) + 1L
                use_pixel = lonarr(ntwi)
             endif else begin
                infiles =  $
                   djs_filepath(planstr[indx[ipix[0]]].filename $
                                , root_dir = indir) 
                use_illum = 0L
                use_pixel = 0L
             endelse
             splog, 'Do not overwrite existing pixel flats ' $
                    , thispixflatfile
          ENDIF
          ;; Dark files for making superdark
          darkfiles=djs_filepath(planstr[indx[idark]].filename $
                                 , root_dir = indir) 
          superdarkfile = 'dark-' + planstr[indx[idark[0]]].filename
          ;; need one object file for wavelength solution
          objind = WHERE(planstr[indx].FLAVOR EQ 'science', nobj)
          IF nobj NE 0 THEN iobj = indx[objind[0]]  $
          ELSE message, 'Need one object file to run nirspec_superflat'
          objfile =  djs_filepath(planstr[iobj].filename, root_dir = indir)
          nirspec_superflat, infiles, darkfiles $
                             ,superdarkfile,pixflatfile, illumflatfile $
                             , slitfile = slitfile $
                             , filter=planstr[iobj].filter, mxshft=mxshft $
                             , objfile = objfile $
                             , verbose = verbose, indir = indir $
                             , tempdir = tempdir $
                             , use_illum = use_illum $
                             , use_pixel = use_pixel $
                             , npoly = npoly, CHK = CHK $
                             , _EXTRA = extra 
          
       ENDIF ELSE BEGIN
          splog, 'Do not overwrite existing pixel flats ' $
                 , thispixflatfile
          splog, 'Do not overwrite existing illumination flat ' $
                 , thisillumflatfile
       ENDELSE
    endif else begin
       pixflatfile = ''
       illumflatfile=''
       splog, 'No input pixel flats or illum flats for GROUP=', $
              group_list[igroup]
    endelse
    
    targ_inds = WHERE(planstr[indx].FLAVOR EQ 'science' $
                      OR planstr[indx].FLAVOR EQ 'tell', ntfile)
    IF ntfile NE 0 THEN BEGIN
        targ_list = planstr[indx[targ_inds]].TARGET 
        targ_list = targ_list[uniq(targ_list, sort(targ_list))]
        ntarg = n_elements(targ_list)
        
        FOR itarg = 0L, ntarg-1L DO BEGIN
            targdir =  scidir + '/' + strcompress(targ_list[itarg], /rem)  $
              + '_' + strcompress(string(group_list[igroup]), /rem) 
            spawn, '\mkdir -p '+targdir
;           indices of science files for this target
            inow = where(planstr[indx].TARGET EQ targ_list[itarg] $
                         AND (planstr[indx].FLAVOR EQ 'science' $
                              OR planstr[indx].FLAVOR EQ 'tell'),nsci)
            jndx = indx[inow]
            FOR jsci = 0L, nsci-1L DO BEGIN
               skyfile = nirspec_parse_skyfiles(jsci, planstr[jndx] $
                                                , skyinds = skyinds)
               nsky = n_elements(skyfile)
               scifile1 = 'sci-' + $
                          gnirs_fileprefix(planstr[jndx[jsci]].filename)  $
                          + '-' + strcompress(string(skyinds[0]), /REM) $
                          + '.fits'
               scifile = djs_filepath(scifile1, root_dir = targdir)
               thisfile = findfile(scifile, count = ct)
               if (ct EQ 0 OR keyword_set(clobber)) then begin
                  splog, 'Reducing science frames ', prelog = scifile
                  IF planstr[jndx[jsci]].FLAVOR EQ 'tell' THEN BEGIN
;                           set telluric to file science file which will 
;                           give wavelengths. We choose the image that best 
;                           matches the airmass of telluric sequence
                     airmass_tell = planstr[jndx[jsci]].AIRMASS
;                           between telluric and science frame???
                     objind = WHERE(planstr[indx].FLAVOR EQ 'science')
                     airmass = planstr[indx[objind]].AIRMASS
                     min_diff = min(abs(airmass_tell-airmass), kk)
                     sky_tell = $
                        nirspec_parse_skyfiles(kk, planstr[indx[objind]] $
                                               , skyinds = skyinds)
                     telluric = djs_filepath(sky_tell $
                                             , root_dir = indir) 
;                              djs_filepath(planstr[indx[objind[kk]]].FILENAME $
;                                           , root_dir = indir)
                  ENDIF ELSE TELLURIC = 0
                  nirspec_reduce1 $
                     , djs_filepath(planstr[jndx[jsci]].filename $
                                    , root_dir = indir) $
                     , djs_filepath(skyfile $
                                    , root_dir = indir) $
                     , scifile, slitfile $
                     , pixflatfile = pixflatfile $
                     , simple_sub = simple_sub, MAXOBJ=maxobj $
                     , illumflatfile=illumflatfile $
                     , targdir = targdir, MXSHFT=mxshft $
                     , TELLURIC = TELLURIC, CHK = CHK $
                     , planstr=planstr[jndx[jsci]]
                  splog, prelog = ''
               ENDIF ELSE BEGIN
                  splog, 'Do not overwrite existing science frame ' $
                         , scifile
               ENDELSE
            ENDFOR              ; End loop over science frames
         ENDFOR                 ; End loop over targets
     ENDIF                      ; End if for targets
 ENDFOR                         ; End loop over groups

if (keyword_set(plotfile)) then begin
    dfpsclose
endif

splog, /close

return
end
;------------------------------------------------------------------------------
