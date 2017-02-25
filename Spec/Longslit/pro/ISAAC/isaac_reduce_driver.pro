PRO isaac_reduce1, filename, skyfile, scifile $
                  , pixflatfile = pixflatfile, illumflatfile = illumflatfile $
                  , darkfile = darkfile $
                  , TELLURIC = TELLURIC, CHK = CHK, WVCHK = WVCHK $
                  , VERBOSE = VERBOSE, TARGDIR = TARGDIR

  
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


sky_model = isaac_skysub(filename, skyfile, scifile, tset_slits = tset_slits $
                         , pixflatfile = pixflatfile $
                         , illumflatfile = illumflatfile $
                         , sciimg = sciimg, ivar_diff = ivar_diff $
                         , objstruct = objstruct, HDR = SCIHDR $
                         , WAVEIMG = WAVEIMG, SKYIMG = SKYIMG  $
                         , SLITMASK = SLITMASK $
                         , TELLURIC = TELLURIC, CHK = CHK, WVCHK = WVCHK $
                         , VERBOSE = VERBOSE, targdir = targdir $
                         , LOCAL_SKY = LOCAL_SKY)
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
                           , sky_model, objstruct[iobj], plate_scale)
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
       isaac_reduce, planfile[i], clobber = clobber, verbose = verbose
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
    targ_inds = WHERE(planstr[indx].FLAVOR EQ 'SCIENCE' $
                      OR planstr[indx].FLAVOR EQ 'TELL', ntfile)
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
               skyfile = planstr[WHERE(planstr.FILEINDX EQ $
                                       planstr[jndx[jsci]].SKYINDX)].FILENAME
               nsky = n_elements(skyfile)
               IF KEYWORD_SET(TELLURIC) THEN pref = 'tel-' ELSE pref = 'sci-'
               scifile1 = pref + $
                          gnirs_fileprefix(planstr[jndx[jsci]].filename)  $
                          + '-' + gnirs_filepreifx(skyfile) + '.fits'
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
                     telluric = djs_filepath(planstr[indx[objind[kk]]].FILENAME $
                                             , root_dir = indix)
                  ENDIF ELSE TELLURIC = 0
                  isaac_reduce1 $
                     , djs_filepath(planstr[jndx[jsci]].filename $
                                    , root_dir = indir) $
                     , djs_filepath(skyfile $
                                    , root_dir = indir) $
                     , scifile $
                     , pixflatfile = pixflatfile $
                     , illumflatfile = illumflatfile $
                     , darkfile = darkfile $
                     , targdir = targdir $
                     , TELLURIC = TELLURIC, CHK = CHK
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
