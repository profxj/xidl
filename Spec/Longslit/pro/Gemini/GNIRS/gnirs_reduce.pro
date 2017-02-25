; BUGS:
;
;------------------------------------------------------------------------------
; Main reduction routine.
PRO gnirs_reduce1, filenames, scifile, acqfiles, flatfile = flatfile $
                   , ZAP = ZAP, TELLURIC = TELLURIC, CHK = CHK $
                   , VERBOSE = VERBOSE, TARGDIR = TARGDIR, WVCHK = WVCHK $
                   , BOX_RAD = BOX_RAD1
;----------
; Set defaults
IF KEYWORD_SET(BOX_RAD1) THEN BOX_RAD = BOX_RAD1 $
ELSE IF KEYWORD_SET(TELLURIC) THEN BOX_RAD = 15L $
ELSE BOX_RAD = 8L

t0 = systime(1)
;------------
; Read in order set structure and create ordermask
tset_slits = mrdfits(flatfile, 1, silent = (keyword_set(verbose) EQ 0))
;------------
; Account for possible flexure between slits (determined from flat 
; field images) and the data. 
; ???? How necessary is this?? 
IF NOT KEYWORD_SET(NOSHIFT) AND NOT KEYWORD_SET(TELLURIC) THEN BEGIN
    gnirs_proc, filenames[0], sciimg
    xshift = long_xcorr_slits(sciimg, tset_slits, /shift)
ENDIF
plate_scale = 0.15D             ; GNIRS plate scale
dimt = size(tset_slits.coeff, /dimen)
norders = dimt[1]
order_vec = [3, 4, 5, 6, 7, 8]
ordermask = long_slits2mask(tset_slits)
ordermask[WHERE(ordermask GT 0)] = ordermask[WHERE(ordermask GT 0)] + 2L
;----------
; Perform ABBA sky-subtraction 
abba = gnirs_skysub(filenames, flatfile, acqfiles, /OBJMASK $
                    , IVAR_ABBA = IVAR_ABBA, SKY_RESIDS = SKY_RESIDS $
                    , WAVEIMG = WAVEIMG, SKY_MODEL = SKY_MODEL $
                    , OBJ_POS = OBJ_POS, OBJ_NEG = OBJ_NEG $
                    , TELLURIC = TELLURIC, ZAP = ZAP, CHK = CHK $
                    , WVCHK = WVCHK, VERBOSE = VERBOSE, hdr = scihdr $
                    , TARGDIR = TARGDIR)
;;, IVAR_POS=IVAR_POS,IVAR_NEG = IVAR_NEG, 
final_struct = 0
;----------
; Loop over each order and extract
FOR iorder = 0L, norders-1L DO BEGIN
    thismask = (ordermask EQ order_vec[iorder])
    ii_pos = where(obj_pos.SLITID EQ order_vec[iorder], npos)
    ii_neg = where(obj_neg.SLITID EQ order_vec[iorder], nneg)
    if npos NE 1 then begin
        message, 'Error with number of objects npos=', npos, ' on order #' $
                 , order_vec[iorder]
    endif
    if nneg NE 1 then begin
        message, 'Error with number of objects nneg=', nneg, ' on order #' $
                 , order_vec[iorder]
     endif
    extract_pos = gnirs_extract(abba-sky_resids, ivar_abba, waveimg $
                                , thismask, sky_model, obj_pos[ii_pos] $
                                , plate_scale, TELLURIC = TELLURIC)
    extract_neg = gnirs_extract(-abba+sky_resids, ivar_abba, waveimg $
                                , thismask, sky_model, obj_neg[ii_neg] $
                                , plate_scale, TELLURIC = TELLURIC)
    final_struct = struct_append(final_struct, extract_pos)
    final_struct = struct_append(final_struct, extract_neg)
ENDFOR

;----------
; Write output file

splog, 'Writing FITS file ', scifile
mwrfits, float(abba), scifile, scihdr, /create
mwrfits, float(sky_resids), scifile
mwrfits, float(ivar_abba), scifile
mwrfits, float(waveimg), scifile
mwrfits, final_struct, scifile

gnirs_plotsci, scifile, hard_ps = repstr(scifile, '.fits', '.ps') $
               , box = keyword_set(TELLURIC)

splog, 'Elapsed time = ', systime(1)-t0, ' sec'

RETURN
END


;------------------------------------------------------------------------------
PRO gnirs_reduce, planfile, clobber = clobber, verbose = verbose $
                  , chk = chk, wvchk = wvchk

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
      gnirs_reduce, planfile[i], clobber = clobber, verbose = verbose
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
shift8_1 = yanny_par(planhdr, 'shift8')
IF KEYWORD_SET(shift8_1) THEN shift8 = float(shift8_1)
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
    ii = where(planstr[indx].flavor EQ 'FLAT', nflat)
    if (nflat GT 0) then begin
        iflat = indx[ii]
        ;; need one object file. avoid the zeroth since it will 
        ;; have the acquisition residuals. use the first
        objind = WHERE(planstr[indx].FLAVOR EQ 'SCIENCE' AND $
                       planstr[indx].SEQ_NUM EQ 0, nobj)
        IF nobj NE 0 THEN iobj = indx[objind]  $
        ELSE message, 'Need one object file to run gnirs_superflat'
        cd, current = cdir
        superflatfile = cdir + '/superflat-' + $
          gnirs_fileprefix(planstr[iflat[0]].FILENAME) + '-' $
          + strcompress(string(planstr[iflat[nflat-1]].GEMINDX), /REM) + '.fits'
        thisfile = findfile(superflatfile, count = ct)
        if (ct EQ 0 OR keyword_set(clobber)) then begin
            splog, 'Generating superflat for group', group_list[igroup]
            gnirs_superflat, djs_filepath(planstr[iflat].filename $
                                          , root_dir = indir) $
                             ,  djs_filepath(planstr[iobj].filename $
                                             , root_dir = indir) $
                             , superflatfile, /xd32, shift8 = shift8 $
                             , verbose = verbose
        endif else begin
            splog, 'Do not overwrite existing superflat ', thisfile
        endelse
    endif else begin
        superflatfile = ''
        splog, 'No input  for superflat for INSTRUMENT=', $
               group_list[igroup]
    endelse
    
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
;           Incides of science files for this target
            jndx    = indx[where(planstr[indx].TARGET EQ targ_list[itarg] $
                                 AND (planstr[indx].FLAVOR EQ 'SCIENCE' $
                                      OR planstr[indx].FLAVOR EQ 'TELL'))]
            seq_list = planstr[jndx].SEQ_NUM
            seq_list = seq_list[uniq(seq_list, sort(seq_list))]
            nseq = n_elements(seq_list)
            FOR iseq = 0L, nseq-1L DO BEGIN
                inow = WHERE((planstr[jndx].FLAVOR EQ 'SCIENCE' OR $
                              planstr[jndx].FLAVOR EQ 'TELL')  AND $
                             planstr[jndx].SEQ_NUM EQ seq_list[iseq] $
                             , nsci)
;-----------------------------------
;        Reduce each science image
;-----------------------------------
                ii = jndx[inow]
                ii = ii[sort(planstr[ii].ABBA_NUM)]
                IF planstr[ii[0]].FLAVOR EQ 'TELL' THEN pre = 'tel-' ELSE pre = 'sci-'
                scifile1 = pre + gnirs_fileprefix(planstr[ii[0]].filename) $
                  + '-' + strcompress(string(planstr[ii[3]].GEMINDX), /REM)  $
                  + '.fits'
                scifile = djs_filepath(scifile1 $
                                       , root_dir = targdir)
;                profile1 = 'profile-' + $
;                  gnirs_fileprefix(planstr[ii[0]].filename) + '-' $
;                  + strcompress(string(planstr[ii[3]].GEMINDX), /REM) + '.fits'
;                profile_filename = djs_filepath(profile1, root_dir = scidir)
                thisfile = findfile(scifile, count = ct)
                if (ct EQ 0 OR keyword_set(clobber)) then begin
                    splog, 'Reducing science frames ' $
                           , prelog = scifile
                    IF planstr[ii[0]].FLAVOR EQ 'TELL' THEN BEGIN
;                       set telluric to file science file which will 
;                       give wavelengths. We choose the image that best 
;                       matches the airmass of telluric sequence
                        airmass_tell = planstr[ii[0]].AIRMASS
;                       between telluric and science frame???
                        objind = WHERE(planstr[indx].FLAVOR EQ 'SCIENCE')
                        airmass = planstr[indx[objind]].AIRMASS
                        min_diff = min(abs(airmass_tell-airmass), kk)
                        seq_kk = planstr[indx[objind[kk]]].SEQ_NUM
                        group_kk =  planstr[indx[objind[kk]]].GROUP
                        target_kk = planstr[indx[objind[kk]]].TARGET
                        tell_inds = WHERE(planstr[indx].FLAVOR EQ 'SCIENCE' $
                                          AND planstr[indx].SEQ_NUM EQ seq_kk $
                                          AND planstr[indx].GROUP EQ group_kk $
                                          AND planstr[indx].TARGET EQ target_kk)
                        TELLURIC =  $
                           djs_filepath(planstr[indx[tell_inds]].FILENAME $
                                        , root_dir = indir)
                        acq_inds = where(planstr[indx].TARGET EQ $
                                         target_kk $
                                         AND planstr[indx].FLAVOR EQ 'ACQ' $
                                         AND strmatch(planstr[indx].SLIT $
                                                      , '*Acq*'))
                        IF acq_inds[0] EQ -1 THEN message, 'Need to specify acquisition frames'
                        acqfiles = $
                           djs_filepath(planstr[indx[acq_inds]].filename $
                                        , root_dir = indir) 
                     ENDIF ELSE BEGIN
                        acq_inds = where(planstr[indx].TARGET EQ $
                                         targ_list[itarg] $
                                         AND planstr[indx].FLAVOR EQ 'ACQ' $
                                         AND strmatch(planstr[indx].SLIT $
                                                      , '*Acq*'))
                        IF acq_inds[0] EQ -1 THEN message, 'Need to specify acquisition frames'
                        acqfiles = $
                           djs_filepath(planstr[indx[acq_inds]].filename $
                                        , root_dir = indir) 
                        TELLURIC = 0
                     ENDELSE
                    gnirs_reduce1, djs_filepath(planstr[ii].filename $
                                                , root_dir = indir), scifile $
                                   , acqfiles, flatfile = superflatfile $
                                   , TELLURIC = TELLURIC, TARGDIR = TARGDIR $
                                   , CHK = CHK, WVCHK = WVCHK
                    splog, prelog = ''
                ENDIF ELSE BEGIN
                    splog, 'Do not overwrite existing science frame ', scifile
                ENDELSE
            ENDFOR              ; End loop over sequences
        ENDFOR                  ; End loop over targets
    ENDIF                       ; End if for targets
ENDFOR                          ; End loop over groups

if (keyword_set(plotfile)) then begin
    dfpsclose
endif

splog, /close

return
end
;------------------------------------------------------------------------------
