; BUGS:
;   Heliocentric corrections!!!???
;
;------------------------------------------------------------------------------
; The code that follows is the science frame reduction code,
; which should get split out to another routine.
FUNCTION niri_parse_skyfiles, indx, planstr, skyinds = skyinds

skyframes = planstr[indx].SKYFRAMES
; remove brackets
split1 = strsplit(skyframes, '[*]', /extract)
indstr = strsplit(split1, ',', /extract)
skyinds = long(indstr)
nsky = n_elements(skyinds)
skyfiles = strarr(nsky)
FOR j = 0L, nsky-1L DO BEGIN
    jnd = WHERE(planstr.GEMINDX EQ skyinds[j], nmatch)
    IF nmatch NE 1 THEN message, 'Error parsing skyfiles' $
    ELSE skyfiles[j] = planstr[jnd].FILENAME
ENDFOR

RETURN, skyfiles
END

PRO niri_reduce1, filename, skyfiles, scifile, flatfile = flatfile $
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
sky_model = niri_skysub(filename, skyfiles, flatfile $
                        , sciimg = sciimg, ivar = ivar $
                        , objstruct = objstruct, HDR = SCIHDR $
                        , WAVEIMG = WAVEIMG, AVG_SKY = AVG_SKY  $
                        , TELLURIC = TELLURIC, CHK = CHK, WVCHK = WVCHK $
                        , VERBOSE = VERBOSE, targdir = targdir)
; Read in order set structure and create ordermask
dims = size(sciimg, /dim)
nx = dims[0]
ny = dims[1]
tset_slits = niri_slitset(nx, ny)
slitmask = long_slits2mask(tset_slits)
plate_scale = 0.1171D

final_struct = 0
;----------
; Loop over objects and extract
IF KEYWORD_SET(objstruct) THEN nobj = n_elements(objstruct) $
ELSE nobj = 0L

FOR iobj = 0L, nobj -1L DO BEGIN
    extract = gnirs_extract(sciimg-sky_model, ivar, waveimg, slitmask $
                            , sky_model, objstruct[iobj], plate_scale)
    final_struct = struct_append(final_struct, extract)
ENDFOR

;----------
; Write output file
splog, 'Writing FITS file ', scifile
mwrfits, float(sciimg), scifile, scihdr, /create
mwrfits, float(sky_model), scifile
mwrfits, float(ivar), scifile
mwrfits, float(waveimg), scifile
mwrfits, final_struct, scifile

IF nobj NE 0 THEN $
  niri_plotsci, scifile, hard_ps = repstr(scifile, '.fits', '.ps') $
  , box = keyword_set(TELLURIC)

splog, 'Elapsed time = ', systime(1)-t0, ' sec'

RETURN
END


;------------------------------------------------------------------------------
PRO niri_reduce, planfile, clobber = clobber, verbose = verbose $
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
       niri_reduce, planfile[i], clobber = clobber, verbose = verbose
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
    
    ii = where(planstr[indx].flavor EQ 'FLAT', nflat)
    jj = wherE(planstr[indx].flavor EQ 'DARK', ndark)
    if (nflat GT 0) then begin
        iflat = indx[ii]
        idark = indx[jj]
;       need one object file
        objind = WHERE(planstr[indx].FLAVOR EQ 'SCIENCE', nobj)
        IF nobj NE 0 THEN iobj = indx[objind[0]]  $
        ELSE message, 'Need one object file to run niri_superflat'
        cd, current = cdir
        superflatfile = cdir + '/superflat-' + $
          gnirs_fileprefix(planstr[idark[0]].FILENAME) + '-' $
          + strcompress(string(planstr[iflat[nflat-1]].GEMINDX), /REM) + '.fits'
        thisfile = findfile(superflatfile, count = ct)
        if (ct EQ 0 OR keyword_set(clobber)) then begin
            splog, 'Generating superflat for group', group_list[igroup]
            niri_superflat, djs_filepath(planstr[iflat].filename $
                                         , root_dir = indir) $
                            , djs_filepath(planstr[idark].filename $
                                         , root_dir = indir) $
                            , superflatfile $
                            , OBJFILE =  djs_filepath(planstr[iobj].filename $
                                                      , root_dir = indir) $
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
;           indices of science files for this target
            jndx = indx[where(planstr[indx].TARGET EQ targ_list[itarg] $
                              AND (planstr[indx].FLAVOR EQ 'SCIENCE' $
                                   OR planstr[indx].FLAVOR EQ 'TELL'))]
            
            seq_list = planstr[jndx].SEQ_NUM
            seq_list = seq_list[uniq(seq_list, sort(seq_list))]
            nseq = n_elements(seq_list)
;           Pick out objects
            FOR iseq = 0L, nseq-1L DO BEGIN
                inow = WHERE((planstr[jndx].FLAVOR EQ 'SCIENCE' OR $
                              planstr[jndx].FLAVOR EQ 'TELL')  AND $
                             planstr[jndx].SEQ_NUM EQ seq_list[iseq], nsci)
                ii = jndx[inow] ; indices for files in this sequence
                FOR jsci = 0L, nsci-1L DO BEGIN
                    skyfiles = niri_parse_skyfiles(jsci, planstr[ii] $
                                                   , skyinds = skyinds)
                    nsky = n_elements(skyfiles)
                    IF planstr[ii[jsci]].FLAVOR EQ 'TELL' THEN scipref = 'tel-' ELSE scipref = 'sci-'
                    scifile1 = scipref + $
                      gnirs_fileprefix(planstr[ii[jsci]].filename)  $
                      + '-' + strcompress(string(skyinds[0]), /REM) $
                      + '-' + strcompress(string(skyinds[nsky-1L]), /REM) $
                      + '.fits'
                    scifile = djs_filepath(scifile1, root_dir = targdir)
                    thisfile = findfile(scifile, count = ct)
                    if (ct EQ 0 OR keyword_set(clobber)) then begin
                        splog, 'Reducing science frames ', prelog = scifile
                        IF planstr[ii[jsci]].FLAVOR EQ 'TELL' THEN BEGIN
;                           set telluric to file science file which will 
;                           give wavelengths. We choose the image that best 
;                           matches the airmass of telluric sequence
                            airmass_tell = planstr[ii[jsci]].AIRMASS
;                           between telluric and science frame???
                            objind = WHERE(planstr[indx].FLAVOR EQ 'SCIENCE')
                            airmass = planstr[indx[objind]].AIRMASS
                            min_diff = min(abs(airmass_tell-airmass), kk)
                            sky_tell = $
                              niri_parse_skyfiles(kk, planstr[indx[objind]] $
                                                  , skyinds = skyinds)
                            telluric = djs_filepath(sky_tell $
                                                    , root_dir = indir) 
;                              djs_filepath(planstr[indx[objind[kk]]].FILENAME $
;                                           , root_dir = indir)
                        ENDIF ELSE TELLURIC = 0
                        niri_reduce1, djs_filepath(planstr[ii[jsci]].filename $
                                                   , root_dir = indir) $
                                      , djs_filepath(skyfiles $
                                                     , root_dir = indir) $
                                      , scifile, flatfile = superflatfile $
                                      , targdir = targdir $
                                      , TELLURIC = TELLURIC, CHK = CHK
                        splog, prelog = ''
                    ENDIF ELSE BEGIN
                        splog, 'Do not overwrite existing science frame ' $
                               , scifile
                    ENDELSE
                ENDFOR          ; End loop over science frames
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
