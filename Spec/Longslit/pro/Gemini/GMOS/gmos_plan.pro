;+
; NAME:
;   gmos_plan
;
; PURPOSE:
;   Create plan file(s) for running the low-redux pipeline.  This code
;   parses headers, does image stats, etc.
;
; CALLING SEQUENCE:
;   gmos_plan, [ fileexpr, indir, planfile= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   fileexpr   - File names in the input directory; default to '*.fits*'
;   indir      - Input directory(s) for reading files;
;                default to current directory
;   planfile   - Output plan file; default to 'plan.par'.
;                This file is put in the same directory as the raw data files.
;
; OUTPUT:
;
; COMMENTS:
;   One plan file is made for each input directory.
;
;   The following flavors of images are listed:
;     bias
;     domeflat
;     iflat (internal flat)
;     twiflat
;     arc
;     science
;
; EXAMPLES:
; gmos_plan,'*.fits','/b/martell/data_arx/09072005/Raw/',planfile='plan-master.par'
; BUGS:
;
; PROCEDURES CALLED:
;   fileandpath()
;   headfits()
;   idlutils_version()
;   splog
;   sxpar()
;   yanny_write
;
; INTERNAL SUPPORT ROUTINES:
;   gmos_plan_struct()
;
;  To Do: for twilight images taken with a grism in but no mask,
;  assume they are slitless, and treat them as pixel flat field images 
;  not twiflats. 
;
;
;
; REVISION HISTORY:
;   13-Mar-2005  Written by David Schlegel, LBL.
;-
;------------------------------------------------------------------------------
function gmos_plan_struct, nfile
   planstr = create_struct(name='lexp', $
    'FILENAME'    , '', $
    'FLAVOR'      , '', $
    'TARGET'      , '', $
    'EXPTIME'     , 0., $
    'INSTRUMENT'  , '', $
    'GRATING'     , '', $
    'WAVE'        , '', $                       
    'MASKNAME'    , '', $
    'DTAX'        , 0., $
    'ACNT'        , 0L, $
    'BCNT'        , 0L, $
    'AIRMASS'     , 0.)                          
   return, replicate(planstr, nfile)
end
;------------------------------------------------------------------------------
pro gmos_plan, fileexpr, indir, planfile=planfile

   COMMON SITE, lat, lng, tzone
   DRADEG = 180.d0/!dpi
   
   if (NOT keyword_set(indir)) then indir = '.'
   if (NOT keyword_set(fileexpr)) then fileexpr = '*.fits*'
   if (NOT keyword_set(planfile)) then planfile = 'plan.par'

   if not keyword_set(INDIR) then begin
       spawn, '\ls -d '+indir, dirlist
       if (NOT keyword_set(dirlist)) then begin
           splog, 'No input directories found'
           return
       endif
       ndir = n_elements(dirlist)
   endif else begin
       dirlist = [indir]
       ndir = n_elements(dirlist)
   endelse
   ;;----------
   ;; Loop over each input directory
   for idir = 0L, ndir-1L do begin
      splog, 'Working on directory ', dirlist[idir]
      cd, dirlist[idir], current = olddir
      if (idir EQ 0) then origdir = olddir
      filenames = findfile(fileexpr, count = nfile)
      splog, 'Number of FITS files found: ', nfile
      
      if (nfile GT 0) then begin
         planstr = gmos_plan_struct(nfile)
         planstr.filename = repstr(filenames, '.gz', '')
         for i = 0L, nfile-1L do begin
            hdr = xheadfits(filenames[i])
            IF (size(hdr, /tname) EQ 'STRING') then begin
               IF strmatch(sxpar(hdr, 'INSTRUME'), 'GMOS-N*') OR $
               strmatch(sxpar(hdr, 'INSTRUME'), 'GMOS-S*') THEN BEGIN
                  ;;---------------------------
                  ;; This is GMOS on Gemini-N
                  ;;---------------------------
                  exptime = strcompress(sxpar(hdr, 'EXPOSURE'), /rem)
                  obstype = strcompress(sxpar(hdr, 'OBSTYPE'), /rem)
                  object =  strcompress(sxpar(hdr, 'OBJECT'), /rem)
                  instrument = strcompress(sxpar(hdr, 'INSTRUME'), /rem)
                  filter1 = strcompress(sxpar(hdr, 'FILTER1'), /rem)
                  filter2 = strcompress(sxpar(hdr, 'FILTER2'), /rem)
                  grating =  strcompress(sxpar(hdr, 'GRATING'), /rem)
                  IF obstype EQ 'FLAT' THEN BEGIN
                     IF strmatch(grating, '*MIRROR*') THEN $
                        planstr[i].flavor = 'maskimg' $
                     ELSE planstr[i].flavor = 'domeflat' 
                  ENDIF else if obstype EQ 'ARC' THEN BEGIN
                     IF strmatch(filter1, '*open*') THEN $
                        planstr[i].flavor = 'arc' $
                     ELSE planstr[i].flavor = 'filtarc'
                  ENDIF else if obstype EQ 'BIAS'then $
                     planstr[i].flavor = 'bias' $
                  else if (obstype EQ 'OBJECT' AND object EQ 'Twilight') $
                  then planstr[i].flavor = 'twiflat' $
                  else if obstype EQ 'MASK'then planstr[i].flavor $
                     = 'mask' $
                  else if obstype EQ 'DARK'then planstr[i].flavor $
                     = 'dark' $
                  ELSE IF (obstype EQ 'OBJECT' AND object NE 'Twilight') THEN $
                     BEGIN
                     IF strmatch(filter2, '*open*') THEN BEGIN
                        IF exptime LT 120.0 THEN planstr[i].flavor = 'std' $
                        ELSE  planstr[i].flavor = 'science'
                     ENDIF ELSE planstr[i].flavor = 'acq'
                  ENDIF ELSE message,  'Unrecognized obstype'
                  planstr[i].exptime = exptime
                  planstr[i].target = object
                  planstr[i].maskname = strtrim(sxpar(hdr, 'MASKNAME'))
                  instrument = strtrim(sxpar(hdr, 'INSTRUME'))
                  planstr[i].instrument = instrument
                  planstr[i].grating = grating
                   planstr[i].wave = strcompress(sxpar(hdr, 'GRWLEN'), /rem)
                  planstr[i].DTAX = strcompress(sxpar(hdr, 'DTAX'), /rem)
                  planstr[i].AIRMASS = strcompress(sxpar(hdr, 'AIRMASS'), /rem)
                  planstr[i].ACNT = strcompress(sxpar(hdr, 'ANODCNT'), /rem)
                  planstr[i].BCNT = strcompress(sxpar(hdr, 'BNODCNT'), /rem)
               ENDIF ELSE $
                  splog, 'WARNING: Unknown instrument for ', filenames[i]
            ENDIF
        ENDFOR                  ; End loop over files
       
       logfile = repstr(planfile, '.par', '') + '.log'
       plotfile = repstr(planfile, '.par', '') + '.ps'
       hdr = ''
       hdr = [hdr, ' ']
       hdr = [hdr, "logfile '" + logfile + "'   # Log file"]
       hdr = [hdr, "plotfile '" + plotfile + "'   # Plot file"]
       hdr = [hdr, "indir '" + dirlist[idir] + "'   # Raw data directory"]
       hdr = [hdr, "scidir  Science  # Science output directory"]
       hdr = [hdr, "idlutilsVersion '" + idlutils_version() $
              + "'  # Version of idlutils when building plan file"]
       hdr = [hdr, "LongslitVersion '" + longslit_version() $
              + "'  # Version of Longslit when building plan file"]
       ;; Write this plan file
       cd, olddir
       yanny_write, planfile, ptr_new(planstr), hdr = hdr, /align
    endif
endfor                          ; End loop over directories

cd, origdir

return
end
;------------------------------------------------------------------------------
