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
;   22-Jul-2008  Written by J. F. Hennawi Berkeley
;-  
;-----------------------------------------------------------------------------
PRO gmos_ns_reduce, planfile, clobber = clobber, verbose = verbose $
                    , CHK = CHK, NOZAP = NOZAP, NOTRAPZAP = NOTRAPZAP $
                    , _EXTRA = extra, usage = usage
  if  KEYWORD_SET(USAGE) THEN BEGIN
     print, 'Syntax - ' + $
            'gmos_ns_reduce, planfile, /CLOBBER, /VERBOSE, /NOZAP, NOHELIO' 
     return
  endif 
  
  if (NOT keyword_set(planfile)) then planfile = findfile('plan*.par')
  
  ;;----------
  ;; If multiple plan files exist, then call this script recursively
  ;; for each such plan file.
  
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
        gmos_ns_reduce, planfile[i], clobber = clobber, verbose = verbose $
                        , CHK = CHK, NOZAP = NOZAP, NOTRAPZAP = NOTRAPZAP $
                        , _EXTRA = extra
     ENDFOR
     RETURN
  endif

  ;;----------
  ;; Read the plan file
  planstr = yanny_readone(planfile, hdr = planhdr, /anonymous)
  if (NOT keyword_set(planstr)) then begin
     splog, 'Empty plan file ', planfile
     return
  endif
  
  ;; Truncate the gz stuff
  nfil = n_elements(planstr)
  for qq = 0L, nfil-1 do begin
     slen = strlen(planstr[qq].filename)
     if strmid(planstr[qq].filename, slen-3) EQ '.gz' then $
        planstr[qq].filename = strmid(planstr[qq].filename, 0, slen-3)
  endfor
  logfile = yanny_par(planhdr, 'logfile')
  plotfile = yanny_par(planhdr, 'plotfile')
  indir = yanny_par(planhdr, 'indir')
  scidir  = yanny_par(planhdr, 'scidir')
  edg_pix1 = yanny_par(planhdr, 'edg_pix')
  IF NOT KEYWORD_SET(STD) THEN STD = yanny_par(planhdr, 'std') 
  IF KEYWORD_SET(slity1_1) THEN slity1 = long(slity1_1)
  IF KEYWORD_SET(slity2_1) THEN slity2 = long(slity2_1)
  IF KEYWORD_SET(edg_pix1) THEN edg_pix = long(edg_pix1)
  IF KEYWORD_SET(minslit1) THEN minslit = long(minslit1)

  ;;----------
  ;; Create science dir
  IF keyword_set(scidir) THEN spawn, 'mkdir -p '+scidir

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
  
  ;;----------
  ;; Loop over each unique MASK 
  inmask=where(planstr.FLAVOR NE 'dark' AND planstr.FLAVOR NE 'bias')
  mask_list = planstr[inmask].maskname 
  iuniq = uniq(mask_list, sort(mask_list))
  mask_list  = mask_list[iuniq]
  nmask = n_elements(mask_list)

  ;;-----------------------------
  ;; Make superdark, darkmask
  ;;-----------------------------
  ibias = where(planstr.flavor EQ 'bias', nbias)
  idark = where(planstr.flavor EQ 'dark', ndark)
  if (ndark GT 0) then begin
     superdarkfile = 'superdark-' + planstr[idark[0]].filename
     thisfile = findfile(superdarkfile, count = ct)
     if (ct EQ 0 OR keyword_set(clobber)) then begin
        IF nbias EQ 0 THEN message, 'No inpute biases cannot create superdark'
        superbiasfile = 'superbias-' + planstr[ibias[0]].filename
        splog, 'Generating superdark for Gemini GMOS with N=' $
               , (planstr[idark[0]].ACNT + planstr[idark[0]].BCNT)/2 $
               , ' nod and shuffle cycles'
        gmos_superdark $
           , djs_filepath(planstr[ibias].filename $
                          , root_dir = indir) $
           ,  djs_filepath(planstr[idark].filename $
                           , root_dir = indir) $
           , superbiasfile, superdarkfile, sigrej = sigrej
     endif else begin
        splog, 'Do not overwrite existing superdark ', thisfile
     endelse
  endif else begin
     superdarkfile = ''
     message, 'No input darks in this plan file'
  endelse
  
  
  FOR imask = 0L, nmask-1L DO BEGIN
     ;; Grab all the relevant files for this mask
     jndx = where(planstr.maskname EQ mask_list[imask])
     ;;---------------------------
     ;; Create slitmask file 
     ;;---------------------------
     id = where(planstr[jndx].flavor EQ 'domeflat', nflat)
     im = where(planstr[jndx].flavor EQ 'maskimg', nmask1)
     IF nflat GT 0 AND nmask1 GT 0 THEN BEGIN
        islit = jndx[id[0]] 
        imask1 = jndx[im[0]]
        slitfile = 'slits-' + planstr[islit].filename
        indxfile = 'indx-' + repstr(planstr[imask1].filename, '.fits', '.txt')
        thisfile = findfile(slitfile, count = ct)
        if (ct EQ 0 OR keyword_set(clobber)) then begin
           splog, 'Generating slits for MASK+WAVE=', mask_list[imask] $
                  , planstr[islit].WAVE
           thisfile = findfile(slitfile, count = ct)
           gmos_slitmask $
              , djs_filepath(planstr[islit].filename $
                             , root_dir = indir) $
              ,  djs_filepath(planstr[imask1].filename $
                              , root_dir = indir) $
              , indxfile, slitfile, CHK = CHK
        endif else begin 
           splog, 'Do not overwrite existing slitmask file ', $
                  thisfile
        endelse
     ENDIF ELSE BEGIN
        message, 'No input flats for slitmask for MASK=' + mask_list[imask]
     ENDELSE
     ;;---------------------------
     ;; Make a wavelength solution
     ;;---------------------------
     iarc = where(planstr[jndx].flavor EQ 'arc' OR $
                  planstr[jndx].flavor EQ 'filtarc', narc)
     IF (narc GT 0) then begin
        iarc = jndx[iarc[0]]
        wavefile = 'wave-' + planstr[iarc[0]].filename
        if strpos(wavefile, '.fits') LT 0 then wavefile = wavefile+'.fits'
        thisfile = findfile(wavefile, count = ct)
        IF (ct EQ 0 OR keyword_set(clobber)) then begin
           splog, 'Generating wavelengths for MASK+WAVE=', mask_list[imask] $
                  , planstr[iarc[0]].WAVE
           gmos_wavesolve $
              , djs_filepath(planstr[iarc].filename $
                             , root_dir = indir), slitfile, wavefile 
        ENDIF ELSE BEGIN splog, 'Do not overwrite existing wavelength file ' $
           , wavefile
        ENDELSE
     ENDIF ELSE BEGIN
        wavefile = ''
        message, 'No input arcs for MASK=' +  mask_list[imask]
     ENDELSE
     ;;-----------------------------------
     ;; Finally, reduce each science image
     ;;-----------------------------------
     ii = where(planstr[jndx].flavor EQ 'science', nsci)
     for isci = 0L, nsci-1 do begin
        j = jndx[ii[isci]]
        scifile = djs_filepath('sci-' + planstr[j].filename $
                               , root_dir = scidir)
        if strpos(scifile, '.fits') LT 0 then scifile = scifile+'.fits'
        thisfile = findfile(scifile+'*', count = ct)
        if (ct EQ 0 OR keyword_set(clobber)) THEN BEGIN
           splog, 'Reducing science frame ', prelog = planstr[j].filename
           splog, 'for MASK+WAVE=',  mask_list[imask], planstr[j].WAVE
           gmos_ns_extract $
              , djs_filepath(planstr[j].filename, root_dir = indir) $
              , superdarkfile, slitfile, wavefile, scifile $
              , edg_pix = edg_pix $
              , NOZAP = NOZAP, NOTRAPZAP = NOTRAPZAP
           splog, prelog = ''
        endif else begin
           splog, 'Do not overwrite existing science frame ', scifile
        endelse
     endfor
  endfor                         ; End loop over GRATING+MASK

  splog, /close
   
   return
end
;------------------------------------------------------------------------------
