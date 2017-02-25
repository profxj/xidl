PRO GMOS_COADD, infiles, slitno, outfil, CHK = CHK, DEBUG = DEBUG $
                , sigrej = sigrej, medscale = medscale $
                , LAM_MASK_RANGE = LAM_MASK_RANGE
  
;objno = [32] ;; The object numbers to be coadded
;; the relevant science files
;infiles  = ['sci-N20080407S0007.fits'] ;['sci-N20080407S0007.fits', 

nfiles = n_elements(infiles)
;outfil = 'ARC1.fits'

head0 = xheadfits(infiles[0])
;; Compile one master object structure containing everything
FOR ifile = 0L, nfiles-1L DO BEGIN
   objnow = mrdfits(infiles[ifile], 4, hdr)
   IF NOT KEYWORD_SET(OBJSTRUCT) THEN objstruct = objnow $
   ELSE objstruct = [objstruct, objnow]
ENDFOR
;; Find the indices which contain the objects of interest
nstruct = n_elements(objstruct)
nobj = n_elements(slitno)
FOR iobj = 0L, nobj-1L DO BEGIN
   ind = where(objstruct.SLIT EQ slitno[iobj], nthis)
   IF nthis EQ 0 THEN message, 'Problem with input slitno. Could not find' $
   ELSE BEGIN
      IF NOT KEYWORD_SET(ALLIND) THEN ALLIND = ind $
      ELSE ALLIND = [allind, ind]
   ENDELSE
ENDFOR
exptime = total(objstruct[allind].EXPTIME)
nimgs = n_elements(allind)
splog, 'Searched for spectra of slitno #', slitno
splog, 'Found ' + strcompress(string(nimgs, format = '(I3)'), /rem) $
       + ' for a total exptime = ' + $
       strcompress(string(exptime, format = '(F7.1)'), /rem) + ' sec'
influx = objstruct[allind].FLUX_OPT
inivar = objstruct[allind].IVAR_OPT
wave_all = objstruct[allind].WAVE_OPT
;; Mask junk wavelengths
bad_wave = WHERE(wave_all LE 3000.0 OR wave_all GT 1.2d4 $
                 OR finite(influx) EQ 0 OR finite(inivar) EQ 0, nbad)
IF nbad GT 0 THEN BEGIN
   inivar[bad_wave] = 0.0
   influx[bad_wave] = 0.0
ENDIF
inloglam =  alog10(wave_all)
;; data are extracted onto a grid so no shifting required
newloglam = inloglam[*, 0]

FOR kk = 0L, nimgs-1L DO BEGIN
   IF kk EQ 0 THEN splot, wave_all[*, kk], influx[*, kk] $
   ELSE soplot, wave_all[*, kk], influx[*, kk], col = (kk MOD 7)
   wait, 1
ENDFOR


IF KEYWORD_SET(LAM_MASK_RANGE) THEN BEGIN
   loglammask = newloglam
   masklam = lonarr(n_elements(newloglam))
   indx = WHERE(10.0d^loglammask GE LAM_MASK_RANGE[0] AND $
                10.0d^loglammask LE LAM_MASK_RANGE[1], nrange)
   IF nrange GT 0 THEN masklam[indx] = 1 $
   ELSE message, 'ERROR: No wavelengths in specified range'
ENDIF


long_combspec, influx, inivar, inloglam $
               , newloglam = newloglam, newflux = newflux $
               , newivar = newivar, newmask = newmask $
               , iref = iref, SIGREJ = SIGREJ, CHECK = CHK $
               , /NOSHIFT, /NOSHARP, NOREJ = NOREJ $
               , SN2 = SN2, YMULT = YMULT, DEBUG = DEBUG $
               , MEDSCALE = MEDSCALE $
               , LOGLAMMASK = LOGLAMMASK, MASKLAM = MASKLAM

;; Write combined spectrum out to a file
newlam = 10.0D^newloglam
IF keyword_set(OUTFIL) THEN BEGIN
   sxaddpar, head0, 'NEXP', nimgs
   sxaddpar, head0, 'EXPTIME_TOT', exptime
   sxaddpar, head0, 'BITPIX', -32
   sxaddpar, head0, 'NAXIS', 1
   sxaddpar, head0, 'NAXIS1', n_elements(newflux)
;    IF KEYWORD_SET(SKYFILE) AND nimgs GT 1 THEN $
;      sxaddpar, head0, 'FLX_SHFT_WAV', flx_shft_wav
    sxdelpar, head0, 'NAXIS2'
    sxdelpar, head0, 'BZERO'
    sxdelpar, head0, 'BSCALE'
    mwrfits, newflux, outfil, head0, /create
    giv = where(newivar GT 0., ngiv)
    sig = 0*newivar - 1.0D
    sig[giv] = 1./sqrt(newivar[giv])
    mwrfits, sig, outfil
    mwrfits, 10.0d^newloglam, outfil
;    gniv = where(newnivar GT 0., ngniv)
    mwrfits, sn2, outfil
    mwrfits, ymult, outfil
    print, 'long_coadd: Final file is ', outfil
ENDIF

RETURN
END







