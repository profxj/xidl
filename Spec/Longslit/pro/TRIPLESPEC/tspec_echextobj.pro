;+ 
; NAME:
; fire_echextobj   
;     Version 1.1
;
; PURPOSE:
;    Extract flux from 2D image to create ten 1D spectra (1 per order)
;    Output is written to the object structure (e.g. Extract/Obj_fire0024.fits)
;    The code only does boxcar extraction for now.
;
; CALLING SEQUENCE:
;   
;  fire_echextobj, fire, obj_id, [exp], /DEBUG, /CHK, /STD, APER=,
;  RADIUS=
;
; INPUTS:
;   fire  -  FIRE structure
;   indx  -  Indices of objects to process
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   /CHK    - Show intermidate diagnostics (e.g. profile)
;   /SPCCHK - Show extracted spectra
;   /STD    - Extraction should be set for a standard star
;   /DEBUG  - Stop within extraction routine to check stuff
;   APER=   - Set aperture by hand (e.g. [5., 7.] )
;   RADIUS= - Size of window for setting aperture size (default: 20L)
;   ORDRS=  - Orders to extract (default: [0L,9L])
;   /OPTIMAL -  Perform an optimal extraction
;   INPGAU  - Gaussian sigam to use for optimal [default: Code
;             calculates the best value]
;   ABSMAXGAU - maximum gaussian sigma allowed for any order 
;               (default: 7.0 pix or 2.4")
;   ABSMINGAU - maximum gaussian sigma allowed for any order 
;               (default: 1.0 or 0.35t")
;   LOWSNR   - Minimum S2N for performing profile fitting
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;  1)  The program begins extracting in order 9L (physical 6) and will
;  automatically calculate an aperture for that order.  If there is
;  insufficient flux in the following orders, it will adopt the value
;  from the next higher order. 
;
; EXAMPLES:
;   fire_echextobj, fire, 1L, [0L]
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   15-Aug-2002 Written by JXP
;   04-Feb-2003 Polished (JXP)
;   04-June-2008 Modificatios to improve extraction fits JFH
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO tspec_echextobj, AB, sky_resids, ivar_AB, scihdr, sky_model $
                     , piximg, waveimg, tset_slits, objstr $
                     , outfil = outfil, BOXCAR = BOXCAR $
                     , BOX_RAD = BOX_RAD1, STD = STD, CHK = CHK $
                     , AIRVAC_CORR = AIRVAC_CORR, HELIO_CORR = HELIO_CORR $
                     , SPCCHK = SPCCHK, DEBUG = DEBUG, VERBOSE = verbose

  plate_scale = 0.234D          ; TSPEC plate scale
  t1 = systime(1)
  spcchk = 1
  chk = 1
  ;; Added by JFH to prevent rejection from standard stars
  IF KEYWORD_SET(STD) THEN BEGIN
     REJ_CR = 0
     NOCRREJ = 0                ; We'll try leaving this on as an experiment
     BOXCAR = 1
     BOX_RAD1 = 14L
     REJSIG = 5.0
  ENDIF ELSE BEGIN
     REJ_CR = 1
     NOCRREJ = 0
  ENDELSE
  ;;  Optional Keywords
  ;; The gaussian S/N ratio is for whether we use bspline profile
  ;; or a gaussian with the right FWHM 
  IF NOT KEYWORD_SET(SN_GAUSS) THEN SN_GAUSS = 3.5d 
  IF NOT keyword_set(LOWSNR) then LOWSNR = 2.0
  ;; The LOWSNR determines if we fit a model from the other orders 
  ;; and also for using the brightest object on the slit 
  ;; to fit the fainter object 
  if n_elements(ORDRS1) EQ 0 THEN ordrs = lindgen(5) $
  ELSE ordrs = ordrs1
  if not keyword_set( CBIN ) then cbin = 1
  if not keyword_set( RBIN ) then rbin = 1
  if not keyword_set( LSLITE ) then lslite = round(3./cbin)
  if not keyword_set( RSLITE ) then rslite = round(3./cbin)
  IF NOT KEYWORD_SET(PROF_LIM) THEN PROF_LIM = 0.40 
  ;; fraction of good profile pixels required to not be masked
  IF NOT KEYWORD_SET(BOX_RAD1) THEN box_rad = 10L/cbin $
  ELSE BOX_RAD = BOX_RAD1
  SN_BOX_RAD=BOX_RAD/2.0 
  ;; SN_BOX_RAD used only for estimating SN ratio core only 1.5"
  alphabet = ['a', 'b', 'c', 'd', 'e', 'f']
  ;central wavelength of each order
  wave_ordr = [9350.0d, 10900.0d, 13070.0d, 16300.0d, 21700.0d]
  ;; make sure orders is sorted
  nord = n_elements(ordrs)
  outgau = fltarr(nord)
  icheck = WHERE(sort(ordrs) - lindgen(nord) NE 0, nbad)
  IF nbad GT 0 THEN message, 'The numbers in ordrs must be sorted'
  IF KEYWORD_SET(INPGAU) THEN BEGIN
     IF n_elements(inpgau) EQ nord THEN inpgau_vec = inpgau $
     ELSE IF n_elements(inpgau) EQ 1 THEN inpgau_vec = replicate(inpgau, nord) $
     ELSE message, 'invalid nubmer of elements in inpgau'
  ENDIF
  IF NOT KEYWORD_SET(ABSMAXGAU) THEN ABSMAXGAU = 7.0d
  IF NOT KEYWORD_SET(ABSMINGAU) THEN ABSMINGAU = 1.0d
  ;;  Find all relevant obj

  ;; This assumes FOWLER sampling. Update to be correct for TELLURICS??????
  dnsamp = double((strsplit(sxpar(scihdr, 'FSAMPLE'), /extract))[0])
  rdnoise = 14.0D/sqrt(dnsamp)

  dim = size(AB, /dimensions)
  nx = dim[0]
  ny = dim[1] 
  dim_tset=size(tset_slits[0].COEFF,/dim)
  norders=dim_tset[1]
  slit_edg = fltarr(ny, norders, 2)
  traceset2xy, tset_slits[0], rows, left_edge
  traceset2xy, tset_slits[1], rows, right_edge
  slit_edg[*, *, 0] = left_edge
  slit_edg[*, *, 1] = right_edge
  ;; ???? This is just a guess for now???????
  cdelt = 40.0d/299792.458d/alog(10.0d)
  npix = 3000 
  tot_wave = 10^(alog10(8000.d) + dindgen(8550)*cdelt)
  all_crval1 = dblarr(norders)
  all_crval1[0] = 8050.00
  all_crval1[1] = 9400.00
  all_crval1[2] = 11250.00
  all_crval1[3] = 14050.00
  all_crval1[4] = 18700.00

  for qq=0L,norders-1L do begin
     mn = min(abs(all_crval1[qq]-tot_wave),imn)
     all_crval1[qq] = alog10(tot_wave[imn])
  endfor
  all_coll = dblarr(norders,2)
  all_coll[0, *] = [8070., 10650.]
  all_coll[1, *] = [9400., 12400.]
  all_coll[2, *] = [11250., 14850.]
  all_coll[3, *] = [14050., 18550.]
  all_coll[4, *] = [18700., 25000.]

  max_hwidth = fltarr(norders)
  FOR kk = 0L, norders-1L DO $
    max_hwidth[kk] = max(slit_edg[*, kk, 1]-slit_edg[*, kk, 0])
  rnd_edg = round(slit_edg)

  ;; create some quantitites we will need
  ordermask = tspec_ordermask(tset_slits, order_vec = order_vec) 
  ;; this ordermask explicitly masks problematic regions for extraction
  ;; where orders fall off of CCD 

  ; Trimming by 3 pixels from each slit
  ximg = long_slits2x(tset_slits, edgmask = edgmask, nslit = nslit,TOL_EDG=3)
  ;;bsp = 0.4d                    ; Default 0.2, changes RAS for FIRE
  nccd = 1
  yarr = 0
  
  IF NOT KEYWORD_SET(BOXCAR) THEN BEGIN
     ;window, 0, xsize = 450, ysize = 370
     ;window, 1, xsize = 450, ysize = 370
     window, 0, xsize = 585, ysize = 481
     window, 1, xsize = 585, ysize = 481
     wset, 0
  ENDIF
  DEFGAU = (3.0d < ABSMAXGAU) > ABSMINGAU 
  IF KEYWORD_SET(OBJ_NM) THEN BEGIN
     nobj = n_elements(obj_nm)
     alpha_redux = obj_nm
  ENDIF ELSE BEGIN
     nobj = n_elements(objstr)/norders
     alpha_redux = alphabet[lindgen(nobj)]
  ENDELSE
  var = (ivar_AB GT 0.0)/(ivar_AB + (ivar_AB EQ 0.0))
  outmask = fltarr(nx, ny)
  modelivar = ivar_AB*(ordermask GT 0.0)
  objimage = fltarr(nx, ny)
  prof_img = fltarr(nx, ny)

  IF djs_median(waveimg[where(waveimg GT 0.0)]) LE 10.0 THEN BEGIN
     iwv=WHERE(waveimg GT 0.0)
     waveimg[iwv] = 10.0d^waveimg[iwv]
  ENDIF

  ;; AIR TO VAC CORRECTION
  if keyword_set(AIRVAC_CORR) then begin
     iwv=WHERE(waveimg GT 0)
     ;; NO AIR TO VAC CORRECTION NEEDED
     ;; WE ARE USING VACUUM OH WAVELENGTHS
     ;; AND FIRE"S GRATING IS IN VACUUM
     tmp_waveimg=waveimg[iwv]
     airtovac, tmp_waveimg
     waveimg[iwv]=tmp_waveimg
  endif
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;  HELIO correction
  if keyword_set(HELIO_CORR) then begin
  		;; Determine the Julian date
     jd = tspec_get_jd(scihdr, /REDUCED) ; Geocentric JD-2400000 (ie, 'reduced' Julian date)
     ra = sxpar(scihdr, "RA")
     dec = sxpar(scihdr, "DEC")
     x_radec, ra, dec, ra_deg, dec_deg
     jdhelio = helio_jd(jd, ra_deg, dec_deg) ;Returns heliocentric JD
     jdhelio += 2400000.d
    ; Heliocentric velocity correction in km/s

     ;!!!!!!!!!! This sign is definitely correct !!!!!!!!!
     helio_corr = (-1.0)*x_keckhelio(ra_deg, dec_deg, jd=jdhelio, obs='lco')
     gamma = sqrt( (1.d + helio_corr/299792.458d) / $
                   (1.d - helio_corr/299792.458d))
     waveimg = waveimg * gamma
     
     msg = 'Helio corr (NO vacuum) applied for '+sxpar(scihdr, 'OBJECT')+', helio = '+strtrim(helio_corr, 2)+' km/s'
     print, msg
     sxaddpar, scihdr, "HISTORY", msg
     sxaddpar, scihdr, "HELIO", strtrim(helio_corr, 2), 'km/s'

  endif


;; PREPROCESSING TO PRIORITIZE WHICH ORDERS TO EXTRACT FIRST
;; GO IN ORDER OF DECREASING SNR
;;
  print, systime()
  ;; Calculate S2N of each object for each order 
  ordr_s2n = fltarr(nord, nobj)
  FOR iobj = 0L, nobj-1L DO BEGIN
     FOR iorder = 0L, nord-1L DO BEGIN
        qq = ordrs[iorder]
        sci11 = where(objstr.obj_id EQ alpha_redux[iobj] $
                     AND objstr.order EQ qq)
        sci1=sci11[0]
        objcen   = objstr[sci1].trace[0:ny-1L]  
        ;; rename the objects 1 and 2 for each slit to simplify QA
        ;; use real-order number for QA. 
        objstr[sci1].OBJID = iobj + 1L
        objstr[sci1].SLITID = order_vec[qq]
        ;; extract asymmetric boxcars for S/N calculation
        left_edge   = (objcen - sn_box_rad/cbin) > (rnd_edg[*, qq, 0]+LSLITE)
        right_edge  = (objcen + sn_box_rad/cbin) < (rnd_edg[*, qq, 1]-RSLITE)
        ;; mask
        snmask = ordermask EQ order_vec[qq] AND ivar_AB GT 0.0 $
                 AND waveimg GT 0.0
        ;; boxcar extraction
        fx    = extract_asymbox2((AB-sky_resids)*snmask $
                                 , left_edge, right_edge)
        fvar  = extract_asymbox2(var*snmask, left_edge, right_edge)
        pixtot = extract_asymbox2(float(var*0 + 1), left_edge, right_edge)
        mask_box = extract_asymbox2(float(snmask EQ 0) $
                                    , left_edge, right_edge) NE pixtot
        fi = mask_box/(fvar + (fvar EQ 0))
        box_denom = extract_asymbox2((waveimg GT 0.0) $
                                     , left_edge, right_edge)
        wave  = extract_asymbox2(waveimg, left_edge, right_edge) $
                /(box_denom + (box_denom EQ 0))
        ;; median filter the resulting fx and fi to remove hot pixels
        flux_sm = djs_median(fx, width = 5, boundary = 'reflect')
        fluxivar_sm =  djs_median(fi, width = 5, boundary = 'reflect')
        fluxivar_sm = fluxivar_sm*(fi GT 0.0)
        indsp = WHERE(wave GT 9000.0 AND wave LT 2.46d4 $
                      AND NOT (wave GE 1.36d4 AND wave LE 1.496d4) $
                      AND NOT (wave GE 1.80d4 AND wave LE 1.995d4) $
                      AND finite(flux_sm) AND flux_sm LT 5.0d5 $
                      AND flux_sm GT -1000.0d $
                      AND fluxivar_sm GT 0.0, nsp)
        ;; spline fit the flux 
        b_answer = bspline_iterfit(wave[indsp], flux_sm[indsp] $
                                   , everyn = 1.5, yfit = spline_flux $
                                   , invvar = fluxivar_sm[indsp] $
                                   , upper = 5 $
                                   , lower = 5, /groupbadpix, maxrej = 1 $
                                   , outmask = bmask, /silent, /relative)
        b_answer = bspline_iterfit(wave[indsp], flux_sm[indsp] $
                                   , everyn = 1.5, yfit = spline_flux  $
                                   , invvar = fluxivar_sm[indsp]*bmask $
                                   , upper = 5 $
                                   , lower = 5, /groupbadpix, maxrej = 1 $
                                   , outmask = bmask2, /silent, /relative)
        sn2 = ((spline_flux*sqrt(fluxivar_sm[indsp] >  0)*bmask2) > 0)^2
        ind_nonzero = where(sn2 GT 0, nzero)
        IF nzero GT 0 THEN djs_iterstat, sn2[ind_nonzero] $
                                         , median = med_sn2 $
        ELSE med_sn2 = 0.0
        ordr_s2n[iorder, iobj] = med_sn2
        x_splot, wave, flux_sm, xtwo = wave[indsp], ytwo = flux_sm[indsp] $
                 , ythr = spline_flux, xthr = wave[indsp], /block
     ENDFOR
  ENDFOR
  ;; Compute the average S/N and find the brightest object
  s2n_bar = total(ordr_s2n, 1)
  srt_obj = reverse(sort(s2n_bar))
  ibright = srt_obj[0]
  ;; extract orders starting with highest S/N for brightest object
  srt_s2n = reverse(sort(ordr_s2n[*, ibright]))
  prof_fwhm = replicate(defgau, nord)
  fwhm_here = fltarr(nord)
  fwhm_was_fit = lonarr(nord)


  ;; Loop according to S2N
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;LOOP ON ORDERS
  splog, '*************************************************'
  splog, 'Extracting spectra in order of S/N:'
  print,  '   order#   S/N' 
  forprint, order_vec[ordrs[srt_s2n]], sqrt(ordr_s2n[srt_s2n, 0]) $
            , format = '(I6,F10.2)', textout = 1
  splog, '*************************************************'
  FOR jorder = 0L, nord-1L DO BEGIN
     ;; S2N
     iorder = srt_s2n[jorder]
     ;; Fun packing
     qq = ordrs[iorder]
     ;; Grab objects on this order
     ordind = where(objstr.order EQ qq)
     objstr[ordind].MASKWIDTH = max_hwidth[qq]/2.0d
     thismask = (ordermask EQ order_vec[qq]) AND waveimg GT 0.0
     ;; If S/N is lower than LOWSNR and at least two FWHM have 
     ;; been determined then fit FWHM vs order and use for profile
     other_orders = where(fwhm_here GT 0 AND fwhm_was_fit EQ 0, nother)
     other_fit    = where(fwhm_here GT 0 AND fwhm_was_fit EQ 1, nfit)
     
     ;; PRELIMINARIES FOR THIS ORDER: PLOT UP THE PROFILE, FIT PROFILE.
     FOR isrt = 0L, nobj-1L DO BEGIN
        iobj = srt_obj[isrt]
        IF ordr_s2n[iorder, iobj] LE LOWSNR AND nother GT 2 THEN BEGIN
           ;; If this is the brightest object, then we have to 
           ;; extrapolate the FWHM from a fit. 
           IF isrt EQ 0 THEN BEGIN
              splog, '*************************************************'
              splog, 'Doing linear fit for FWHM of object=' $
                     + alpha_redux[iobj] + ' in order# '$
                     + strcompress(string(order_vec[qq]), /rem)
              print,  '   order#   FWHM' 
              forprint, order_vec[ordrs[other_orders]] $
                        , fwhm_here[other_orders] $
                        , format = '(I6,F10.2)', textout = 1
              fwhm_coeffs = ladfit(order_vec[ordrs[other_orders]] $
                                   , fwhm_here[other_orders])
              def_fwhm = poly(order_vec[qq], fwhm_coeffs)
              def_fwhm = def_fwhm[0]
              fwhm_was_fit[iorder] = 1
              splog, 'FIT: ' + string(order_vec[qq], format = '(I4)') $
                     + string(def_fwhm, format = '(F10.2)')
              splog, '*************************************************'
              wset, 1
              plot, order_vec[ordrs[other_orders]], fwhm_here[other_orders] $
                    , psym = 4, symsize = 2, xr = [10.0, 33.0] $
                    , yrang = [0.7*min([fwhm_here, def_fwhm]) $
                               , 1.4*max([fwhm_here, def_fwhm])] $
                    , /xsty, /ysty
              clr = getcolor(/load)
              IF nfit GT 0 THEN oplot, order_vec[ordrs[other_fit]] $
                                       , fwhm_here[other_fit], psym = 6 $
                                       ,  symsize = 2 $
                                       , col = clr.cyan
              oplot, [order_vec[qq]], [def_fwhm], psym = 5, symsize = 2 $
                     , col = clr.green
              oplot, order_vec[ordrs], poly(order_vec[ordrs], fwhm_coeffs) $
                     , col = clr.red
              wset, 0
              sci1 = where(objstr.obj_id EQ alpha_redux[iobj] $
                           AND objstr.order EQ qq)
              objstr[sci1].FWHM = def_fwhm
           ENDIF ELSE BEGIN 
              sci1 = where(objstr.obj_id EQ alpha_redux[iobj] $
                           AND objstr.order EQ qq)
              sci1a = where(objstr.obj_id EQ alpha_redux[ibright] $
                            AND objstr.order EQ qq)
              objstr[sci1].FWHM =  objstr[sci1a].FWHM
           ENDELSE
        ENDIF
     ENDFOR ;;

     
     ;;;;;;;; OPTIMAL EXTRACTION CASE ;;;;;;;;;

     IF NOT KEYWORD_SET(BOXCAR) THEN BEGIN
        extract_struct = gnirs_extract(AB-sky_resids, ivar_AB, waveimg $
                                       , thismask, sky_model, objstr[ordind] $
                                       , plate_scale $
                                       , profile_struct = profile_struct)
        ;extract_struct = fire_localskysub(AB, ivar_AB, sky_model $
        ;                                  , piximg, waveimg, ximg, objstr $
        ;                                  , thismask $
        ;                                  , slit_edg[*, qq, 0] $
        ;                                  , slit_edg[*, qq, 1] $
        ;                                  , edgmask, bsp = bsp $
        ;                                  , objimage = objimage $
        ;                                  , modelivar = modelivar $
        ;                                  , outmask = outmask $
        ;                                  , indx = ordind, nccd = nccd $
        ;                                  , box_rad = box_rad $
        ;                                  , profile_struct =  $
        ;                                  profile_struct $
        ;                                  , scihdr = scihdr $
        ;                                  , SN_GAUSS = SN_GAUSS $
        ;                                  , CHK = LCHK, /SKYSAMPLE)
        ;; update the FWHM fitting vector for brightest object
        sci3 = where(objstr.obj_id EQ alpha_redux[ibright] $
                     AND objstr.order EQ qq)
        fwhm_here[iorder] = djs_median(objstr[sci3].FWHMFIT)
        ;; Did the FWHM get updated in the routine? 
        ;; If so, include this value in future fits. 
        IF abs(fwhm_here[iorder]-objstr[sci3].FWHM) GE 0.01 THEN $
           fwhm_was_fit[iorder] = 0
        ;; For cases where there are multiple objects on the slit 
        ;; loop over the objects that have low S/N ratio and 
        ;; if their fwhmfit is signifinicantly different from that 
        ;; measured from the brightest object, then rerun with this
        ;; new FWHM 
        IF nobj GT 1 THEN BEGIN
           RERUN_LOCAL = 0L
           FOR isrt = 1L, nobj-1L DO BEGIN
              iobj = srt_obj[isrt]
              sci4 = where(objstr.obj_id EQ alpha_redux[iobj] $
                           AND objstr.order EQ qq)
              fwhmfit = djs_median(objstr[sci4].FWHMFIT)
              fwhm_diff = abs((fwhmfit - fwhm_here[iorder]) $
                              /fwhm_here[iorder])
              ;; if brightest object has high S/N and fainter objects
              ;; have low S/N and there fwhm differs significantly 
              ;; from brightest then use brightest FWHM and rerun
              IF (profile_struct[ibright].MED_SN2 GT SN_GAUSS^2) AND $
                 (profile_struct[iobj].MED_SN2 LE LOWSNR^2) AND $
                 (fwhm_diff GT 0.2*fwhm_here[iorder]) THEN BEGIN
                 sci1 = where(objstr.obj_id EQ alpha_redux[iobj] $
                              AND objstr.order EQ qq)
                 objstr[sci1].FWHM = fwhm_here[iorder]
                 splog, '----------------------------------------'
                 splog, 'RE-RUNNING local sky subtraction because' 
                 splog, 'object=' + alpha_redux[iobj] + ' has'
                 splog, 'FWHM = ' + $
                        strcompress(string(fwhmfit, format = '(F4.2)') $
                                    , /rem) + $
                        '  which differs signifncantly from brightest' $
                        + ' object which has FWHM = ' + $
                        strcompress(string(fwhm_here[iorder] $
                                           , format = '(F4.2)'), /rem) 
                 splog, '----------------------------------------'
                 RERUN_LOCAL = 1L
              ENDIF
           ENDFOR
           IF RERUN_LOCAL EQ 1 THEN  $
              extract_struct = gnirs_extract(AB-sky_resids, ivar_AB, waveimg $
                                             , thismask, sky_model $
                                             , objstr[ordind] $
                                             , plate_scale $
                                             , profile_struct = profile_struct)

          ; extract_struct = fire_localskysub(AB, ivar_AB, sky_model $
          ;                                      , piximg, waveimg, objstr $
          ;                                      , thismask $
          ;                                      , slit_edg[*, qq, 0] $
          ;                                      , slit_edg[*, qq, 1] $
          ;                                      , edgmask, bsp = bsp $
          ;                                      , objimage = objimage $
          ;                                      , modelivar = modelivar $
          ;                                      , outmask = outmask $
          ;                                      , indx = ordind $
          ;                                      , nccd = nccd $
          ;                                      , box_rad = box_rad $
          ;                                      , profile_struct = $
          ;                                      profile_struct $
          ;                                      , scihdr = scihdr $
          ;                                      , SN_GAUSS = SN_GAUSS $
          ;                                      , CHK = LCHK)
        ENDIF
     ENDIF ;; END OF OPTIMAL EXTRACTION PROFILE FITTING

     ;; Now loop over the objects and extract each using profile
     FOR isrt = 0L, nobj-1L DO BEGIN
        iobj = srt_obj[isrt]
        print, '-----------------------------------------------'
        print, 'fire_echextobj: Extracting order ' $
               , string(order_vec[qq], FORMAT = '(i3)') +  $
               ' with  wave_cen = ' +  $
               string(wave_ordr[qq], FORMAT = '(I5)')
        print, '               Object = '  + alpha_redux[iobj]
        IF NOT KEYWORD_SET(BOXCAR) THEN BEGIN
           
           mincol = profile_struct[iobj].MINCOL
           maxcol = profile_struct[iobj].MAXCOL
        ENDIF ELSE BEGIN
           lhs = (rnd_edg[*, qq, 0]+LSLITE) > 0L 
           rhs = (rnd_edg[*, qq, 1]-RSLITE) < (nx-1L)
           mincol = min(round(lhs))
           maxcol = max(round(rhs))
        ENDELSE
        sci2 = where(objstr.obj_id EQ alpha_redux[iobj] $
                     AND objstr.order EQ qq)
        ;; Transpose
        timg = transpose((AB[mincol:maxcol, *]-sky_resids[mincol:maxcol, *])* $
                         thismask[mincol:maxcol, *])
        tarc = transpose(waveimg[mincol:maxcol, *]* $
                         thismask[mincol:maxcol, *])
        modelvar = (modelivar GT 0.0)/(modelivar + (modelivar LE 0.0))
        ;; explicit masking
        tsky = transpose(sky_model[mincol:maxcol, *]* $
                         thismask[mincol:maxcol, *])

        ; The last order has a negative y
        ; value for the trace where it falls
        ; off the chip

        IF NOT KEYWORD_SET(BOXCAR) THEN BEGIN ;; OPTIMAL EXTRACTION
           modelvar = (modelivar GT 0.0)/(modelivar + (modelivar LE 0.0))
           ;; explicit masking
           tvar = transpose(modelvar[mincol:maxcol, *])
           ;tvar = transpose(modelvar[mincol:maxcol, *]* $
           ;                 outmask[mincol:maxcol, *])
           spots = where(tarc EQ 0.0, nspots)
           if nspots NE 0 then tvar[ spots ] = -1.0

           ;;if (31L-qq EQ 31) then begin
           ;;   ytrace = (objstr[sci2].trace[0:ny-1L]-mincol > 0)
           ;;   wv_trace = tarc[indgen(2048),ytrace]
           ;;   min_wv = min(wv_trace[where(ytrace GT 0)])
           ;;   max_wv = max(wv_trace[where(ytrace GT 0)])
           ;;   trc_mask = where(tarc LT min_wv OR tarc GT max_wv,ntrcmsk)
           ;;   if (ntrcmsk GT 0) then begin
           ;;      tvar[trc_mask] = -1
           ;;   endif
           ;;endif

           tprof = $
              transpose(profile_struct[iobj].PROFILE[mincol:maxcol, *]* $
                        thismask[mincol:maxcol, *]) 
           prof_img = prof_img + profile_struct[iobj].PROFILE
           ;; Extract
           x_extobjopt, timg, tarc, tvar, tsky, PROFILE = tprof $
                        , [objstr[sci2].xcen, objstr[sci2].ycen-mincol] $
                        , fin_spec, COLLMNX = all_coll[qq, *] $
                        , CRVAL1 = all_crval1[qq], CDELT = cdelt $
                        , NPIX = npix $
                        , TOT_TRC = (objstr[sci2].trace[0:ny-1L]-mincol > 0) $
                        , /REBINC, READNO = rdnoise $
                        , REDBLUE = 0L, PROF_LIM = PROF_LIM $
                        , CHK = chk, DEBUG = debug
           ;; PROF_LIM added by JFH 06/08. This explicitly masks
           ;; any region for which the fraction of good pixels 
           ;; on the 2-d object profile frac_mask < PROF_LIM =
           ;; 0.40
        ENDIF ELSE BEGIN ;; BOXCAR EXTRACTION
           tvar = transpose(var[mincol:maxcol, *])
           spots = where(tarc EQ 0.0, nspots)
           if nspots NE 0 then tvar[ spots ] = -1.0

           ;;if (31L-qq EQ 31) then begin
           ;;   ytrace = (objstr[sci2].trace[0:ny-1L]-mincol > 0)
           ;;   wv_trace = tarc[indgen(2048),ytrace]
           ;;   min_wv = min(wv_trace[where(ytrace GT 0 AND tarc NE 0.0)])
           ;;   max_wv = max(wv_trace[where(ytrace GT 0 AND tarc NE 0.0)])
           ;;   trc_mask = where(tarc LT min_wv OR tarc GT max_wv,ntrcmsk)
           ;;   if (ntrcmsk GT 0) then begin
           ;;      tvar[trc_mask] = -1
           ;;   endif
           ;;endif
           f_extobjbox, timg-tsky, tarc, VAR = tvar $
                        , [objstr[sci2].xcen, objstr[sci2].ycen-mincol] $
                        , fin_spec, msk = msk $
                        , DEBUG = debug, wvmnx = [8000., 26000.] $
                        , COLLMNX = all_coll[qq, *] $
                        , CRVAL1 = all_crval1[qq], CDELT = cdelt $
                        , NPIX = npix, REJ_CR = REJ_CR $
                        ;; The lowest orders often hit the bottom of
                        ;; the image, so we need to make sure that
                        ;; tot_trc doesn't dip below 0.0 outside the
                        ;; range of the order
                        , TOT_TRC = ((objstr[sci2].trace[0:ny-1L]-mincol) $
                                     > 0.0 ) $
                        , REJSIG = rejsig, /REBINC, RADIUS = box_rad $
                        , NOCRREJ = NOCRREJ $
                        , BKAPER = svaper, CHK = chk 
           ;; Save the aperture
           svaper=fin_spec.aper
        ENDELSE
        if keyword_set(SPCCHK) AND fin_spec.npix NE 0 then begin
           x_splot, fin_spec.wv, fin_spec.fx $
                    , YTWO = sqrt(fin_spec.var > 0.), ymnx=[-1000,10000], /block
        endif
        ;; Write to structure
        if fin_spec.npix NE 0 then begin
           indf=lindgen(fin_spec.npix)
           objstr[sci2].npix = fin_spec.npix
           objstr[sci2].wave[indf] = fin_spec.wv
           objstr[sci2].fx[indf] = fin_spec.fx
           objstr[sci2].var[indf] = fin_spec.var
           igood=where(fin_spec.var GT 0.0,ngood,comp=ibad,ncomp=nbad)
           IF ngood GT 0 THEN $
              objstr[sci2].sig[indf[igood]]= sqrt(fin_spec.var[igood])
           IF nbad GT 0 THEN objstr[sci2].sig[indf[ibad]]=0.0
           objstr[sci2].novar[indf] = fin_spec.novar
           objstr[sci2].sky[indf] = fin_spec.sky
           objstr[sci2].sky_cts[indf] = fin_spec.sky ; This to hold sky in ADU/gain
           ;; Box
           if KEYWORD_SET(BOXCAR) THEN BEGIN
              objstr[sci2].box_wv[0:fin_spec.npix-1] = fin_spec.wv
              objstr[sci2].box_fx[0:fin_spec.npix-1] = fin_spec.fx
              objstr[sci2].box_var[0:fin_spec.npix-1] = fin_spec.var
           endif
           ;; Aper
           objstr[sci2].aper = fin_spec.aper
           objstr[sci2].gauss_sig = fin_spec.gauss_sig
           ;; Flag
           objstr[sci2].flg_anly = 1
	   ;; Convert '' to ' ' in order to avoid an error in the saved fits file (mrwfits should do this, but fails)
	   if objstr[sci2].img_fil EQ '' then begin 
		objstr[0].img_fil = ' '
	   endif
	   if objstr[sci2].arc_fil EQ '' then begin
		objstr[0].arc_fil = ' '
	   endif
	   if objstr[sci2].ut EQ '' then begin
		objstr[0].ut = ' '
	   endif
        endif else objstr[sci2].flg_anly = 0
     endfor
  ENDFOR ;; END ORDER LOOP
  print, systime()
  stop
  ;; Write reduction outputs
  if NOT keyword_set( objfile ) then begin
     objstrfil = strtrim(reduxdir,2)+"/Object/Obj_"+(strsplit(outfil, "_", /extract, count=count))[count-1]
     ;; Determine the directory name, and make the directory
     objstrdir = strtrim(reduxdir,2)+"/Object/"
     FILE_MKDIR, objstrdir
  endif else begin
     objstrfil = objfile
  endelse
  print, 'fire_echextobj: Writing object structure to: ', objstrfil
  
  objhdr = scihdr
  sxdelpar, objhdr, "BITPIX"
  sxdelpar, objhdr, "BITPIX"
  sxdelpar, objhdr, "NAXIS"
  sxdelpar, objhdr, "NAXIS1"
  sxdelpar, objhdr, "NAXIS2"
  mwrfits, objstr, objstrfil, objhdr, /create

	if NOT keyword_set( ECHEXTFILE ) AND keyword_set( OUTFIL ) then begin
	  print, 'fire_echextobj: Writing output to: ', outfil
	  ;; Determine the directory name, and make the directory
	  split = strsplit(outfil, '/', /EXTRACT, count=count)
	  if count NE 1 then begin
			outdir = strjoin(split[0:count-2], '/') + "/"
	  endif
	  FILE_MKDIR, outdir
	endif else if keyword_set( ECHEXTFILE ) then begin
		outfil = echextfile
	endif

  IF keyword_set(CHK) THEN BEGIN
     print, 'fire_echextobj:  Displaying sky subtracted image'
     xatv, (AB - sky_model)*float(ordermask GT 0.0) $
           , WVIMG = waveimg, /block, min = -50, max = 200
     print, 'fire_echextobj:  Displaying residual image'
     xatv, (AB - sky_model - objimage)*sqrt(modelivar)* $
           float(ordermask GT 0.0), WVIMG = waveimg, /block $
           , min = -6.0, max = 6.0
  ENDIF
  IF KEYWORD_SET(OUTFIL) THEN BEGIN
     mwrfits, float(AB), outfil, scihdr, /create, /silent
     mwrfits, float(ivar_AB), outfil, /silent
     mwrfits, float(sky_model), outfil, /silent
     mwrfits, float(piximg), outfil, /silent
     mwrfits, float(modelivar), outfil, /silent
     mwrfits, float(objimage), outfil, /silent
     mwrfits, float(outmask), outfil, /silent
     mwrfits, float(prof_img), outfil, /silent
     ;; Compress
     print, 'fire_echextobj: Compressing...'
     spawn, 'gzip -f '+ outfil
  ENDIF
  ;  DONE

  print, 'fire_echextobj: Elapsed time: '+strtrim((systime(1)-t1)/60.0,2)+' minutes'

  if keyword_set(DEBUG ) then begin
  	print, 'fire_echextobj: Placing final stop (just in case the file did not save...)'
  	stop
  endif

  print, 'fire_echextobj: All done! '

;  return,objstr
  RETURN
END
