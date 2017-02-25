;+ 
; NAME:
; esi_echextobj   
;     Version 1.1
;
; PURPOSE:
;    Extract flux from 2D image to create ten 1D spectra (1 per order)
;    Output is written to the object structure (e.g. Extract/Obj_esi0024.fits)
;    The code only does boxcar extraction for now.
;
; CALLING SEQUENCE:
;   
;  mage_echextobj, esi, obj_id, [exp], /DEBUG, /CHK, /STD, APER=,
;  RADIUS=
;
; INPUTS:
;   esi   -  ESI structure
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
;   mage_echextobj, esi, 1L, [0L]
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

PRO mage_echextobj,sciimg,sciivar,scihdr,skyimage $
                   , piximg,waveimg,tset_slits,objstr $
                   , outfil=outfil,BOXCAR=BOXCAR $
                   , BOX_RAD=BOX_RAD1,STD=STD,CHK=CHK $
                   , NOVACHELIO=NOVACHELIO

  ;; Added by JFH to prevent rejection from standard stars
  IF KEYWORD_SET(STD) THEN BEGIN
      REJ_CR = 0
      NOCRREJ = 1
      BOXCAR = 1
      BOX_RAD1 = 14L
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
  if n_elements(ORDRS1) EQ 0 THEN ordrs = lindgen(15) $
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
  wave_ordr = [3150.0,3200.0,3450.0,3650.0,3850.0,4150.0,4400.0,4750.0,5150 $
               ,5600,6150,6850.0,7700.0,8800.0,9750.0]
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
  rdnoise = float(sxpar(scihdr,'ENOISE'))

  dim = size(sciimg, /dimensions)
  nx = dim[0]
  ny = dim[1] 
  dim_tset=size(tset_slits[0].COEFF,/dim)
  norders=dim_tset[1]
  slit_edg = fltarr(ny, 15, 2)
  traceset2xy, tset_slits[0], rows, left_edge
  traceset2xy, tset_slits[1], rows, right_edge
  slit_edg[*, *, 0] = left_edge
  slit_edg[*, *, 1] = right_edge
;;  Wavelength solution (put in a file)
  cdelt = 22.0d/299792.458d/alog(10.0d)
  npix = 2500 
;;  cdelt = 0.00001447624d
  tot_wave = 10^(alog10(3000.d) + dindgen(17720L)*cdelt)
  all_crval1 = dblarr(15)
  all_crval1[0] = 3000.0
  all_crval1[1] = 3025.0
  all_crval1[2] = 3140.0
  all_crval1[3] = 3330.0
  all_crval1[4] = 3540.0
  all_crval1[5] = 3780.0
  all_crval1[6] = 4050.0
  all_crval1[7] = 4365.0
  all_crval1[8] = 4730.0
  all_crval1[9] = 5160.0
  all_crval1[10] = 5680.0
  all_crval1[11] = 6310.0
  all_crval1[12] = 7105.0
  all_crval1[13] = 8125.0
  all_crval1[14] = 9490.0
  for qq=0L,norders-1L do begin
     mn = min(abs(all_crval1[qq]-tot_wave),imn)
     all_crval1[qq] = alog10(tot_wave[imn])
  endfor
  all_coll = dblarr(15,2)
  all_coll[0,*] = [3100., 3200.]
  all_coll[1,*] = [3200., 3400.]
  all_coll[2,*] = [3300., 3575.]
  all_coll[3,*] = [3500., 3800.]
  all_coll[4,*] = [3700., 4050.]
  all_coll[5,*] = [3950., 4350.]
  all_coll[6,*] = [4150., 4650.]
  all_coll[7,*] = [4500., 5050.]
  all_coll[8,*] = [4800., 5500.]
  all_coll[9,*] = [5250., 6000.]
  all_coll[10,*] = [5750., 6600.]
  all_coll[11,*] = [6350., 7350.]
  all_coll[12,*] = [7150., 8250.]
  all_coll[13,*] = [8150., 9450.]
  all_coll[14,*] = [9500., 9980.]


  max_hwidth = fltarr(norders)
  FOR kk = 0L, norders-1L DO $
    max_hwidth[kk] = max(slit_edg[*, kk, 1]-slit_edg[*, kk, 0])
  rnd_edg = round(slit_edg)

  ;; create some quantitites we will need
  ordermask = mage_ordermask(tset_slits) 
  ;; this ordermas explicitly masks problematic regions for extraction
  ;; where orders fall off of CCD or edge of order is too red (order 6)
  ximg = long_slits2x(tset_slits, edgmask = edgmask, nslit = nslit,TOL_EDG=0)
  bsp = 0.6d
  nccd = 1
  yarr = 0
  
  IF NOT KEYWORD_SET(BOXCAR) THEN BEGIN
     window, 0
     window, 1
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
  var = (sciivar GT 0.0)/(sciivar + (sciivar EQ 0.0))
  outmask = fltarr(nx, ny)
  modelivar = sciivar*(ordermask GT 0.0)
  objimage = fltarr(nx, ny)
  prof_img = fltarr(nx, ny)

  IF djs_median(waveimg[where(waveimg GT 0.0)]) LE 10.0 THEN BEGIN
     iwv=WHERE(waveimg GT 0.0)
     waveimg[iwv] = 10.0d^waveimg[iwv]
  ENDIF
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;  HELIO correction
  IF NOT keyword_set(NOVACHELIO) then begin
     iwv=WHERE(waveimg GT 0)
     tmp_waveimg=waveimg[iwv]
     airtovac, tmp_waveimg
     waveimg[iwv]=tmp_waveimg
     date_day = float(strsplit(sxpar(scihdr,"UT-DATE"),"-", /extract))
     date_time = float(strsplit(sxpar(scihdr,"UT-TIME"),":", /extract))
     date = [date_day[0],date_day[1],date_day[2],date_time[0],date_time[1],date_time[2]]
     ra_deg = sxpar(scihdr, "RA-D")
     dec_deg = sxpar(scihdr, "DEC-D")
     juldate, date, jd                     ; Geocentric JD-2400000
     jdhelio = helio_jd(jd, ra_deg, dec_deg) ;Returns heliocentric JD
     jdhelio += 2400000.d

    ; Heliocentric velocity correction in km/s
     helio_corr = (-1.0)*x_keckhelio(ra_deg, dec_deg, jd=jdhelio, obs='lco')
     gamma = sqrt( (1.d + helio_corr/299792.458d) / $
                   (1.d - helio_corr/299792.458d))
     waveimg = waveimg * gamma

     msg = 'Vac/helio applied for '+sxpar(scihdr,"FILENAME")+', helio = '+strtrim(helio_corr,2)+' km/s'
     print, msg
     sxaddpar, scihdr, "COMMENT", msg
  endif

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
        objstr[sci1].SLITID = 20L-qq
        ;; extract asymmetric boxcars for S/N calculation
        left_edge   = (objcen - sn_box_rad/cbin) > (rnd_edg[*, qq, 0]+LSLITE)
        right_edge  = (objcen + sn_box_rad/cbin) < (rnd_edg[*, qq, 1]-RSLITE)
        ;; mask
        snmask = ordermask EQ (20-qq) AND sciivar GT 0.0 AND waveimg GT 0.0
        ;; boxcar extraction
        fx    = extract_asymbox2((sciimg-skyimage)*snmask $
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
        indsp = WHERE(wave GT 3200.0 AND wave LT 1.0d4 $
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
  forprint, 20L-ordrs[srt_s2n], sqrt(ordr_s2n[srt_s2n, 0]) $
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
     thismask = (ordermask EQ (20-qq)) AND waveimg GT 0.0
     ;; If S/N is lower than LOWSNR and at least two FWHM have 
     ;; been determined then fit FWHM vs order and use for profile
     other_orders = where(fwhm_here GT 0 AND fwhm_was_fit EQ 0, nother)
     other_fit    = where(fwhm_here GT 0 AND fwhm_was_fit EQ 1, nfit)
     FOR isrt = 0L, nobj-1L DO BEGIN
        iobj = srt_obj[isrt]
        IF ordr_s2n[iorder, iobj] LE LOWSNR AND nother GT 2 THEN BEGIN
           ;; If this is the brightest object, then we have to 
           ;; extrapolate the FWHM from a fit. 
           IF isrt EQ 0 THEN BEGIN
              splog, '*************************************************'
              splog, 'Doing linear fit for FWHM of object=' $
                     + alpha_redux[iobj] + ' in order# '$
                     + strcompress(string(20L-qq), /rem)
              print,  '   order#   FWHM' 
              forprint, 20L-ordrs[other_orders] $
                        , fwhm_here[other_orders] $
                        , format = '(I6,F10.2)', textout = 1
              fwhm_coeffs = ladfit(20-ordrs[other_orders] $
                                   , fwhm_here[other_orders])
              def_fwhm = poly(20L-qq, fwhm_coeffs)
              def_fwhm = def_fwhm[0]
              fwhm_was_fit[iorder] = 1
              splog, 'FIT: ' + string(20L-qq, format = '(I4)') $
                     + string(def_fwhm, format = '(F10.2)')
              splog, '*************************************************'
              wset, 1
              plot, 20L-ordrs[other_orders], fwhm_here[other_orders] $
                    , psym = 4, symsize = 2, xr = [5.0, 21.0] $
                    , yrang = [0.7*min([fwhm_here, def_fwhm]) $
                               , 1.4*max([fwhm_here, def_fwhm])] $
                    , /xsty, /ysty
              clr = getcolor(/load)
              IF nfit GT 0 THEN oplot, 20L-ordrs[other_fit] $
                                       , fwhm_here[other_fit], psym = 6, symsize = 2 $
                                       , col = clr.cyan
              oplot, [20L-qq], [def_fwhm], psym = 5, symsize = 2 $
                     , col = clr.green
              oplot, 20L-ordrs, poly(20L-ordrs, fwhm_coeffs) $
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
     ENDFOR
     IF NOT KEYWORD_SET(BOXCAR) THEN BEGIN
        extract_struct = long_localskysub(sciimg, sciivar, skyimage $
                                          , piximg, waveimg, ximg, objstr $
                                          , thismask $
                                          , slit_edg[*, qq, 0] $
                                          , slit_edg[*, qq, 1] $
                                          , edgmask, bsp = bsp $
                                          , objimage = objimage $
                                          , modelivar = modelivar $
                                          , outmask = outmask $
                                          , indx = ordind, nccd = nccd $
                                          , box_rad = box_rad $
                                          , profile_struct =  $
                                          profile_struct $
                                          , scihdr = scihdr $
                                          , SN_GAUSS = SN_GAUSS $
                                          , CHK = LCHK)
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
              extract_struct = long_localskysub(sciimg, sciivar, skyimage $
                                                , piximg, waveimg, objstr $
                                                , thismask $
                                                , slit_edg[*, qq, 0] $
                                                , slit_edg[*, qq, 1] $
                                                , edgmask, bsp = bsp $
                                                , objimage = objimage $
                                                , modelivar = modelivar $
                                                , outmask = outmask $
                                                , indx = ordind $
                                                , nccd = nccd $
                                                , box_rad = box_rad $
                                                , profile_struct = $
                                                profile_struct $
                                                , scihdr = scihdr $
                                                , SN_GAUSS = SN_GAUSS $
                                                , CHK = LCHK)
        ENDIF
     ENDIF
     ;; Now loop over the objects and extract each using profile
     FOR isrt = 0L, nobj-1L DO BEGIN
        iobj = srt_obj[isrt]
        print, '-----------------------------------------------'
        print, 'mage_echextobj: Extracting order ' $
               , string(20L-qq, FORMAT = '(i3)') +  $
               ' with  wave_cen = ' +  $
               string(wave_ordr[qq], FORMAT = '(I4)')
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
        timg = transpose(sciimg[mincol:maxcol, *]* $
                         thismask[mincol:maxcol, *])
        tarc = transpose(waveimg[mincol:maxcol, *]* $
                         thismask[mincol:maxcol, *])
        modelvar = (modelivar GT 0.0)/(modelivar + (modelivar LE 0.0))
        ;; explicit masking
        tvar = transpose(modelvar[mincol:maxcol, *]* $
                         outmask[mincol:maxcol, *])
        tsky = transpose(skyimage[mincol:maxcol, *]* $
                         thismask[mincol:maxcol, *])
        IF NOT KEYWORD_SET(BOXCAR) THEN BEGIN ;; OPTIMAL EXTRACTION
           modelvar = (modelivar GT 0.0)/(modelivar + (modelivar LE 0.0))
           ;; explicit masking
           tvar = transpose(modelvar[mincol:maxcol, *]* $
                            outmask[mincol:maxcol, *])
           tprof = $
              transpose(profile_struct[iobj].PROFILE[mincol:maxcol, *]* $
                        thismask[mincol:maxcol, *]) 
           prof_img = prof_img + profile_struct[iobj].PROFILE
           ;; Extract
           x_extobjopt, timg-tsky, tarc, tvar, tsky, PROFILE = tprof $
                        , [objstr[sci2].xcen, objstr[sci2].ycen-mincol] $
                        , fin_spec, COLLMNX = all_coll[qq, *] $
                        , CRVAL1 = all_crval1[qq], CDELT = cdelt $
                        , NPIX = npix $
                        , TOT_TRC = objstr[sci2].trace[0:ny-1L]-mincol $
                        , /REBINC, READNO = rdnoise $
                        , REDBLUE = 0L, PROF_LIM = PROF_LIM $
                        , CHK = chk, DEBUG = debug
           ;; PROF_LIM added by JFH 06/08. This explicitly masks
           ;; any region for which the fraction of good pixels 
           ;; on the 2-d object profile frac_mask < PROF_LIM =
           ;; 0.40
        ENDIF ELSE BEGIN ;; BOXCAR EXTRACTION
           tvar = transpose(var[mincol:maxcol, *])
           x_extobjbox, timg-tsky, tarc, VAR = tvar $
                        , [objstr[sci2].xcen, objstr[sci2].ycen-mincol] $
                        , fin_spec $
                        , APER = aper, DEBUG = debug $
                        , COLLMNX = all_coll[qq, *] $
                        , CRVAL1 = all_crval1[qq], CDELT = cdelt $
                        , NPIX = npix, REJ_CR = REJ_CR $
                        , TOT_TRC = objstr[sci2].trace[0:ny-1L]-mincol $
                        , REJSIG = rejsig, /REBINC, RADIUS = box_rad $
                        , NOCRREJ = NOCRREJ $
                        , BKAPER = svaper, CHK = chk 
           ;; Save the aperture
           svaper=fin_spec.aper
        ENDELSE
        if keyword_set(SPCCHK) AND fin_spec.npix NE 0 then begin
           x_splot, fin_spec.wv, fin_spec.fx $
                    , YTWO = sqrt(fin_spec.var > 0.), /block
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
        endif else objstr[sci2].flg_anly = 0
     endfor
  endfor
  print, systime()
  ;; Write reduction outputs
  print, 'mage_echskysub: Writing output to: ', outfil
  IF keyword_set(CHK) THEN BEGIN
     print, 'mage_echextobj:  Displaying sky subtracted image'
     xatv, (sciimg - skyimage)*float(ordermask GT 0.0) $
           , WVIMG = waveimg, /block, min = -50, max = 200
     print, 'mage_echextobj:  Displaying residual image'
     xatv, (sciimg - skyimage - objimage)*sqrt(modelivar)* $
           float(ordermask GT 0.0), WVIMG = waveimg, /block $
           , min = -6.0, max = 6.0
  ENDIF
  IF KEYWORD_SET(OUTFIL) THEN BEGIN
     mwrfits, float(sciimg), outfil, scihdr, /create, /silent
     mwrfits, float(sciivar), outfil, /silent
     mwrfits, float(skyimage), outfil, /silent
     mwrfits, float(piximg), outfil, /silent
     mwrfits, float(modelivar), outfil, /silent
     mwrfits, float(objimage), outfil, /silent
     mwrfits, float(outmask), outfil, /silent
     mwrfits, float(prof_img), outfil, /silent
     ;; Compress
     print, 'mage_echskysub: Compressing...'
     spawn, 'gzip -f '+outfil
  ENDIF
  ;  DONE
  print, 'mage_echextobj: All done! '
;  return,objstr
  RETURN
END
