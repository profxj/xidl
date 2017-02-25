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
;  esi_echextobj, esi, obj_id, [exp], /DEBUG, /CHK, /STD, APER=,
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
;   /OPTIMAL -  Perform an optimal extraction -- Deprecated
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
;   esi_echextobj, esi, 1L, [0L]
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   15-Aug-2002 Written by JXP
;   04-Feb-2003 Polished (JXP)
;   04-June-2008 Modificatios to improve extraction fits JFH
;   19-Jan-2015 added ximg to RERUN_LOCAL long_localskysub call, KLC
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echextobj, esi, obj_id, exp, DEBUG = debug, CHK = chk $
                   , STD = std, APER = aper, BOX_RAD = BOX_RAD $
                   , BOXCAR = boxcar $
                   , ORDRS = ordrs1, AIMG=aimg $
                   , SEDG_FIL = sedg_fil, OPTIMAL = optimal, INPGAU = inpgau $
                   , FIT_APER = FIT_APER $
                   , SPCCHK = spcchk, ABSMAXGAU = ABSMAXGAU $
                   , ABSMINGAU = ABSMINGAU $
                   , _EXTRA = extra $
                   , CBIN = cbin, RBIN = rbin, NITER = niter $
                   , LSLITE = LSLITE, RSLITE = RSLITE $
                   , OBJ_NM = OBJ_NM, OUTGAU = OUTGAU, SN_GAUSS = SN_GAUSS $
                   , DEFRINGE = DEFRINGE, FRINGEMASK = FRINGEMASK $
                   , NOFRINGEMASK = NOFRINGEMASK

  print, systime()

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echextobj, esi, obj_id, [exspr], RADIUS=, APER=, /DEBUG, /CHK'
      print, '          /STD, ORDRS=, /SPCCHK, /OPTIMAL [v1.1]'
      return
  endif 

  ;; size of mask for defringed images in  units of FWHM
  IF NOT KEYWORD_SET(FRINGEMASK) THEN FRINGEMASK = 2.0
  ;; Added by JFH to prevent rejection from standard stars
  IF KEYWORD_SET(STD) THEN BEGIN
      REJ_CR = 0
      NOCRREJ = 1
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
  if n_elements(ORDRS1) EQ 0 THEN ordrs = lindgen(10) $
  ELSE ordrs = ordrs1
  if not keyword_set( CBIN ) then cbin = 1
  if not keyword_set( RBIN ) then rbin = 1
  if not keyword_set( LSLITE ) then lslite = round(17./cbin)
  if not keyword_set( RSLITE ) then rslite = round(9./cbin)
  IF NOT KEYWORD_SET(PROF_LIM) THEN PROF_LIM = 0.40 
  IF NOT KEYWORd_SET(BOX_RAD) THEN box_rad = 10L/cbin
  ;; fraction of good profile pixels required to not be masked JFH 05/08
  alphabet = ['a', 'b', 'c', 'd', 'e', 'f']
  plate_scale = reverse([0.168, 0.163, 0.158, 0.153, 0.149, 0.144, 0.137 $
                         , 0.134, 0.127, 0.120])*cbin
  wave_ordr   = reverse([9800.0, 8750.0, 7800.0, 6850.0, 6200.0, 5650.0 $
                         , 5200.0, 4800.0, 4500.0, 4200.0])
  ;; Maximum FWHM = 1.5" seeing
  IF NOT KEYWORD_SET(MAXFWHM) THEN MAXFWHM = 1.5D/plate_scale
  ;; Minimum FWHM is 1 pixel
  IF NOT KEYWORD_SET(MINFWHM) THEN MINFWHM = replicate(2.0D, 10L)
    
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
  if not keyword_set(STD) then begin
      indx = where(esi.flg_anly NE 0 AND esi.mode EQ 2 AND $
                   esi.rbin EQ rbin AND esi.cbin EQ cbin AND $
                   esi.obj_id EQ obj_id AND strtrim(esi.type,2) EQ 'OBJ', nindx)
      if nindx EQ 0 then begin
          print, 'esi_echextobj: No images to find obj for!', obj_id
          return
      endif
  endif else begin  ; STD star
      indx = obj_id[0]
      nindx = 1L
      BOXCAR = 1L
      box_rad = 40L
  endelse

;  Exposures
  if n_elements(exp) EQ 0 then exp = lindgen(nindx)

;  Wavelength solution (put in a file)
  cdelt = 0.00001447624d
  tot_wave = 10^(alog10(3900.d) + dindgen(33000L)*cdelt)
  all_crval1 = dblarr(10)
  all_crval1[0] = 3900.
  all_crval1[1] = 4200.
  all_crval1[2] = 4450.
  all_crval1[3] = 4675.
  all_crval1[4] = 5100.
  all_crval1[5] = 10^3.74818803
  all_crval1[6] = 10^3.79448805
  all_crval1[7] = 7000.
  all_crval1[8] = 8000.
  all_crval1[9] = 9300.
  for qq=0L,9 do begin
      mn = min(abs(all_crval1[qq]-tot_wave),imn)
      all_crval1[qq] = alog10(tot_wave[imn])
  endfor
  all_coll = dblarr(10,2)
  all_coll[0,*] = [4100., 4250.]
  all_coll[1,*] = [4300., 4490.]
  all_coll[2,*] = [4600., 4850.]
  all_coll[3,*] = [5100., 5450.]
  all_coll[4,*] = [5500., 5850.]
  all_coll[5,*] = [5900., 6290.]
  all_coll[6,*] = [6600., 6850.]
  all_coll[7,*] = [7390., 7650.]
  all_coll[8,*] = [8650., 8900.]
  all_coll[9,*] = [9500., 9900.]
  all_mnxwv = dblarr(10,2)
  all_mnxwv[0,*] = [3900., 4380.]
  all_mnxwv[1,*] = [4200., 4700.]
  all_mnxwv[2,*] = [4450., 5060.]
  all_mnxwv[3,*] = [4660., 5500.]
  all_mnxwv[4,*] = [5055., 5980.]
  all_mnxwv[5,*] = [5618., 6572.]
  all_mnxwv[6,*] = [6230., 7300.]
  all_mnxwv[7,*] = [7000., 8210.]
  all_mnxwv[8,*] = [8000., 9380.]
  all_mnxwv[9,*] = [9330., 10200.]

; Open Slit file
  if not keyword_set( SEDG_FIL ) then $
    sedg_fil = esi_getfil('sedg_fil', SLIT = esi[indx[0]].slit, cbin = cbin, $
                          rbin = rbin, /name)
  if x_chkfil(sedg_fil+'*') EQ 0 then begin
      print, 'esi_echextobj: Slit edge file doesnt exist: ', sedg_fil
      return
  endif
  print, 'esi_echechextobj: Grabbing slit edges from: ', sedg_fil
  slit_edg = xmrdfits(sedg_fil, /silent)
  dim1 = size(slit_edg, /dim)
  norders = dim1[1] 
  tset_slits = xmrdfits(sedg_fil, 1)
  nx = tset_slits[0].dims[0]
  ny = tset_slits[0].dims[1]
  yarr = replicate(1.0, nx) # findgen(ny)
  max_hwidth = fltarr(norders)
  FOR kk = 0L, norders-1L DO $
    max_hwidth[kk] = max(slit_edg[*, kk, 1]-slit_edg[*, kk, 0])
  rnd_edg = round(slit_edg)

  ;; This flag prevents updating of noise model from sky image in
  ;; long_localskysub.pro. Give this a more intuitive name
  DEFRINGE_FLAG = lonarr(norders) 

  
  ;; create some quantitites we will need
  ordermask = esi_ordermask(tset_slits)
  ximg = long_slits2x(tset_slits, edgmask = edgmask, nslit = nslit)
  ;; create a fake header with the read noise
  mkhdr, scihdr, 4, [nx, ny]
  sxaddpar, scihdr, 'RDNOISE', esi[indx[0]].READNO
  bsp = 0.6d
  nccd = 1
  yarr = 0

  window, 0
  window, 1
  wset, 0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Loop over exposures
  FOR q = 0L, n_elements(exp)-1 do begin
      print, 'esi_echextobj: Reading files...'
      ;; Set the default gaussian sigma to be 3.0 pixels (about
      ;; 0.8" seeing) for just the first order (6). This will be
      ;; adjusted as we extract the orders
      DEFGAU = (3.0d < ABSMAXGAU) > ABSMINGAU 
      ;;;;;;;;;;;;;;
      ;; Open Obj file
      objfil = esi[indx[exp[q]]].obj_fil
      if x_chkfil(objfil+'*') EQ 0 then begin
          print, 'esi_echextobj: No Obj file! ', objfil, ' Skipping...'
          continue
       endif
      ;; I don't see the point of these named structures JFH 01-13-2016
      ;;objstr = xmrdfits(objfil, 1, STRUCTYP = 'dblsobjstrct', /silent)
      objstr = xmrdfits(objfil, 1, /silent)
      IF KEYWORD_SET(OBJ_NM) THEN BEGIN
          nobj = n_elements(obj_nm)
          alpha_redux = obj_nm
      ENDIF ELSE BEGIN
          nobj = n_elements(objstr)/norders
          alpha_redux = alphabet[lindgen(nobj)]
      ENDELSE
      ;;;;;;;;;;;;;;
      ;; SKY SUB Fil 
      imgfil = esi_getfil('fin_fil', subfil = esi[indx[exp[q]]].img_root, /name)
      if x_chkfil(imgfil+'*') EQ 0 then begin
          print, 'esi_echextobj: No Image file!  Returning...'
          return
      endif
      print, 'esi_echextobj: Image file -- ', imgfil
      head = xheadfits(imgfil)
      sciimg = xmrdfits(imgfil, 0, /silent)
      var = xmrdfits(imgfil, 1, /silent)
      skyimage = xmrdfits(imgfil, 2, /silent)
      skyimage_orig = skyimage
      piximg = xmrdfits(imgfil, 3, /silent)
      ;; the 3* is because otherwise var=-1 causes division by zero. 
      sciivar = (var GT 0.0)/(var + 3*(var LE 0.0))
      outmask = fltarr(nx, ny)
      modelivar = sciivar*(ordermask GT 0.0)

      objimage = fltarr(nx, ny)
      prof_img = fltarr(nx, ny)
      ;; Read ARC
      if keyword_set(AIMG) then arc_fil = AIMG else begin
         arc_fil = strtrim(esi[indx[exp[q]]].arc_fil, 2)
         print, 'esi_echextobj: Arc -- ', arc_fil
         if x_chkfil(arc_fil+'*') EQ 0 then begin
            print, 'esi_echextobj: No Arc file!  Returning...', arc_fil
            return
         endif
      endelse
      waveimg = xmrdfits(arc_fil, /silent) 
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;;  HELIO correction  ;; -1 fixed on 3/25/09 JXP
      ;; RA and DEC need to be deg  12/15/2009 JXP
      x_radec, esi[indx[exp[q]]].ra, esi[indx[exp[q]]].dec, radeg, decdeg 
      helio = (-1.)*x_keckhelio(radeg, decdeg $
                          , esi[indx[exp[q]]].equinox $
                          , jd = esi[indx[exp[q]]].date)
      hel_corr = sqrt( (1.d + helio/299792.458) / (1.d - helio/299792.458) )
      waveimg = waveimg*hel_corr
      sxaddpar, head, 'HELIO', helio
      print, 'esi_echextobj: Helio correction applied -- ', helio, hel_corr
      ;; This is mainly used for the defringing where negative objects
      ;; are masked in the sky-subtraction. 
      objmask = lonarr(nx, ny)
      objmask[where(ordermask)] = 1L
      ndefr = n_elements(DEFRINGE)
      IF ndefr GT 0 THEN BEGIN
         ;; Defringe option assumes code has been run once, and so a
         ;; modelvar image exist. Use this as the sciivar, since this
         ;; has no pixel noise in it and is thus more robust.
         modelvar = xmrdfits(imgfil, 4)
         outmask  = xmrdfits(imgfil, 6)
         ;; Since we won't be updating the noise model on the
         ;; defringe pass, just set both sciivar and modelvar to the
         ;; noiseless noise model. 
         sciivar = outmask*(modelvar GT 0.0)/(modelvar + 3*(modelvar LE 0.0))
         modelivar = sciivar
         var =  outmask*modelvar
         ;; Now process the second negative image
         imgfil2 = esi_getfil('fin_fil', subfil = esi[indx[exp[q]]].img_root_fringe, /name)
         objfil2 = esi[indx[exp[q]]].obj_fil_fringe
         objstr2 = xmrdfits(objfil2, 1, /silent)
         nobj2 = n_elements(objstr2)/norders

         sciimg2 = xmrdfits(imgfil2, 0, /silent)
         modelvar2 = xmrdfits(imgfil2, 4)
         outmask2  = xmrdfits(imgfil2, 6)
         skyimage2 = xmrdfits(imgfil2, 2, /silent)
         ;;piximg2 = xmrdfits(imgfil2, 3, /silent)
         ;; the 3* is because otherwise var=-1 causes division by
         ;; zero.
         sciivar2 = outmask2*(modelvar2 GT 0.0)/(modelvar2 + 3*(modelvar2 LE 0.0))
         modelivar2 = sciivar2
         sciimg_diff = sciimg - sciimg2
         skyimage_diff = skyimage - skyimage2
         ;; ???? Think more carefully about masking here and CR rejection (later)
         var_diff = modelvar + modelvar2
         mask_diff = (outmask GT 0.0) AND (outmask2 GT 0.0)
         sciivar_diff = float(mask_diff GT 0.0)/(var_diff + 3*(var_diff LE 0.0))
         modelivar_diff = sciivar_diff
         FOR iorder = 0L, ndefr-1L DO BEGIN
            qq = DEFRINGE[iorder]
            thismask = (ordermask EQ (15-qq)) AND waveimg GT 0.0
            ithis = WHERE(thismask)
            sciimg[ithis] = sciimg_diff[ithis]
            skyimage[ithis] = skyimage_diff[ithis]
            var[ithis] = var_diff[ithis]
            sciivar[ithis] = sciivar_diff[ithis]
            modelivar[ithis] = modelivar_diff[ithis]
         ENDFOR
         ;; Create an object mask to mask negative traces in difference image
         trcmask = lonarr(ny, norders) + 1L
         ymin_ord15 = 1500
         ymax_ord15 = 3700
         ymax_ord06 = 2000.0 ;; 2170 ;; JFH Should this bet set to 2000 to fix order 6??
         ;; Mask the bad tracee (and off edge of CCD)
         trcmask = lonarr(ny, norders) + 1L
         trcmask[0:ymin_ord15, 0] = 0
         trcmask[ymax_ord15:*, 0] = 0
         trcmask[ymax_ord06:*, 9] = 0
         IF NOT KEYWORD_SET(NOFRINGEMASK) THEN BEGIN 
            FOR iobj = 0L, nobj2-1L DO BEGIN
               FOR iorder = 0L, ndefr-1L DO BEGIN
                  qq = DEFRINGE[iorder]
                  sci1 = where(objstr2.obj_id EQ alpha_redux[iobj] $
                               AND objstr2.order EQ qq)
                  objcen   = objstr2[sci1].trace[0:ny-1L]
                  obj_fwhm = (objstr2[sci1].FWHM > MINFWHM[qq]) < MAXFWHM[qq]
                  ;; extract asymmetric boxcars for S/N calculation
                  left_edge   = (objcen - FRINGEMASK/2.0*obj_fwhm)
                  right_edge  = (objcen + FRINGEMASK/2.0*obj_fwhm)
                  FOR j = 0, ny-1L DO BEGIN
                     IF trcmask[j, qq] EQ 1 THEN BEGIN 
                        xmin_obj = floor(left_edge[j]) >  0
                        xmax_obj = ceil(right_edge[j]) < (nx-1L)
                        objmask[xmin_obj:xmax_obj, j] = 0L
                     ENDIF
                  ENDFOR
               ENDFOR
            ENDFOR
         ENDIF
         DEFRINGE_FLAG[DEFRINGE] = 1 ;; This prevents updating of noise model from sky image
         outfil = esi_getfil('fringe_fil', subfil = esi[indx[exp[q]]].img_root, /name)
         obj_outfil = 'Extract/FringeObj_'+esi[indx[exp[q]]].img_root
      ENDIF ELSE BEGIN
         outfil = imgfil
         obj_outfil = objfil
      ENDELSE
      ;; Calculate S2N of each object for each order 
      ordr_s2n = fltarr(nord, nobj)
      FOR iobj = 0L, nobj-1L DO BEGIN
          FOR iorder = 0L, nord-1L DO BEGIN
              qq = ordrs[iorder]
              sci1 = where(objstr.obj_id EQ alpha_redux[iobj] $
                           AND objstr.order EQ qq)
              objcen   = objstr[sci1].trace[0:ny-1L]  
              ;; rename the objects 1 and 2 for each slit to simplify QA
              ;; use real-order number for QA. 
              objstr[sci1].OBJID = iobj + 1L
              objstr[sci1].SLITID = 15L-qq
              ;; extract asymmetric boxcars for S/N calculation
              left_edge   = (objcen - box_rad/cbin) > (rnd_edg[*, qq, 0]+LSLITE)
              right_edge  = (objcen + box_rad/cbin) < (rnd_edg[*, qq, 1]-RSLITE)
              ;; mask
              snmask = ordermask EQ (15-qq) AND var GT 0.0 AND waveimg GT 0.0
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
              indsp = WHERE(wave GT 2900.0 AND wave LT 1.3d4 $
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
              ordr_s2n[iorder, iobj] = sqrt(abs(med_sn2))
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
      forprint, 15L-ordrs[srt_s2n], ordr_s2n[srt_s2n, 0] $
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
          thismask = (ordermask EQ (15-qq)) AND waveimg GT 0.0
          ;; If S/N is lower than LOWSNR and at least two FWHM have 
          ;; been determined then fit FWHM vs order and use for profile
          other_orders = where(fwhm_here GT 0 AND fwhm_was_fit EQ 0, nother)
          other_fit    = where(fwhm_here GT 0 AND fwhm_was_fit EQ 1, nfit)
          FOR isrt = 0L, nobj-1L DO BEGIN
             iobj = srt_obj[isrt]
             IF ordr_s2n[iorder, iobj] LE LOWSNR AND nother GE 2 THEN BEGIN
                 ;; If this is the brightest object, then we have to 
                 ;; extrapolate the FWHM from a fit. 
                 IF isrt EQ 0 THEN BEGIN
                    splog, '*************************************************'
                    splog, 'Doing linear fit for FWHM of object=' $
                           + alpha_redux[iobj] + ' in order# '$
                           + strcompress(string(15L-qq), /rem)
                    print,  '   order#   FWHM' 
                    forprint, 15L-ordrs[other_orders] $
                              , fwhm_here[other_orders] $
                              , format = '(I6,F10.2)', textout = 1
                    fwhm_coeffs = ladfit(15L-ordrs[other_orders] $
                                         , fwhm_here[other_orders])
                    def_fwhm = poly(15L-qq, fwhm_coeffs) 
                    def_fwhm = ((def_fwhm[0] > MINFWHM[qq]) < MAXFWHM[qq])
                    fwhm_was_fit[iorder] = 1
                    splog, 'FIT: ' + string(15L-qq, format = '(I4)') $
                           + string(def_fwhm, format = '(F10.2)')
                    splog, '*************************************************'
                    wset, 1
                    plot, 15L-ordrs[other_orders], fwhm_here[other_orders]*plate_scale[other_orders] $
                          , psym = 4, symsize = 2, xr = [5.0, 16.0] $
                          , yrang = [0.7*min([fwhm_here*plate_scale, def_fwhm]*plate_scale[qq]) $
                                     , 1.4*max([fwhm_here*plate_scale, def_fwhm*plate_scale[qq]])] $
                          , xtitle = 'Order Number', ytitle = 'FWHM (arcsec)' $
                          , /xsty, /ysty
                    clr = getcolor(/load)
                    IF nfit GT 0 THEN oplot, 15L-ordrs[other_fit] $
                                             , fwhm_here[other_fit]*plate_scale[other_fit] $
                                             , psym = 6, symsize = 2 $
                                             , col = clr.cyan
                    oplot, [15L-qq], [def_fwhm]*plate_scale[qq], psym = 5, symsize = 2 $
                           , col = clr.green
                    oplot, 15L-ordrs, (poly(15L-ordrs, fwhm_coeffs) > 2.0/cbin)*plate_scale $
                           , col = clr.red
                    wset, 0
                    sci1 = where(objstr.obj_id EQ alpha_redux[iobj] $
                                 AND objstr.order EQ qq)
                    objstr[sci1].FWHM = def_fwhm
                 ENDIF ELSE IF ordr_s2n[iorder, iobj] LE LOWSNR AND nother EQ 1 THEN BEGIN
                    ;; If just one order is good, then use that FWHM throughout
                    def_fwhm = fwhm_here[other_orders]
                    sci1 = where(objstr.obj_id EQ alpha_redux[iobj] $
                                 AND objstr.order EQ qq)
                    objstr[sci1].FWHM = def_fwhm  > 2.0/cbin 
                 ENDIF ELSE BEGIN 
                    sci1 = where(objstr.obj_id EQ alpha_redux[iobj] $
                                 AND objstr.order EQ qq)
                    sci1a = where(objstr.obj_id EQ alpha_redux[ibright] $
                                  AND objstr.order EQ qq)
                    objstr[sci1].FWHM =  objstr[sci1a].FWHM  > 2.0/cbin 
                 ENDELSE
              ENDIF
           ENDFOR
          IF NOT KEYWORD_SET(BOXCAR) THEN BEGIN
             extract_struct = long_localskysub(sciimg, sciivar*(objmask GT 0), skyimage $
                                               , piximg, waveimg, ximg, objstr $
                                               , thismask $
                                               , slit_edg[*, qq, 0] $
                                               , slit_edg[*, qq, 1] $
                                               , edgmask, bsp = bsp $
                                                , /SKYSAMPLE $
                                                , niter = niter $
                                                , objimage = objimage $
                                                , modelivar = modelivar $
                                                , varnoobj = varnoobj $
                                                , outmask = outmask $
                                                , indx = ordind, nccd = nccd $
                                                , box_rad = box_rad $
                                                , profile_struct =  $
                                                profile_struct $
                                                , scihdr = scihdr $
                                                , SN_GAUSS = SN_GAUSS $
                                                , CHK = LCHK, COADD_2D = DEFRINGE_FLAG[qq])
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
                                                     , piximg, waveimg, ximg, objstr $ 
                                                     , thismask $
                                                     , slit_edg[*, qq, 0] $
                                                     , slit_edg[*, qq, 1] $
                                                     , edgmask, bsp = bsp $
                                                     , /SKYSAMPLE $
                                                     , niter = niter $
                                                     , objimage = objimage $
                                                     , modelivar = modelivar $
                                                     , varnoobj = varnoobj $
                                                     , outmask = outmask $
                                                     , indx = ordind $
                                                     , nccd = nccd $
                                                     , box_rad = box_rad $
                                                     , profile_struct = $
                                                     profile_struct $
                                                     , scihdr = scihdr $
                                                     , SN_GAUSS = SN_GAUSS $
                                                     , CHK = LCHK, COADD_2D = DEFRINGE_FLAG[qq])
             ENDIF
          ENDIF
         ;; Now loop over the objects and extract each using profile
         FOR isrt = 0L, nobj-1L DO BEGIN
             iobj = srt_obj[isrt]
             print, '-----------------------------------------------'
             print, 'esi_echextobj: Extracting order ' $
                    , string(15L-qq, FORMAT = '(i3)') +  $
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
           IF NOT KEYWORD_SET(BOXCAR) THEN tnovar = transpose(varnoobj[mincol:maxcol, *]* $
                                                              outmask[mincol:maxcol, *])
           tsky = transpose(skyimage[mincol:maxcol, *]* $
                            thismask[mincol:maxcol, *])
           tsky_orig = transpose(skyimage_orig[mincol:maxcol, *]* $
                                 thismask[mincol:maxcol, *])
            IF NOT KEYWORD_SET(BOXCAR) THEN BEGIN ;; OPTIMAL EXTRACTION
               tprof = $
                  transpose(profile_struct[iobj].PROFILE[mincol:maxcol, *]* $
                            thismask[mincol:maxcol, *]) 
               prof_img = prof_img + profile_struct[iobj].PROFILE
               ;; 
               ;; Extract
               ;; 
               ;; This optimal extraction routine re-synthesizes
               ;; variance. For Defringing with differenced frames
               ;; that noise won't be correct but our noise was
               ;; correct. ???
               ;; The optimal extractino routine now uses my input
               ;; variance image and does not try to re-calculate. 
               ;; Think about if Xavier's code was doing
               ;; anything deep (I think not). 
               IF DEFRINGE_FLAG[qq] EQ 1 THEN tsky_opt = tsky_orig ELSE tsky_opt = tsky
               ;; For defrining, the tsky is the sky residual imag,
               ;; for that reason, it needs the actual original sky value in tsky_opt.
               x_extobjopt, timg-tsky, tarc, tvar, tsky_opt, PROFILE = tprof $
                            , NOVAR = tnovar $
                            , [objstr[sci2].xcen, objstr[sci2].ycen-mincol] $
                            , fin_spec, COLLMNX = all_coll[qq, *] $
                            , CRVAL1 = all_crval1[qq], CDELT = cdelt $
                            , NPIX = round(5000./rbin) $
                            , TOT_TRC = objstr[sci2].trace[0:ny-1L]-mincol $
                            , /REBINC, READNO = esi[indx[exp[q]]].readno $
                            , REDBLUE = 0L, PROF_LIM = PROF_LIM $
                            , CHK = chk, DEBUG = debug, /USE_INPUT_VAR
               ;; PROF_LIM added by JFH 06/08. This explicitly masks
               ;; any region for which the fraction of good pixels 
               ;; on the 2-d object profile frac_mask < PROF_LIM =
               ;; 0.40
            ENDIF ELSE BEGIN ;; BOXCAR EXTRACTION
                tvar = transpose(var[mincol:maxcol, *])
                x_extobjbox, timg-tsky, tarc, VAR = tvar $
                             , [objstr[sci2].xcen, objstr[sci2].ycen-mincol] $
                             , fin_spec, WVMNX = all_mnxwv[qq, *] $
                             , APER = aper, DEBUG = debug $
                             , COLLMNX = all_coll[qq, *] $
                             , CRVAL1 = all_crval1[qq], CDELT = cdelt $
                             , NPIX = 5000L, REJ_CR = REJ_CR $
                             , TOT_TRC = objstr[sci2].trace[0:ny-1L]-mincol $
                             , REJSIG = rejsig, /REBINC, RADIUS = box_rad $
                             , NOCRREJ = NOCRREJ, SKY=tsky $
                             , BKAPER = svaper, CHK = chk 
                ;; Save the aperture
              svaper=fin_spec.aper
           ENDELSE
            ;;spcchk = 1
            if keyword_set(SPCCHK) AND fin_spec.npix NE 0 then begin
               x_splot, fin_spec.wv, fin_spec.fx $
                        , YTWO = sqrt(fin_spec.var > 0.), /block
           endif
            ;; Write to structure
            if fin_spec.npix NE 0 then begin
               objstr[sci2].npix = fin_spec.npix
               objstr[sci2].wave[0:fin_spec.npix-1] = fin_spec.wv
               objstr[sci2].fx[0:fin_spec.npix-1] = fin_spec.fx
               objstr[sci2].var[0:fin_spec.npix-1] = fin_spec.var
               objstr[sci2].novar[0:fin_spec.npix-1] = fin_spec.novar
               objstr[sci2].sky[0:fin_spec.npix-1] = fin_spec.sky
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
      profil = esi_getfil('fin_fil', subfil = esi[indx[exp[q]]].img_root, /name)
      print, systime()
      ;; Write Spectra
      print, 'esi_echextobj: Output spectrum in -- ', obj_outfil
      mwrfits, objstr, obj_outfil, /create, /silent
      spawn, 'gzip -f '+obj_outfil
      ;; Write reduction outputs
      print, 'esi_echextobj: Writing output to: ', outfil
      IF keyword_set(CHK) THEN BEGIN
          print, 'esi_echextobj:  Displaying sky subtracted image'
          xatv, (sciimg - skyimage)*float(ordermask GT 0.0) $
                , WVIMG = waveimg, /block, min = -50, max = 200
          print, 'esi_echextobj:  Displaying residual image'
          xatv, (sciimg - skyimage - objimage)*sqrt(modelivar)* $
                float(ordermask GT 0.0), WVIMG = waveimg, /block $
                , min = -6.0, max = 6.0
      ENDIF
      mwrfits, float(sciimg), outfil, head, /create, /silent
      mwrfits, float(var), outfil, /silent
      mwrfits, float(skyimage), outfil, /silent
      mwrfits, float(piximg), outfil, /silent
      mwrfits, float(modelvar), outfil, /silent
      mwrfits, float(objimage), outfil, /silent
      mwrfits, float(outmask), outfil, /silent
      mwrfits, float(prof_img), outfil, /silent
      ;; Compress
      print, 'esi_echextobj: Compressing...'
      spawn, 'gzip -f '+outfil
   endfor
  
  
;  DONE
  print, 'esi_echextobj: All done! '
  return
end
