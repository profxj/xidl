;+ 
; NAME:
; esi_echfltsct   
;     Version 1.1
;
; PURPOSE:
;    Subtracts a scattered light model from the flat
;    Normalizes the Image and Outputs to same file with
;    updated headers
;
; CALLING SEQUENCE:
;   
;  esi_echfltsct, esi, /IFLAT, GAIN_RTO=, /CHK
;
; INPUTS:
;   esi     -  ESI structure
;   slit    -  Slit width  (0.5, 0.75, 1.0)
;
; RETURNS:
;
; OUTPUTS:
;  Image with scattered light removed (e.g. Flats/FlatECH##_D[I].fits)
;
; OPTIONAL KEYWORDS:
;   /IFLAT    - Use internal flats 
;   GAIN_RTO= - Value for gain ratio for two amp mode (default:
;      calculate from image)
;   /CHK      - Display final image and a few other steps
;   MINWID=   - Minimum width between orders for measuring scattered light
;   /PINH     - Indicates pinhole slit (used in early 2000)
;   /SCAT     - Do scattered light subtraction (default is not to)
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Only set for 1x1 binning -- Many 'hard wired' numbers
;
; EXAMPLES:
;   esi_echfltsct, esi, 0.75, /CHK
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   11-Aug-2002 Written by JXP
;   01-Feb-2003 Polished (JXP)
;   17-Feb-2004 Added kludge for order #2
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echfltsct, esi, slit, IFLAT=iflat, GAIN_RTO=gain_rto, CHK=chk $
                   , MINWID = minwid, PINH = pinh, CBIN = cbin, RBIN = rbin $
                   , KLUDGEOR = kludgeor, IFU = ifu, BSPLINE = BSPLINE $
                   , NRMMED = NRMMED, SCAT = SCAT, TWIFLAT = TWIFLAT $
                   , TWEAK = TWEAK , OUTFIL=outfil, ARC_FILE=arc_file $
                   , ORDER_VEC=order_vec

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echfltsct, esi, slit, /IFLAT, /CHK, GAIN_RTO=, MINWID= '
      print, '       /PINH [v1.1]'
      return
  endif 
  
;  Optional Keywords

  if keyword_set(IFLAT) THEN ftype = 'I' $
  ELSE IF KEYWORD_SET(TWIFLAT) THEN ftype = 'T' $
  ELSE ftype = 'D'
  if not keyword_set(SMSHROW) then smshrow = 2150L
  if not keyword_set( FITFIL ) then fitfil = 'Maps/hole_fit.idl'
  if not keyword_set( MINWID ) then minwid = 6.
  if not keyword_set( MED_SMTH ) then med_smth = 5L
  if not keyword_set( XSEP ) then xsep = 2L
  if not keyword_set( CBIN ) then cbin = 1
  if not keyword_set( RBIN ) then rbin = 1
  if not keyword_set( NRMMED ) then NRMMED = 20L;/cbin
  IF NOT KEYWORD_SET(SPECSAMP) THEN BEGIN 
      IF KEYWORD_SET(TWIFLAT) THEN SPECSAMP = 2.5 $
     ELSE SPECSAMP = 5L 
  ENDIF
  IF NOT KEYWORD_SET(SPATSAMP) THEN SPATSAMP = 0.02D
    
  ;; Open Flat
  IF ftype EQ 'I' then mode = 0 $
  ELSE IF ftype EQ 'D' THEN mode = 1 $
  ELSE IF ftype EQ 'T' THEN mode = 2
  flatfil = esi_getfil('echflat_fil', mode, SLIT=slit, $
                       cbin=cbin, rbin=rbin, /name)
          
  print, 'esi_echfltsct: Opening flat ', flatfil
  if x_chkfil(flatfil+'*') EQ 0 then begin
      print, 'esi_echfltsct: File ', flatfil, ' does not exist.  Returning!'
      return
  endif
  flat = xmrdfits(flatfil, 0, head, /silent)
  sz_flat = size(flat, /dimensions)
  nx = sz_flat[0]
  ny = sz_flat[1] 
  ;; Open Slit file
  slit_file = esi_getfil('sedg_fil', SLIT = slit, cbin = cbin, rbin = rbin $
                         , /name)
  slit_edg = xmrdfits(slit_file, 0, /silent)
  tset_slits = xmrdfits(slit_file, 1)
  dslit = tset_slits[1].COEFF[0, *]-tset_slits[0].COEFF[0, *]
  ;; Open wavelength map
  if not keyword_set(ARC_FILE) then $
    arc_file = esi_getfil('arc_img', SLIT = slit, cbin = cbin, rbin = rbin, /name)
  waveimg = xmrdfits(arc_file, 0)
  ;; create auxilary images needed for flat field fits
  ordermask = esi_echordermask(tset_slits)
  ximg = long_slits2x(tset_slits, edgmask = edgmask, TOL_EDG = 3L $
                      , nslit = norders)

  ;; Truncate orders (usually the reddest)
  if not keyword_set(ORDER_VEC) then begin
     order_vec = lindgen(10) + 6
  endif else norders=n_elements(order_vec)

  xarr = findgen(nx) # replicate(1.0, ny)
  yarr = findgen(ny)## replicate(1.0, nx)
  ;; dispersion per order
  print, 'esi_echfltsct: Normalizing the Gain'

  ;; these are the spatial regions used for the illumination function
  ;spat_regions = fltarr(10, 2)
  ;buf = 4L
  ;spat_regions[0, *] = [buf, 2170]
  ;spat_regions[1, *] = [buf, ny-buf]
  ;spat_regions[2, *] = [buf, ny-buf]
  ;spat_regions[3, *] = [buf, ny-buf]
  ;spat_regions[4, *] = [240, ny-buf]
  ;spat_regions[5, *] = [430, ny-buf]
  ;spat_regions[6, *] = [785, ny-buf]
  ;; these are pixels where flat is dominated by scattered light
  ;; or where the flat is bad and set to unity
  scat1 = (ordermask EQ 15 AND waveimg GT 4350) $
    OR  (ordermask EQ 14 AND waveimg LT 4230) $
    OR  (ordermask EQ 13 AND waveimg LT 4500) $
    OR  (ordermask EQ 12 AND waveimg LT 4840) $
    OR  (ordermask EQ 11 AND waveimg LT 5200) $
    OR  (ordermask EQ 10 AND waveimg LT 5670) $
    OR  (ordermask EQ 6  AND waveimg GT 10200)


  IF NOT KEYWORD_SET(HOT_THRESH) THEN HOT_THRESH = 0.2
  ;; hot spot masking for flats
  sky_sec = WHERE(xarr GT 3 $
                  AND yarr GE 3800 $
                  AND yarr LE 3860 $
                  AND NOT EDGMASK $
                  AND ordermask EQ 15, nsec)
  IF nsec GT 0 THEN BEGIN
      djs_iterstat, flat[sky_sec], median = avg_hot_sky, sigrej = 3.0
      img_bad = WHERE(yarr GE 3860 AND yarr LE 3950 AND $
                      flat GE (avg_hot_sky + HOT_THRESH*sqrt(avg_hot_sky)) $
                      AND ordermask EQ 15, nhot)
      flat[img_bad] = -1.0
  ENDIF
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Normalize Gain
  case cbin of
      1: midx = 1023L
      2: midx = 512L
      else: stop
  endcase

  IF NOT KEYWORD_SET(TWIFLAT) THEN BEGIN
     if not keyword_set(GAIN_RTO) THEN BEGIN
        ;; Centers
        cen_4 = (slit_edg[*, 4, 0]+slit_edg[*, 4, 1])/2.
        ;; Running median
        lft4 = round(cen_4)-7L
        rgt4 = round(cen_4)+7L
        arr4 = fltarr(sz_flat[1], 15)
        for i = 0L, sz_flat[1]-1 do arr4[i, *] = flat[lft4[i]:rgt4[i], i]
        ;; Smash
        smsh4 = djs_median(arr4, 2)
        ;; Fit RHS
        gainr_fit = x_setfitstrct(NITER = 1L, NORD = 2L, FLGREJ = 1L $
                                  , HSIG = 3., LSIG = 3. $
                                  , FUNC = 'POLY')
        gainl_fit = gainr_fit
        ;; Binning
        if rbin NE 1 then stop
        ;; Fit
        fitrhs = x_fitrej(2700L+findgen(221), smsh4[2700L:2920L] $
                          , FITSTR = gainr_fit)
        fitlhs = x_fitrej(3040L+findgen(361), smsh4[3040L:3400L] $
                          , FITSTR = gainl_fit)
        mn = min(abs(cen_4[1500L:sz_flat[1]-1]-(midx-0.5)), imn)
        ;; CORRECT RHS
        val_rhs = x_calcfit(float(imn+1500L), FITSTR = gainr_fit)
        val_lhs = x_calcfit(float(imn+1500L), FITSTR = gainl_fit)
        gain_rto = val_lhs/val_rhs
        print, 'esi_echfltsct: Gain Ratio LHS/RHS = ', gain_rto
     endif
     sxaddpar, head, 'GAINFIX', gain_rto
  ENDIF ELSE BEGIN
     ;; Don't measure the gain ratio from twilight flats
     flatfil = $
        esi_getfil('finflat_fil', SLIT = slit, cbin = cbin, rbin = rbin, /name)
     IF file_test(flatfil+'*') THEN BEGIN
        headflt = xheadfits(flatfil, /silent)
        gain_rto = sxpar(headflt, 'GAINFIX', count = nval)
     ENDIF ELSE BEGIN
        print, '    *********************************'
        print, '    Flat file does not exist which is'
        print, '    needed for the gain ratio.'
        print, '    Please process the dome or internal'
        print, '    flats first with esi_echfltsct'
        print, '    *********************************'
        stop
     ENDELSE
  ENDELSE
  flat[midx:*, *] = flat[midx:*, *]*gain_rto
  ;; Now create a model of the scattered light
  ;; This is commented out since the mage_scattlight.pro file no
  ;; longer exists. commented out by MR 3/4/2013
  newflt = flat
;  IF NOT KEYWORD_SET(SCAT) THEN newflt = flat $
;  ELSE BEGIN
;      ;; need to convince myself scattered light correction is 
;      ;; benefitting the flats. 
;      img_scatt = mage_scattlight(flat, tset_slits, scatt_shift = 3L)
;      ;; set scattered light to zero for order 15
;      if keyword_set(CHK) then xatv, img_scatt, /block
;      img_scatt[where(ordermask EQ 15 OR ordermask EQ 14)] = 0.0D
;      newflt = flat - img_scatt 
;  ENDELSE
  sxaddpar, head, 'SCATTER', 'T'
  nrm_flat = fltarr(sz_flat[0], sz_flat[1])
  illum_flat = fltarr(sz_flat[0], sz_flat[1]) + 1.0
  print, 'esi_echfltsct: Normalizing.. ' 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NORMALIZE
  FOR iorder = 0L, norders-1L DO BEGIN
;      iorder = norders -2L
      inorder = where(ordermask EQ order_vec[iorder], nord)
      fitpix = where(ordermask  EQ order_vec[iorder]  $
                     AND finite(flat)  $
                     AND flat GT 1.0 $
                     AND abs(flat) LE 6.0d4 $
                     AND waveimg GT 0.0 $
                     AND xarr GT 3 $
                     AND xarr LE ny-2L $
                     AND NOT EDGMASK, nfit)
      ;;ximg GT 0.05  $
      ;;AND ximg LT 0.95, nfit)
      psort = sort(waveimg[fitpix])
      med_spec_width = round(double(nfit)*double(specsamp)/double(ny))
      normspec1 = djs_median(newflt[fitpix[psort]], Width = med_spec_width $
                             , bound = 'reflect')
      ;; now smooth with a gausian 
      sig_res_spec = med_spec_width/5L
      nhalf =  long(sig_res_spec)*4L
      xkern = dindgen(2*nhalf+1)-nhalf
      kernel = gauss1(xkern, [0.0, sig_res_spec, 1.0])
      normspec = convol(normspec1, kernel, /edge_truncate)
      everyn = round(double(med_spec_width)/20.0)
      ;fullbkpt = (waveimg[fitpix[psort]])[lindgen(nfit/everyn)*everyn]
      ;fullbkpt = bspline_bkpts(waveimg[fitpix[psort]], nord = 4 $
      ;                         , nbkpts = round(double(nfit)/everyn), /silent)
      sset_spec = bspline_iterfit(waveimg[fitpix[psort]], (normspec >  0.0) $
                                  , upper = 3, lower = 3, yfit = yfit $
                                  , nbkpts = round(double(nfit)/everyn) $
                                  , maxiter = 20, maxrej = 10 $
                                  , requiren = 3)
      nrm_flat[inorder] = $
        newflt[inorder]/bspline_valu(waveimg[inorder], sset_spec)
      ;; now fit out average illuimination function across the order
      fitpix2 = where(ordermask EQ order_vec[iorder]  $
                      AND finite(flat)  $
                      AND flat GT 1.0 $
                      AND abs(flat) LE 6.0d4 $
                      AND waveimg GT 0.0 $
                      AND xarr GT 3 $
                      AND xarr LE ny-2L $
                      AND nrm_flat GT 0.0 AND nrm_flat LE 3.0 $
                      AND randomu(seed, nx, ny) GT 0.7 $
                      AND NOT scat1, nfit2) ;; exclude scattered light regions
      ;; randomly downsample for speed
      psort2 = sort(ximg[fitpix2])
      med_spat_width = round(double(nfit2)*double(spatsamp))
      normspat1 = djs_median(nrm_flat[fitpix2[psort2]], width = med_spat_width $
                             , bound = 'reflect')
      ;; now smooth with a gausian 
      sig_res_spat = med_spat_width/15L
      nhalf =  long(sig_res_spat)*4L
      xkern = dindgen(2*nhalf+1)-nhalf
      kernel = gauss1(xkern, [0.0, sig_res_spat, 1.0])
      normspat = convol(normspat1, kernel, /edge_truncate)
      statinds = WHERE(ximg[fitpix2[psort2]] GT 0.1 AND $
                       ximg[fitpix2[psort2]] LT 0.9)
      moms = moment(normspat[statinds])
      mean = moms[0]
      normspat = normspat/mean
      nrm_flat_orig = nrm_flat[fitpix2[psort2]]/mean
      sset_spat = bspline_iterfit(ximg[fitpix2[psort2]], normspat $
                                  , upper = 3, lower = 3, yfit = yfit2 $
                                  , bkspace = spatsamp/50.0d $
                                  , maxiter = 20, maxrej = 10)
      ;; don't allow super big corrections to protect edges
      normbsp = (bspline_valu(ximg[inorder], sset_spat)) < 3.0
      normbsp = normbsp > 0.2
      nrm_flat[inorder] = nrm_flat[inorder]/normbsp
      illum_flat[inorder] = normbsp
      leftind  = WHERE(xarr LT nx/2 AND ordermask EQ order_vec[iorder], nleft)
      rightind = WHERE(xarr GE nx/2 AND ordermask EQ order_vec[iorder], nright)
      ;; this fixes amplifier difference across order 10
      IF nleft GT 0 AND nright GT 0 THEN BEGIN
          leftgood = WHERE(xarr LT nx/2 $
                           AND ordermask EQ order_vec[iorder] $
                           AND waveimg GT 0.0 $
                           AND finite(flat) $ 
                           AND flat GT -100.0 $
                           AND finite(nrm_flat) $
                           AND nrm_flat GT 0.0 $
                           AND nrm_flat LT 3.0 $
                           AND abs(flat) LE 6.0d4 $
                           AND yarr GT 1200 $
                           AND NOT EDGMASK, nleft1)
          rightgood = WHERE(xarr GE nx/2 $
                            AND ordermask EQ order_vec[iorder] $
                            AND waveimg GT 0.0 $
                            AND finite(flat) $ 
                            AND flat GT -100.0 $
                            AND finite(nrm_flat) $
                            AND nrm_flat GT 0.0 $
                            AND nrm_flat LT 3.0 $
                            AND abs(flat) LE 6.0d4 $
                            AND yarr GT 1200.0 $
                            AND NOT EDGMASK, nright1)
          djs_iterstat, nrm_flat[leftgood], median = left_med
          if nright1 GT 0 then $
            djs_iterstat, nrm_flat[rightgood], median = right_med $
          else begin
              print, 'esi_echfltsct:  On the edge here..'
              print, 'esi_echfltsct:  Setting right_med = 1.'
              print, 'esi_echfltsct:  Continue as you wish'
              right_med = 1.
              stop
          endelse
          nrm_flat[leftind] = nrm_flat[leftind]*right_med/left_med
       ENDIF
      IF KEYWORD_SET(TWEAK) THEN BEGIN
         ;; now tweak the slit position so that the slits cut off no less 
         ;; than 80% from the boundary
         step = 0.005
         xleft_tweak = -1.0d
         inorm = WHERE(ximg[fitpix2[psort2]] GE 0.2 AND $
                       ximg[fitpix2[psort2]] LE 0.8)
         maxnorm = max(yfit2[inorm])
         FOR xleft = 0.5D, 0.0, -step DO BEGIN
            ynow = bspline_valu(xleft, sset_spat)
            IF ynow LE 0.9d*maxnorm THEN BEGIN 
               xleft_tweak = xleft
               BREAK
            ENDIF
         ENDFOR
         IF xleft_tweak NE -1.0 THEN BEGIN
            tset_slits[0].COEFF[0, iorder] = $
               tset_slits[0].COEFF[0, iorder] + xleft_tweak*dslit[iorder]
            splog, 'Tweaking position of LEFT boundary for order', 15L-iorder
         ENDIF
         xright_tweak = -1.0d
         FOR xright = 0.5D, 1.0, step DO BEGIN
            ynow = bspline_valu(xright, sset_spat)
            IF ynow LE 0.9d*maxnorm THEN BEGIN 
               xright_tweak = xright
               BREAK
            ENDIF
         ENDFOR
         IF xright_tweak NE -1.0 THEN BEGIN
             tset_slits[1].COEFF[0, iorder] = tset_slits[1].COEFF[0, iorder] - $
               (1.0d - xright_tweak)*dslit[iorder]
             splog, 'Tweaking position of RIGHT boundary for order', 15L-iorder
         ENDIF
         IF KEYWORD_SET(CHK) THEN BEGIN
            yr = [0.2, 1.3]
            clr=getcolor(/load)
            plot, ximg[fitpix2[psort2]], nrm_flat_orig $
                  , psym = 3, yr = yr, xr = [-0.1, 1.1], /xsty, /ysty $
                  , charthick = 2.0, thick = 3.0
            oplot, ximg[fitpix2[psort2]], yfit2, col = clr.red $
                   ,  thick = 2.0
            IF xleft_tweak NE -1.0 THEN $
               oplot, [xleft_tweak, xleft_tweak], yr, linestyle = 2 $
                      , col = clr.green $
                      , thick = 2.0
            IF xright_tweak NE -1.0 THEN $
               oplot, [xright_tweak, xright_tweak], yr, linestyle = 2 $
                      , col = clr.green $
                      , thick = 2.0
         ENDIF
      ENDIF
      ;;  xatv, double(ordermask EQ order_vec[iorder])*nrm_flat $
      ;;  , min = 0.9, max = 1.1, /block
   ENDFOR
  
  ;; Remake edgmask and ordermask for setting things to 1.0 if slit edges 
  ;; were tweaked
  If KEYWORD_SET(TWEAK) THEN BEGIN
      ximg = long_slits2x(tset_slits, edgmask = edgmask, TOL_EDG = 3L $
                          , nslit = norders)
      ordermask = esi_echordermask(tset_slits)
  ENDIF
  ;; these are edge pixels or pixels with no wavelengths
  unit1 = (waveimg LE 0.0 OR edgmask OR YARR LT 5 $
           OR xarr LT 3 OR xarr GT nx-2L $
           OR nrm_flat GT 3.0 OR nrm_flat LE 0.0 $
           OR ordermask EQ 0)
  unitind1 = WHERE(unit1 OR scat1, nunit1)
  IF nunit1 GT 0 THEN nrm_flat[unitind1] = 1.0d
  ;; make sure bad pixels are masked in final result
  badind1 = WHERE(flat LT -100.0, nbad1)
  IF nbad1 GT 0 THEN nrm_flat[badind1] = -1.0
  ;; these are edge pixels or pixels with no wavelengths
  unit2 = (waveimg LE 0.0 $
           OR xarr LT 3 OR xarr GT nx-2L $
           OR illum_flat GT 3.0 OR illum_flat LE 0.2 $
           OR ordermask EQ 0)
  unitind2 = WHERE(unit2, nunit2)
  IF nunit2 GT 0 THEN illum_flat[unitind2] = 1.0d
  ;; make sure bad pixels are masked in final result
  badind2 = WHERE(illum_flat LT -100.0, nbad2)
  IF nbad2 GT 0 THEN illum_flat[badind2] = -1.0

  if keyword_set(CHK) then xatv, nrm_flat, min = 0.9, max = 1.1, /block
  sxaddpar, head, 'NORM', 'T'
      
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; OUTPUT
;  outfil = 'Flats/FlatECH'+c_s+'N.fits'
  IF ftype EQ 'I' OR ftype EQ 'D' THEN $
     outfil = esi_getfil('finflat_fil', SLIT = slit $
                         , cbin = cbin, rbin = rbin, /name) $
  ELSE outfil = esi_getfil('twiflat_fil', SLIT = slit $
                           , cbin = cbin, rbin = rbin, /name) 
  print, 'esi_echfltsct: Creating Normalized Flat ', outfil
  mwrfits, nrm_flat, outfil, head, /create, /silent
  mwrfits, illum_flat, outfil, head, /silent
  ;; COMPRESS
  spawn, 'gzip -f '+outfil
  
  ;; update slit edges with tweaked values
  IF KEYWORD_SET(TWEAK) THEN BEGIN
     sedg_fil = $
        esi_getfil('sedg_fil', SLIT = slit, cbin = cbin, rbin = rbin, /name)
     print, 'esi_echfltsct: Writing slit edge info to: ', sedg_fil
     slit_edg = fltarr(sz_flat[1], 10, 2)
     traceset2xy, tset_slits[0], rows, left_edge
     traceset2xy, tset_slits[1], rows, right_edge
     slit_edg[*, *, 0] = left_edge
     slit_edg[*, *, 1] = right_edge
     mwrfits, slit_edg, sedg_fil, /create
     mwrfits, tset_slits, sedg_fil
     ;; CHK (Plot tweaked slit edges)
     if keyword_set(CHK) then begin
         tmp = illum_flat
         FOR i = 0, 9 do begin
             FOR j = 0, 1 do begin
                 rnd_trc2 = round(slit_edg[*, i, j])
                 trc_msk = rnd_trc2 + lindgen(sz_flat[1])*sz_flat[0]
                 tmp[trc_msk] = -10000
             ENDFOR
         ENDFOR
         xatv, tmp, /block, min = 0.7, max = 1.2
     endif
 ENDIF
 
 ;; Update header cards
; objstd = where(esi.mode EQ 2 AND esi.flg_anly NE 0 AND $
;                esi.slit EQ slit AND $
;                (strtrim(esi.type, 2) EQ 'STD' OR $
;                 strtrim(esi.type, 2) EQ 'OBJ' OR $ 
;                 strtrim(esi.type, 2) EQ 'ARC' ), nobj)
; IF nobj GT 0 THEN BEGIN
;     IF KEYWORD_SET(TWIFLAT) THEN esi[objstd].twiflat_fil = outfil $
;     ELSE esi[objstd].flat_fil = outfil
; ENDIF
 print, 'esi_echfltsct: All done!'
 
 return 
end



