FUNCTION SINFONI_WAVEIMG, skyfiles, tset_slits $
                          , darkfile = darkfile $
                          , pixflatfile = pixflatfile $
                          , illumflatfile = illumflatfile $
                          , scifile = scifile $
                          , QAFILE = QAFILE $
                          , piximg = piximg, pixset = pixset $
                          , MXSHIFT = MXSHIFT $
                          , SIGREJ = SIGREJ, WAVEFILE = WAVEFILE $
                          , BOX_RAD = BOX_RAD, CHK = CHK $
                          , CALIB = CALIB, NOSAVE = NOSAVE $
                          , outfile = outfile, savefile = savefile
  
  ;chk = 1
t0 = systime(1)

if not keyword_set(MXSHFT) then MXSHFT = 15L
if not keyword_set(SIGREJ) then sigrej = 2.
if (NOT keyword_set(box_rad)) then box_rad = 5L
arc_ncoeff = 3 ;; ???more curvature in arc lines use higher order?????

dims = tset_slits[0].DIMS
nx = dims[0]
ny = dims[1]
nsky = n_elements(skyfiles)
IF nsky EQ 1 THEN BEGIN
   sinfoni_proc, skyfiles[0], sky_bar, ivar_bar $
                 , hdr = hdr_arc, darkfile = darkfile, verbose = verbose $
                 , pixflatfile = pixflatfile, illumflatfile = illumflatfile
ENDIF ELSE BEGIN 
   sky_stk = fltarr(nx, ny, nsky)
   ivar_stk = fltarr(nx, ny, nsky)
   FOR ii = 0L, nsky-1L DO BEGIN 
      sinfoni_proc, skyfiles[ii], arcimg, arcivar $
                    , hdr = hdr_arc, darkfile = darkfile, verbose = verbose $
                    , pixflatfile = pixflatfile, illumflatfile = illumflatfile
      sky_stk[*, *, ii] = arcimg
      ivar_stk[*, *, ii] = arcivar
   ENDFOR
   msk_stk = ivar_stk GT 0.0
   sky_avs = djs_avsigclip(sky_stk, 3 $
                           , sigrej = sigrej, maxiter = maxiter $
                           , inmask = (msk_stk EQ 0) $
                           , outmask = outmask)
   weights = float(outmask EQ 0)
   w_sum = total(weights, 3)
   sky_bar = (w_sum GT 0)*total(weights*sky_stk, 3)/(w_sum + (w_sum EQ 0.0))
   var_stk = (ivar_stk GT 0.0)/(ivar_stk + (ivar_stk LE 0.0))
   var_bar = (w_sum GT 0.0)*total(weights^2*var_stk, 3)/(w_sum + $
                                                         (w_sum EQ 0.0))^2
   ivar_bar = (var_bar GT 0.0)/(var_bar + (var_bar LE 0.0))
ENDELSE
badpix = WHERE(finite(sky_bar) NE 1, nbad)
IF nbad NE 0 THEN sky_bar[badpix] = 0.0D

wstruct = long_wstruct(hdr_arc)
;;savefile = 'SINFONI_0.025_K.sav' ;; Don't write out solutions
IF NOT KEYWORD_SET(savefile) THEN savefile = 0

; generate left and right edge of slits
traceset2xy, tset_slits[0], rows, left_edge
traceset2xy, tset_slits[1], rows, right_edge
dim_slit = size(tset_slits[0].coeff, /dim)
IF n_elements(dim_slit) EQ 1 THEN nslit = 1 ELSE nslit = dim_slit[1]


;;KLUDGE_WAVE = 1
IF NOT KEYWORD_SET(KLUDGE_WAVE) THEN BEGIN
   xfit = long_waveimg(sky_bar, ivar_bar, tset_slits, wstruct, savefile $
                       , fwhmset = fwhmset, qafile = qafile $
                       , arc_inter = arc_inter, tweak_arc = tweak_arc $
                       , ARCTRC_POS = arctrc_pos, CALIB = CALIB)
   pixset = long_wavepix(sky_bar, tset_slits, FWHM = fwhmset.median $
                         , box_radius = wstruct.radius $
                         , sig_thresh = wstruct.sig_wpix $
                         , pkwdth = wstruct.pkwdth $
                         , TOLER = wstruct.TOLER $
                         , nsig = 7.0 $
                         , arc_ncoeff = arc_ncoeff, CHK = CHK) 
   piximg = long_wpix2image(pixset, tset_slits, XFIT = xfit $
                            , waveimg = waveimg)
ENDIF ELSE BEGIN
   pixset = long_wavepix(sky_bar, tset_slits, FWHM = wstruct.pkwdth/1.3d $
                         , box_radius = wstruct.radius $
                         , sig_thresh = wstruct.sig_wpix $
                         , pkwdth = wstruct.pkwdth $
                         , TOLER = wstruct.TOLER $
                         , nsig = 7.0d $
                         , arc_ncoeff = arc_ncoeff, CHK = CHK) 
   piximg = long_wpix2image(pixset, tset_slits)
   IF wstruct.BAND EQ 'H' THEN dlam = 0.079d/double(ny) $
   ELSE IF wstruct.BAND EQ 'K' THEN dlam = 0.122d/double(ny) $
   ELSE message, 'this order not yet supported'
   waveimg = (wstruct.wave_cen - double(ny)/2.0d*dlam) + piximg*dlam

ENDELSE

splog, 'Need to add heliocentric corrections here'

splog, 'Elapsed time = ', systime(1)-t0, ' sec'

IF KEYWORD_SET(OUTFILE) THEN BEGIN
;--------------                                
;;  Write output to wavefile
   splog, 'Writing output file'
   sxdelpar, hdr_arc, 'XTENSION'
   sxdelpar, hdr_arc, 'NAXIS2'
   sxdelpar, hdr_arc, 'NAXIS1'
   sxdelpar, hdr_arc, 'NAXIS'
   sxdelpar, hdr_arc, 'O_BZERO'
   sxaddpar, hdr_arc, 'BITPIX', -64

   mwrfits, waveimg, outfile, hdr_arc, /create 
   mwrfits, pixset, outfile 
   mwrfits, fwhmset, outfile 
ENDIF   

RETURN, waveimg
END

;order = esopar(hdr, 'HIERARCH ESO INS GRAT ORDER')
;wave_cen =  esopar(hdr, 'HIERARCH ESO INS GRAT WLEN')
;CASE setup OF
;   'LONG_SLIT_RED': BEGIN
;      restore, file = getenv('XIDL_DIR') + $
;               '/Spec/Longslit/calib/linelists/sofi_HK.sav'
;      title_string = 'SOFI H+K grism'
;      arc_offset = 0.0
;   END
;   'LONG_SLIT_H': BEGIN 
;      restore, file = getenv('XIDL_DIR') + $
;               '/Spec/Longslit/calib/linelists/sofi_H.sav'
;      title_string = 'SOFI H grism'
;      arc_offset = 30
;   END
;   'LONG_SLIT_K': BEGIN 
;      splog, 'K is not yet supported!! WAVELENGTHS ARE JUNK!!!'


;; Hack to create bogus wavelengths for now,  ;; H-band
;IF order EQ 3.0 THEN dlam = 0.079d4/double(ny) $
;ELSE IF order EQ 2.0 THEN dlam = 0.122d4/double(ny) $
;ELSE message, 'this order not yet supported'
;arc_offset = 0

;;waveimg = (wave_cen*1.0d4 - double(ny)/2.0d*dlam) + piximg*dlam
;;RETURN, waveimg
;; 
      ;; Later use this code for the correct procedure
      ;restore, file = getenv('XIDL_DIR') + $
      ;         '/Spec/Longslit/calib/linelists/sofi_K.sav'
      ;title_string = 'SOFI K grism'
                                ;arc_offset = 0
;   END
;   ELSE: message, 'ERROR: Unknown setup'
;ENDCASE


;if not keyword_set(LINLIST) then $
;   linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/oh_linelist.lst' 

;; Determine parameters for wavelengths. These are the same parameters
;; as used in the definition of the line-list in make_oh_linelist.pro
;slit_str = strcompress(esopar(hdr, 'HIERARCH ESO INS OPTI1 ID'), /rem)
;CASE slit_str OF 
;   'slit_1': slit_width = 1.0d
;   'slit_0.6_tilted': slit_width = 0.6d
;   ELSE: message, 'Unrecognized slit'
;ENDCASE
;plate_scale = 0.147d
;fnslit = slit_width/plate_scale
;pkwdth = 1.3d*fnslit
;toler = fnslit/2.0d > 2.0d
;THIN = 1
;FWEIGHT = 0
;NSIG = 5.0d
;FWHM = fnslit
;FORDR = 9

;order = esopar(hdr, 'HIERARCH ESO INS GRAT ORDER')
;wave_cen =  esopar(hdr, 'HIERARCH ESO INS GRAT WLEN')
;CASE setup OF
;   'LONG_SLIT_RED': BEGIN
;      restore, file = getenv('XIDL_DIR') + $
;               '/Spec/Longslit/calib/linelists/sofi_HK.sav'
;      title_string = 'SOFI H+K grism'
;      arc_offset = 0.0
;   END
;   'LONG_SLIT_H': BEGIN 
;      restore, file = getenv('XIDL_DIR') + $
;               '/Spec/Longslit/calib/linelists/sofi_H.sav'
;      title_string = 'SOFI H grism'
;      arc_offset = 30
;   END
;   'LONG_SLIT_K': BEGIN 
;      splog, 'K is not yet supported!! WAVELENGTHS ARE JUNK!!!'


;; Hack to create bogus wavelengths for now,  ;; H-band
;IF order EQ 3.0 THEN dlam = 0.079d4/double(ny) $
;ELSE IF order EQ 2.0 THEN dlam = 0.122d4/double(ny) $
;ELSE message, 'this order not yet supported'
;arc_offset = 0

;;waveimg = (wave_cen*1.0d4 - double(ny)/2.0d*dlam) + piximg*dlam
;;RETURN, waveimg
;; 
      ;; Later use this code for the correct procedure
      ;restore, file = getenv('XIDL_DIR') + $
      ;         '/Spec/Longslit/calib/linelists/sofi_K.sav'
      ;title_string = 'SOFI K grism'
                                ;arc_offset = 0
;   END
;   ELSE: message, 'ERROR: Unknown setup'
;ENDCASE


;if not keyword_set(LINLIST) then $
;   linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/oh_linelist.lst' 

;; Determine parameters for wavelengths. These are the same parameters
;; as used in the definition of the line-list in make_oh_linelist.pro
;slit_str = strcompress(esopar(hdr, 'HIERARCH ESO INS OPTI1 ID'), /rem)
;CASE slit_str OF 
;   'slit_1': slit_width = 1.0d
;   'slit_0.6_tilted': slit_width = 0.6d
;   ELSE: message, 'Unrecognized slit'
;ENDCASE
;plate_scale = 0.147d
;fnslit = slit_width/plate_scale
;pkwdth = 1.3d*fnslit
;toler = fnslit/2.0d > 2.0d
;THIN = 1
;FWEIGHT = 0
;NSIG = 5.0d
;FWHM = fnslit
;FORDR = 9

