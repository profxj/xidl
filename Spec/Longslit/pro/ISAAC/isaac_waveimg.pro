FUNCTION ISAAC_WAVEIMG, arcimg, arcivar, tset_slits, hdr, scifile $
                        , piximg = piximg, MXSHIFT = MXSHIFT $
                        , SIGREJ = SIGREJ, WAVEFILE = WAVEFILE $
                        , BOX_RAD = BOX_RAD, CHK = CHK $
                        , CALIB = CALIB, pixset = pixset
  
t0 = systime(1)

if not keyword_set(MXSHFT) then MXSHFT = 15L
if not keyword_set(SIGREJ) then sigrej = 2.
if (NOT keyword_set(box_rad)) then box_rad = 5L
arc_ncoeff = 5 ;; ???more curvature in arc lines use higher order?????

badpix = WHERE(finite(arcimg) NE 1, nbad)
IF nbad NE 0 THEN arcimg[badpix] = 0.0D

dims = size(arcimg, /dimens)
nx = dims[0]
ny = dims[1]

wstruct = long_wstruct(hdr)
IF strmatch(scifile, '*tel-*') THEN wavefile = repstr(scifile, 'tel', 'wave') $
ELSE wavefile = repstr(scifile, 'sci', 'wave')

qafile = repstr(wavefile, '.fits', '.ps')
savefile = repstr(wavefile, '.fits', '.sav')

; generate left and right edge of slits
traceset2xy, tset_slits[0], rows, left_edge
traceset2xy, tset_slits[1], rows, right_edge
dim_slit = size(tset_slits[0].coeff, /dim)
IF n_elements(dim_slit) EQ 1 THEN nslit = 1 ELSE nslit = dim_slit[1]


;;KLUDGE_WAVE = 1

IF NOT KEYWORD_SET(KLUDGE_WAVE) THEN BEGIN
   xfit = long_waveimg(arcimg, arcivar, tset_slits, wstruct, savefile $
                       , fwhmset = fwhmset, qafile = qafile $
                       , arc_inter = arc_inter, tweak_arc = tweak_arc $
                       , ARCTRC_POS = arctrc_pos, CALIB = CALIB)
   pixset = long_wavepix(arcimg, tset_slits, FWHM = fwhmset.median $
                         , box_radius = wstruct.radius $
                         , sig_thresh = wstruct.sig_wpix $
                         , pkwdth = wstruct.pkwdth $
                         , TOLER = wstruct.TOLER $
                         , arc_ncoeff = arc_ncoeff, CHK = CHK) 
   piximg = long_wpix2image(pixset, tset_slits, XFIT = xfit $
                            , waveimg = waveimg)
ENDIF ELSE BEGIN
   pixset = long_wavepix(arcimg, tset_slits, FWHM = wstruct.pkwdth/1.3d $
                         , box_radius = wstruct.radius $
                         , sig_thresh = wstruct.sig_wpix $
                         , pkwdth = wstruct.pkwdth $
                         , TOLER = wstruct.TOLER $
                         , arc_ncoeff = arc_ncoeff, CHK = CHK) 
   piximg = long_wpix2image(pixset, tset_slits)
   IF wstruct.BAND EQ 'H' THEN dlam = 0.079d/double(ny) $
   ELSE IF wstruct.BAND EQ 'K' THEN dlam = 0.122d/double(ny) $
   ELSE message, 'this order not yet supported'
   waveimg = (wstruct.wave_cen - double(ny)/2.0d*dlam) + piximg*dlam

ENDELSE

splog, 'Need to add heliocentric corrections!!!"
splog, 'Elapsed time = ', systime(1)-t0, ' sec'


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

