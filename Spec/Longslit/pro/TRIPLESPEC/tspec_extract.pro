;+ 
; NAME:
; tspec_extract
;     Version 1.0
;
; PURPOSE:
;    Extracts the data for a TripleSpec frame
;
; CALLING SEQUENCE:
;   
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:

; OPTIONAL OUTPUTS:
;
; REVISION HISTORY:
;   Written by JFH as GNIRS_EXTRACT
;   26-Jan-2013 Modified by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION TSPEC_EXTRACT, img_minsky, sciivar, waveimg, thismask, sky_model $
                        , objstruct, plate_scale $
                        , MAXFWHM = MAXFWHM, MINFWHM = MINFWHM $
                        , TELLURIC = TELLURIC, PROFILE_STRUCT = PROFILE_STRUCT

IF KEYWORD_SET(TELLURIC) THEN ERRSCL = 0.01 ELSE ERRSCL = 1.0
;; Maximum FWHM is 1.3"
IF NOT KEYWORD_SET(MAXFWHM) THEN MAXFWHM = 1.3D/plate_scale
;; Minimum FWHM is 2 pixels
IF NOT KEYWORD_SET(MINFWHM) THEN MINFWHM = 2.0D
; Boxcar extraction radius
if NOT keyword_set(box_rad) then box_rad = 8L

nx = (size(sciivar))[1]
ny = (size(sciivar))[2]

xmap = findgen(nx) # replicate(1.0, ny)
nobj = n_elements(objstruct)
ntot = n_elements(objstruct)

; Define boundaries of this object for extraction
min_slit = min(xmap[WHERE(thismask)])
max_slit = max(xmap[WHERE(thismask)])
mincols = objstruct.xpos - objstruct.maskwidth - 1
maxcols = objstruct.xpos + objstruct.maskwidth + 1

scope = total(thismask EQ 1, 2)
mincol = floor(min(mincols)) > min(where(scope)) > min_slit
maxcol = ceil(max(maxcols)) < max(where(scope))  < max_slit

mask = (thismask EQ 1) AND (sciivar GT 0) AND (finite(img_minsky) EQ 1) $
  AND (abs(img_minsky) LT 5.0d4)
nc = maxcol - mincol + 1L
obj_profiles = fltarr(nc*ny)
    
; Boxcar fluxes used to normalize profile fit
;flux_center = extract_boxcar(img_minsky*(thismask EQ 1) $
;                             , objstruct.xpos, objstruct.ypos $
;                             , radius = box_rad)
;flux_left  = extract_boxcar(img_minsky*(thismask EQ 1) $
;                            , objstruct.xpos-2.*box_rad $
;                            , objstruct.ypos, radius = box_rad)
;flux_right = extract_boxcar(img_minsky*(thismask EQ 1) $
;                            , objstruct.xpos+2.*box_rad $
;                            , objstruct.ypos, radius = box_rad)

;flux = flux_center - 0.5*(flux_left + flux_right)

;pixtot = extract_boxcar(float(sciivar*0 + 1) $
;                        , objstruct.xpos, objstruct.ypos $
;                        , radius = box_rad) 
;mask_box = extract_boxcar(float(mask EQ 0) $
;                          , objstruct.xpos, objstruct.ypos $
;                          , radius = box_rad) NE pixtot
;fluxvar = extract_boxcar(1./(sciivar + (sciivar EQ 0)) $
;                         , objstruct.xpos, objstruct.ypos $
;                         , radius = box_rad)
;fluxivar = mask_box/(fluxvar + (fluxvar EQ 0))

flux = extract_boxcar(img_minsky*(mask EQ 1) $
                      , objstruct.xpos, objstruct.ypos, radius = box_rad)
varimg = 1.0/(sciivar*errscl + (sciivar EQ 0))
var_box = extract_boxcar(varimg*(mask EQ 1) $
                         , objstruct.xpos, objstruct.ypos, radius = box_rad)
pixtot = extract_boxcar(float(sciivar*0 + 1) $
                        , objstruct.xpos, objstruct.ypos, radius = box_rad) 
mask_box = extract_boxcar(float(mask EQ 0) $
                          , objstruct.xpos, objstruct.ypos $
                          , radius = box_rad) NE pixtot
box_denom = extract_boxcar((waveimg GT 0.0) $
                           , objstruct.xpos, objstruct.ypos, radius = box_rad)
wave = extract_boxcar(waveimg, objstruct.xpos, objstruct.ypos $
                      , radius = box_rad)/(box_denom + (box_denom EQ 0))
fluxivar = mask_box/(var_box + (var_box EQ 0))
; Gaussian profile fit for optimal extraction
;stop ;; Use something else than long_gprofile (i.e. tspec_gprofile)
obj_profiles = long_gprofile(img_minsky[mincol:maxcol, *] $
                             , (sciivar*mask*errscl)[mincol:maxcol, *] $
                             , waveimg[mincol:maxcol, *] $         
                             , objstruct.xpos - mincol $
                             , wave, flux, fluxivar $
                             , objstruct $
                             , hwidth = objstruct.maskwidth $
                             , fwhmfit = fwhmfit $
                             , thisfwhm =  $
                             ((objstruct.fwhm >  MINFWHM) <  MAXFWHM) $
                             , xnew = xnew, nccd = nccd $
                             , SN_GAUSS = SN_GAUSS $ ;1.0d10 $
                             , MED_SN2 = MED_SN2, /silent, wvmnx=[0.9,2.5])
; Update object trace and fwhm 
objstruct.xpos    = xnew + mincol
objstruct.fwhmfit = fwhmfit
objstruct.fwhm    = median(fwhmfit)

ipix = lindgen(nc) # replicate(1, ny) + $
  replicate(1, nc) # lindgen(ny)*nx + mincol

this_profile = reform(obj_profiles, nc, ny)

IF ARG_PRESENT(PROFILE_STRUCT) THEN BEGIN
   prof_proto = {objid: 0L, slitid: 0L, mincol: 0L, maxcol:0L $
                 , MED_SN2: 0.0, profile: fltarr(nx, ny)}
   IF nobj GT 1 THEN message, 'Problem with profile_struct'
   profile_struct = replicate(prof_proto, nobj)
   profile_struct.OBJID = objstruct.objid
   profile_struct.SLITID = objstruct.slitid
   profile_struct.MINCOL = MINCOL
   profile_struct.MAXCOL = MAXCOL
   profile_struct.MED_SN2 = MED_SN2 
   profile_struct.PROFILE[ipix] = obj_profiles
ENDIF

trace = objstruct.xpos ## replicate(1, nc)
objmask =  ((ipix mod nx GE (trace - 2.0*box_rad)) AND $
            (ipix mod nx LE (trace + 2.0*box_rad))) 

; ???? Should we use the model of the sky for the optimal extraction weights 
; rather than using the actual image? This procedure would not include noise
; from the object counts themselves. However for nearly all cases of interest
; the object will be much fainter than the sky. The reason that one might
; use the model for the optimal extraction is that it would mitigate
; the effects of bad pixels and CRs in the weights. The more correct thing 
; to do  would be to use the object profile to reject these bad pixels. For
; now I will just use the ivar image for the weights. 

;; gnirs_extract_optimal uses unit weights instead of ivar weights but 
;; does do optimal extraction using the object profile. This should make
;; things less noisy. 


; Optimal extraction
single_struct = gnirs_extract_optimal((thismask*waveimg)[mincol:maxcol, *] $
                                      , img_minsky[mincol:maxcol, *] $
                                      , (sciivar*thismask)[mincol:maxcol, *] $
                                      , this_profile $
                                      , mask[mincol:maxcol, *]*objmask  $
                                      , sky_model[mincol:maxcol, *] $
                                      , objstruct.xpos - mincol $
                                      ,  box_rad = box_rad)
;single_struct.trace = single_struct.trace + mincol
final_struct = struct_addtags(objstruct, single_struct)


RETURN, final_struct
END    
