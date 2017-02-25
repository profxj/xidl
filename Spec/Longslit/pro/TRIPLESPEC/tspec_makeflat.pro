;+
; NAME:
;        TSPEC_MAKEFLAT
;
;
; PURPOSE:
;        Generates the flat files for Palomar TSPEC spectra.
;
;
; CALLING SEQUENCE:
;
; INPUTS:
;        plan_file
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
; RESTRICTIONS
;
; SIDE EFFECTS:
;
;
; MODIFICATION HISTORY:
;  Written by JXP (April 2011)
;       Version 1.0
;  Modified December 2012
;-
pro tspec_makeflat, plan_file, orderfile, pixflatfile, illumflatfile, pixfile, chk = chk, $
                    flatfiles=flatfiles, illumfiles=illumfiles, scifile=scifile 


;data = data, flatfiles = flatfiles, orders = orderfile, $
;                    arcfile = arcfile, scifile = scifile, illum = illumfiles, $
;                    outfile = outfile, pixflat = pixflat, pixfile = pixfile, $
;                    CLOBBER = CLOBBER, CHK = CHK, VERBOSE = VERBOSE, $
;                    CHECK_INPUTS = check_inputs

  prog_name = 'tspec_makeflat.pro'

  ;; Set defaults
  ;if NOT keyword_set( outfile ) then begin
  ;  outfile = FIRE_GET_FILE_NAMES(illumfiles, /ILLUMFLAT) 
  ;endif
  ;if NOT keyword_set( pixflat ) then begin
  ;  pixflat = FIRE_GET_FILE_NAMES(flatfiles, /PIXFLATFILE) 
  ;endif
  ;if NOT keyword_set( pixfile ) then begin
;		pixfile = FIRE_GET_FILE_NAMES(flatfiles, /PIXIMGFILE) 
;  endif

  ;if keyword_set(VERBOSE) OR keyword_set(CHECK_INPUTS) then begin
;		print, ''
;		print, "++++++++++++++++++++++++++++++++++++"
;	  print, prog_name, ": /VERBOSE flag passed.  Set values:"
;		print, "Flats file(s): ", flatfiles
;		print, "Arcs file(s) (science file used instead!): ", arcfile
;		print, "Science file(s): ", scifile
;		print, "Illumination file(s): ", illumfiles
; 		print, "Order file(s): ", orderfile
; 		print, "Out file: ", outfile
; 		print, "Pixel out file: ", pixflat
; 		print, "Pixel image out file: ", pixfile
;		print, "++++++++++++++++++++++++++++++++++++"
;		print, ''
;		if keyword_set(CHECK_INPUTS) then RETURN
;  endif

  ;; Determine if we are in the scenario where nothing has to be run.  If so, let's just exit immediately
;  if NOT keyword_set(CLOBBER) AND FILE_TEST(outfile) AND FILE_TEST(pixflat) AND FILE_TEST(pixfile) then begin
;    print, prog_name, ": Files ", outfile, ", ", pixflat, ", and ", pixfile, " all already exists and keyword /clobber not passed!  Exiting without further ado..."
;    RETURN
;  endif

  ;; Output files
  if not keyword_set(PIXFLATFILE) then pixflatfile = 'tspec-pixflat.fits'
  if not keyword_set(ILLUMFLATFILE) then illumflatfile = 'tspec-illumflat.fits'
  if not keyword_set(PIXFILE) then pixfile = 'tspec-piximg.fits'

  ;; Read the plan file
  planstr = yanny_readone(plan_file, hdr=planhdr, /anony) 

  ;; Parse on Dome Flats
  indir = yanny_par(planhdr, 'indir')
  dflat = where(planstr.flavor EQ 'dflat', nflat)
  flatfiles = indir+planstr[dflat].filename
  tflat = where(planstr.flavor EQ 'twiflat', ntwi)
  if ntwi EQ 0 then begin
     print, 'tspec_makeflat: No twilight flats'
     if x_chkfil(illumflatfile) NE 1 then begin
        print, 'tspec_makeflat: Using existing '+illumflatfile
        flg_illum = 1
        tflat = 0L
     endif else begin
        print, 'tspec_makeflat:  No twilight files.'
        print, 'tspec_makeflat:  Using dome flats...'
        tflat = dflat
        ntwi = nflat
        flg_illum = 0
     endelse
  endif
  illumfiles = indir+planstr[tflat].filename

  sci = where(planstr.flavor EQ 'object', nsci)
  scifile = indir+planstr[sci[0]].filename

  norders = 5L
  
  if n_elements(orderfile) GT 1 then begin
     stop
     ;fire_siren, prog_name + ": Only one orderfile allowed!  Only using the first order file!"
  endif
  orderfile = orderfile[0] ;; This also converts string arrays of size 1 to string scalars

  tset_slits = xmrdfits(orderfile, 1)
  if size(tset_slits, /type) EQ 2 then begin
                                ;fire_siren, prog_name + ": ERROR reading in tset slits!  (Your order file is " + $
     ;            orderfile + ".  If this does not begin with OStr, you may have a problem...)." + $
                                ;            "  Continuing on, but I'm betting against your chances..."
  endif
  
  
  tspec_proc, scifile, sciimg
  
  ;; This finds the tilt of the arc/sky lines.  Works better with
  ;;   the sky since there are more lines.
  
  print, " "
  print, "tspec_makeflat: Fitting slit tilts..."
  print, " "
  
  piximg = tspec_makescipix(sciimg, tset_slits, pixset = pixset, chk = chk)
  ;;, arc2 = arcimg)
  mwrfits, piximg, pixfile, /create
  mwrfits, pixset, pixfile
    
;  endelse





;; Junk   
;  chk=1
  ;; Run with an illumflat and tweak slit boundaries
;  input = illumfiles
  ;junkflat = fire_superflat(input[0], orderfile, pixfile $
  ;                          , /skyflat, /tweak, chk=chk)
;  junkflat = fire_superflat(input, orderfile, pixfile $
;                            , /skyflat, illum = illumflat, /chk)



  ;;TOL_EDG=[0]
  ;;ximg = long_slits2x(tset_slits, edgmask = edgmask, TOL_EDG = TOL_EDG $
  ;;                    , nslit = norders)
    ordermask = tspec_ordermask(tset_slits) ;;, /fsr)


  ;; Remove EDG_TRIM pixels from each side. 
  EDG_TRIM=[1L,1L]
  IF KEYWORD_SET(EDG_TRIM) THEN BEGIN
     tset_slits=mrdfits(orderfile,1)
     FOR iorder = 0L, norders-1L DO BEGIN
        tset_slits[0].COEFF[0, iorder] = $
           tset_slits[0].COEFF[0, iorder] + EDG_TRIM[0]
        tset_slits[1].COEFF[0, iorder] = $
           tset_slits[1].COEFF[0, iorder] - EDG_TRIM[1]
     ENDFOR
     ordermask = tspec_ordermask(tset_slits) ;;, /fsr) 
;     mwrfits, ordermask, orderfile, /create
;     mwrfits, tset_slits, orderfile
  ENDIF
  tset_slits=mrdfits(orderfile,1)

  if not keyword_set(flg_illum) then begin
     print, " "
     print, "tspec_makeflat: Generating illumination function"
     print, " "
     ;;input = illumfiles
     ;;if (n_elements(input) GT 1) then begin
     illflat    = tspec_superflat(illumfiles, orderfile, pixfile  $
                                  , illum = illumflat, /skyflat $
                                  , chk = 1, tweak = 1)
;endif else begin
;     illflat    = fire_superflat(input, orderfile, pixfile  $
;                                 , illum = illumflat, /skyflat, chk = 1, tweak = 1)
;  endelse
     
     ;; Mask bad pixels
     ;;unit2=(edgmask OR illumflat GT 3.0 OR illumflat LE 0.2 OR ordermask EQ 0)
     unit2 = (illumflat GT 3.0 OR illumflat LE 0.2 OR ordermask EQ 0)
     unitind2 = WHERE(unit2, nunit2)
     IF nunit2 GT 0 THEN illumflat[unitind2] = 1.0d
     
     ;; Plot result
     print, " "
     print, "Displaying illumination flat"
     print, " "
     if keyword_set(CHK) then xatv, illumflat, /block, min = 0.5, max = 1.5
     
     ;; Save result
     mwrfits, illumflat, illumflatfile, /create
  endif

  print, " "
  print, "tspec_makeflat: Generating pixel flat"
  print, " "
  input   =  flatfiles
  ;;if (n_elements(input) GT 1) then begin

  ;stop
  flat    = tspec_superflat(flatfiles, orderfile, pixfile, chk = 1)
;    endif else begin
;       flat    = fire_superflat(input, orderfile, pixfile  $
;                              , chk=1)
;    endelse  

    ;; Mask bad pixels
    ;;unit1=(edgmask OR flat GT 3.0 OR flat LE 0.0 OR ordermask EQ 0)
  unit1 = (flat GT 3.0 OR flat LE 0.0 OR ordermask EQ 0)
  unitind1 = WHERE(unit1, nunit1)
  IF nunit1 GT 0 THEN flat[unitind1] = 1.0d
  
  ;; Plot result
  if keyword_set(CHK) then xatv, flat, min = 0.9, max = 1.1, /block
  
  ;; Save result
  mwrfits, flat, pixflatfile, /create




end
