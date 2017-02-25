; BUGS:
;
;------------------------------------------------------------------------------
; Main reduction routine.
;PRO gnirs_reduce1, filenames, scifile, acqfiles, flatfile = flatfile $
;                   , ZAP = ZAP, TELLURIC = TELLURIC, CHK = CHK $
;                   , VERBOSE = VERBOSE, TARGDIR = TARGDIR, WVCHK = WVCHK $
;                   , BOX_RAD = BOX_RAD1
;----------
; Set defaults
;;; 52-59 ;; Telluric
;;; 60-65 ;; Object 0801+1316

;; orderfile
reduxpath = '/Users/joe/REDUX_OTHER/P200_TSPEC/'
;;scipath = reduxpath + 'Science/'
scipath = reduxpath
IF file_test(scipath, /dir) EQ 0 THEN spawn, 'mkdir ' + scipath
orderfile = reduxpath + 'tspec-orders.fits'
;; define flat files
pixflatfile = reduxpath + 'tspec-pixflat.fits'
illumflatfile = reduxpath + 'tspec-illumflat.fits'
;; Define object files
path = '/Users/joe/DATA/P200_DATA/2011416/'
prefix = 'tspec110416_'
;; This is for the object
objind =  60 + lindgen(6)
obj_str = string(objind, FORMAT = '(I4.4)')
filenames = path + prefix + obj_str + '.fits'
indx = [0, 1]
filenames = filenames[indx]
;; -------------------
;; NOTE ON TELLURICS
;; -------------------
;; If you want to run the Telluric, you still need to pass the regular
;; object files to get the wavelength map. So just uncomment the lines
;; below and pass in both filenames and telluric
;; This is for Telluric
;tellind =  52 + lindgen(8)
;tel_str = string(tellind, FORMAT = '(I4.4)')
;telluric = path + prefix + tel_str + '.fits'
;indx = [4, 5]
;telluric = telluric[indx]
IF KEYWORD_SET(TELLURIC) THEN $
   scifile = scipath + 'tel-' + $
             prefix + $
             strcompress(string(objind[indx[0]], FORMAT = '(I4.4)'), /rem) $
             + '-' $
             + strcompress(string(objind[indx[1]], FORMAT = '(I4.4)'), /rem) $
             + '.fits' $
ELSE scifile = scipath + 'sci-' + $
               prefix + $
               strcompress(string(objind[indx[0]], FORMAT = '(I4.4)'), /rem) $
               + '-' $
               + strcompress(string(objind[indx[1]], FORMAT = '(I4.4)'), /rem) $
               + '.fits' 


IF KEYWORD_SET(BOX_RAD1) THEN BOX_RAD = BOX_RAD1 $
ELSE IF KEYWORD_SET(TELLURIC) THEN BOX_RAD = 15L $
ELSE BOX_RAD = 8L

t0 = systime(1)
;------------
; Read in order set structure and create ordermask
tset_slits = xmrdfits(orderfile, 1)
dimt = size(tset_slits.coeff, /dimen)
norders = dimt[1]
;------------
; Account for possible flexure between slits (determined from flat 
; field images) and the data. 
; ???? How necessary is this?? 
IF NOT KEYWORD_SET(NOSHIFT) AND NOT KEYWORD_SET(TELLURIC) THEN BEGIN
   tspec_proc, filenames[0], sciimg, hdr = schidr
   xshift = long_xcorr_slits(sciimg, tset_slits, /shift)
ENDIF
plate_scale = 0.234D            ; TSPEC plate scale
dimt = size(tset_slits.coeff, /dimen)
ordermask = tspec_ordermask(tset_slits, order_vec = order_vec)
;----------
; Perform AB sky-subtraction 
CHK = 1
wvchk = 1
AB = TSPEC_SKYSUB(filenames, orderfile, pixflatfile, illumflatfile $
                  , tset_slits, /OBJMASK, IVAR_AB = IVAR_AB $
                  , OBJ_POS = OBJ_POS, OBJ_NEG = OBJ_NEG $
                  , WAVEIMG = WAVEIMG, SKY_MODEL = SKY_MODEL $
                  , SCIHDR = SCIHDR $
                  , TELLURIC = TELLURIC, ZAP = ZAP, CHK = CHK $
                  , VERBOSE = VERBOSE, HDR = HDR, targdir = targdir $
                  , MINFWHM = MINFWHM, MAXFWHM = MAXFWHM $
                  , SKY_RESIDS = SKY_RESIDS, WVCHK = WVCHK)
;; Left off here. Still debugging tspec_echextobj
;; Need to make sure S/N stuff is correct and then use the 
;; gnirs_extract in place of more complicated localskysub
tspec_echextobj, AB, sky_resids, ivar_AB, scihdr, sky_model $
                 , piximg, waveimg, tset_slits, obj_pos $
                 , outfil = outfil, boxcar = boxcar $
                 , BOX_RAD = BOX_RAD1, STD = TELLURIC, CHK = CHK $
                 , AIRVAC_CORR = AIRVAC_CORR, HELIO_CORR = HELIO_CORR $
                 , SPCCHK = SPCCHK, DEBUG = DEBUG, VERBOSE = verbose

final_struct = 0
;----------
; Loop over each order and extract
FOR iorder = 0L, norders-1L DO BEGIN
    thismask = (ordermask EQ order_vec[iorder])
    ii_pos = where(obj_pos.SLITID EQ order_vec[iorder], npos)
    ;;ii_neg = where(obj_neg.SLITID EQ order_vec[iorder], nneg)
    if npos NE 1 then begin
       message, 'Error with number of objects npos=', npos, ' on order #' $
                , order_vec[iorder]
    endif
    ;if nneg NE 1 then begin
    ;    message, 'Error with number of objects nneg=', nneg, ' on order #' $
    ;             , order_vec[iorder]
    ;endif
    extract_pos = gnirs_extract(AB-sky_resids, ivar_AB, waveimg $
                                , thismask, sky_model, obj_pos[ii_pos] $
                                , plate_scale, TELLURIC = TELLURIC)
    ;;extract_neg = gnirs_extract(-AB+sky_resids, ivar_AB, waveimg $
    ;;                            , thismask, sky_model, obj_neg[ii_neg] $
    ;;                            , plate_scale, TELLURIC = TELLURIC)
    final_struct = struct_append(final_struct, extract_pos)
    ;;final_struct = struct_append(final_struct, extract_neg)
ENDFOR

;----------
; Write output file

splog, 'Writing FITS file ', scifile
mwrfits, float(AB), scifile, scihdr, /create, /silent
mwrfits, float(sky_resids), scifile, /silent
mwrfits, float(ivar_AB), scifile, /silent
mwrfits, float(waveimg), scifile, /silent
mwrfits, final_struct, scifile, /silent

;gnirs_plotsci, scifile, hard_ps = repstr(scifile, '.fits', '.ps') $
;               , box = keyword_set(TELLURIC)

splog, 'Elapsed time = ', systime(1)-t0, ' sec'

;;RETURN
END

