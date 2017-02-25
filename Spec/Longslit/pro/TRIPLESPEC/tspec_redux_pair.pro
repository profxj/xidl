;; Reduces a pair of TripleSpec science frames using AB subtraction
;;   Default:  Only extracts the first set of data
;;   Can also work on a pair of Telluric frames, but science frames
;;   need to be supplied too
;PRO tspec_redux_pair, filenames, scifile, acqfiles, flatfile = flatfile $
;                   , ZAP = ZAP, TELLURIC = TELLURIC, CHK = CHK $
;                   , VERBOSE = VERBOSE, TARGDIR = TARGDIR, WVCHK = WVCHK $
;                   , BOX_RAD = BOX_RAD1
;  /EX_BOTH -- Extract both objects (no longer the default)
;----------
; Set defaults

pro tspec_redux_pair, in_filenames, telluric=in_telluric, $
                      EX_BOTH=ex_both, CLOBBER=clobber

reduxpath = './'
scipath = reduxpath + 'Science/'
;; orderfile
IF file_test(scipath, /dir) EQ 0 THEN spawn, 'mkdir ' + scipath
orderfile = reduxpath + 'tspec-orders.fits'

;; define flat files
pixflatfile = reduxpath + 'tspec-pixflat.fits'
illumflatfile = reduxpath + 'tspec-illumflat.fits'

;; Define object files
path = 'Raw/'

;; Filenames
filenames = path+in_filenames
if keyword_set(IN_TELLURIC) then telluric = path+in_telluric

;; -------------------
;; NOTE ON TELLURICS
;; -------------------
;; If you want to run the Telluric, you still need to pass the regular
;; object files to get the wavelength map. So just uncomment the lines
;; below and pass in both filenames and telluric

IF KEYWORD_SET(TELLURIC) THEN Begin
   ipos = strpos(in_telluric[0], '.fits')
   scifile = 'Science/tel-'+strmid(in_telluric[0],0,ipos+5) 
Endif ELSE begin
   ipos = strpos(in_filenames[0], '.fits')
   scifile = 'Science/sci-'+strmid(in_filenames[0],0,ipos+5) ;-tspec1072.fits'
endelse

a = file_search(scifile+'*', count=nfil)
if nfil GT 0 and not keyword_set(CLOBBER) then begin
  print, 'tspec_redux_pair: File exists -- ', scifile
  print, 'tspec_redux_pair: Use /CLOBBER to overwrite'
  print, 'tspec_redux_pair: Returning..'
  return
endif

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
   tspec_proc, filenames[0], sciimg, hdr = scihdr
   xshift = long_xcorr_slits(sciimg, tset_slits, /shift)
ENDIF else begin
   if keyword_set(TELLURIC) then scihdr = xheadfits(telluric[0])
endelse
plate_scale = 0.234D            ; TSPEC plate scale
dimt = size(tset_slits.coeff, /dimen)
ordermask = tspec_ordermask(tset_slits, order_vec = order_vec)
;----------
; Perform AB sky-subtraction 
CHK = 0
WVCHK = 0
;;wvchk = 1
AB = TSPEC_SKYSUB(filenames, orderfile, pixflatfile, illumflatfile $
                  , tset_slits, /OBJMASK, IVAR_AB = IVAR_AB $
                  , OBJ_POS = OBJ_POS, OBJ_NEG = OBJ_NEG $
                  , WAVEIMG = WAVEIMG, SKY_MODEL = SKY_MODEL $
                  , TELLURIC = TELLURIC, ZAP = ZAP, CHK = CHK $
                  , VERBOSE = VERBOSE, HDR = HDR, targdir = targdir $
                  , MINFWHM = MINFWHM, MAXFWHM = MAXFWHM $
                  , SKY_RESIDS = SKY_RESIDS, WVCHK = WVCHK $
                  , EX_BOTH=ex_both)
if size(AB,/type) EQ 2 then begin
   print, 'tspec_redux_pair:  No object.  Punting'
   return
endif

final_struct = 0

;----------
; Loop over each order and extract
FOR iorder = 0L, norders-1L DO BEGIN
    thismask = (ordermask EQ order_vec[iorder])
    ii_pos = where(obj_pos.SLITID EQ order_vec[iorder], npos)
    ii_neg = where(obj_neg.SLITID EQ order_vec[iorder], nneg)
    if npos NE 1 then begin
        message, 'Error with number of objects npos=', npos, ' on order #' $
                 , order_vec[iorder]
    endif
    if nneg NE 1 then begin
        message, 'Error with number of objects nneg=', nneg, ' on order #' $
                 , order_vec[iorder]
     endif
    extract_pos = tspec_extract(AB-sky_resids, ivar_AB, waveimg $
                                , thismask, sky_model, obj_pos[ii_pos] $
                                , plate_scale, TELLURIC = TELLURIC)
    final_struct = struct_append(final_struct, extract_pos)
    ;; Extract both frames?
    if keyword_set(EX_BOTH) then begin
       extract_neg = tspec_extract(-AB+sky_resids, ivar_AB, waveimg $
                                   , thismask, sky_model, obj_neg[ii_neg] $
                                   , plate_scale, TELLURIC = TELLURIC)
       final_struct = struct_append(final_struct, extract_neg)
    endif

ENDFOR

final_struct.arc_fil = ' '  ;; Padding to save properly

;----------
; Write output file
splog, 'Writing FITS file ', scifile
mwrfits, float(AB), scifile, scihdr, /create, /silent
mwrfits, float(sky_resids), scifile, /silent
mwrfits, float(ivar_AB), scifile, /silent
mwrfits, float(waveimg), scifile, /silent
mwrfits, final_struct, scifile, /silent

splog, 'Elapsed time = ', systime(1)-t0, ' sec'

;;RETURN
END

