; NAME:
;   long_wavesolve
;
; PURPOSE:
;  Full processing of an arc image.  This includes fitting the lines
;  and generating a wavelength image.
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;  wavefile -- Name of output file for wavelength image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;  long_proc
;  long_waveimg
;   
; REVISION HISTORY:
;   10-Mar-2005  Written by J. Hennawi (UCB) and Scott Burles (MIT)
;-
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
pro long_wavesolve, filename, wavefile, slitfile = slitfile $
                    , biasfile = biasfile, pixflatfile = pixflatfile $
                    , verbose = verbose, LINELIST=linelist, CHK = CHK $
                    , REID_FILE= reid_file, BIN_RATIO=bin_ratio, MED_ERR=med_err $
                    , ARCTRC_POS = arctrc_pos, CORR_PEAK=corr_peak $
                    , ARC_INTER=arc_inter, PROCARC=procarc, TWEAK_ARC=tweak_arc $
                    , calib = calib, NOSHIFT = NOSHIFT1


  IF keyword_set(noshift1) THEN NOSHIFT = NOSHIFT1

t0 = systime(1)

splog, 'Computing wavelength solution from file ', filename
splog, 'Using linelist ', linelist
;;----------
;; Read the arc image
if not keyword_set(PROCARC) then begin
    long_proc, filename, arcimg, arcivar, hdr = hdr $
               , biasfile = biasfile, verbose = verbose, pixflatfile = pixflatfile
endif else begin
    arcimg = xmrdfits(filename, 0 )
    arcivar = xmrdfits(filename, 1, hdr)
;    restore, filename
endelse
dims = size(arcimg, /dimen) 
nx = dims[0]
ny = dims[1]

;; Parse the header information to set paramters for the wavelength structure
wstruct = long_wstruct(hdr, LINELIST = linelist, REID_FILE = reid_file $
                       , BIN_RATIO = bin_ratio)
qafile = repstr(wavefile, '.fits', '.ps')
savefile = repstr(wavefile, '.fits', '.sav')
;----------
; Read in slit structure 
tset_slits = xmrdfits(slitfile, 1, silent = (keyword_set(verbose) EQ 0))

IF NOT KEYWORD_SET(NOSHIFT) THEN $
   xshift = long_xcorr_slits(arcimg, tset_slits, /shift) 
;----------
; Compute wavelength solution
xfit = long_waveimg(arcimg, arcivar, tset_slits, wstruct, savefile $
                    , fwhmset = fwhmset, qafile = qafile, arc_inter=arc_inter,$
                    tweak_arc=tweak_arc, ARCTRC_POS=arctrc_pos, calib = calib, $
                   CORR_PEAK=corr_peak)

; Now construct pixel wavelengths and apply wavelength solutions
;arc_ncoeff=2
pixset = long_wavepix(arcimg, tset_slits, fwhm = fwhmset.MEDIAN $
                      , arctrc_pos=arctrc_pos, med_err=med_err $
                      , box_radius = wstruct.radius $
                      , sig_thresh = wstruct.sig_wpix $
                      , arc_ncoeff = arc_ncoeff $
                      , pkwdth = wstruct.pkwdth $
                      , TOLER = wstruct.TOLER, CHK = CHK)
piximg = long_wpix2image(pixset, tset_slits, XFIT = xfit $
                         , waveimg = waveimg)

;--------------                                
;;  Write output to wavefile
splog, 'Writing output file'
sxdelpar, hdr, 'XTENSION'
sxdelpar, hdr, 'NAXIS2'
sxdelpar, hdr, 'NAXIS1'
sxdelpar, hdr, 'NAXIS'
sxdelpar, hdr, 'O_BZERO'
sxaddpar, hdr, 'BITPIX', -64

mwrfits, waveimg, wavefile, hdr, /create 
mwrfits, pixset, wavefile 
mwrfits, fwhmset, wavefile 

splog, 'Elapsed time = ', systime(1)-t0, ' sec'

return
end
;------------------------------------------------------------------------------
