;+
; NAME:
;   deimos_process_pixflat
;
; PURPOSE:
;   Process pixel-to-pixel response flat to deal with 'halos' around dust grains
;
; CALLING SEQUENCE:
;   deimos_process_pixflat
; 
; INPUTS:
;
; OUTPUTS:
; processed_pix_mult_flat.fits.gz
;
; MODIFICATION HISTORY:
;    code from DEIMOS_MASK_CALIBRATE - significant changes
;-
pro deimos_process_pixflat

	which,'deimos_pixflat

  infilename = 'deimos_pix_mult_flat.fits.gz'
  maskfilename= 'processed_pix_mult_flat.fits'
  mwrfits, junk, maskfilename, /create ;dummy primary

  calib_data = getenv('CALIB_DATA')+'/'

  for i=1,8 do begin

     pixmap = mrdfits(infilename, i, /silent)+1.

; just throw away deviant pixels in pixflat
    pixflathigh = (pixmap GT 1.2)
; tends to be a little halo around dust particles - do no pixflat there
    pixflathalo=dilatemask(pixflathigh,4) AND (pixmap LT 1)
    haloexists=total(pixflathalo) gt 0
    if haloexists then pixmap[where(pixflathalo)]=1.
     whzero=where(pixmap eq 0. OR pixmap eq -1.,zeroct)	
	if zeroct gt 0 then pixmap[whzero]=1.


;	mwrfits, floatcompress(pixmap-1., ndig=9), maskfilename
	mwrfits, pixmap-1., maskfilename

; USING COMPRESSED FILE - REPLACE WITH THIS TO REVERSE
            print, 'wrote chip:', i
  endfor  

; USING COMPRESSED FILE - COMMENT TO REVERSE
  spawn, 'gzip -f '+maskfilename

  return 
end







