pro deimos_hotpix, files, sigma = sigma, path = path
;+
; NAME:  deimos_hotpix
;
; PURPOSE:
;   make a hot pixel mask using DEIMOS dark frames 
;
; CATEGORY:
;   spec2d reduction
;
; CALLING SEQUENCE:
;   deimos_hotpix, files, [sigma=sigma, path=path]
; 
; INPUTS:
;   files - array of string names of unprocessed dark frames to be read
;           ex: ['d0610_0002.fits','d0610_0003.fits','d0610_0004.fits']
;
; KEYWORD PARAMETERS:
;   sigma - # of sigma used to flag a pixel - default=7
;   path - directory path to write the output fits files to - default=CALIB_DATA
;
; OUTPUTS:
;   "deimos_hotpix.fits" - hot pixel mask: 1 where good, 0 where bad
;   "dark.fits" - final processed dark frame (bias subtracted, trimmed, CR rejected, combined)
;
; MODIFICATION HISTORY:
;   alc 12jun02
;-
; files=['d0610_0002.fits','d0610_0003.fits','d0610_0004.fits','d0610_0005.fits']

; define output hotpix file
if NOT keyword_set(path) then path = getenv('CALIB_DATA')+'/'
hotpixfile = path+'deimos_hotpix.fits'

darkimagefile = path+'dark.fits'

mwrfits, junk, hotpixfile, /create ;dummy primary
mwrfits, junk, darkimagefile, /create ;dummy primary

nexp = n_elements(files)

for i = 1, 8 do begin

   ; read in all the exposures, bias subtract, trim and concatenate 
   chip=fltarr(2048,4096,nexp)   
   if NOT keyword_set(path) then path = './'
   for j = 0, nexp-1 do chip(*, *, j) = deimos_read_chip(path+files[j], i)

   ; CR rejection
   invar=chip*0+1./median(reform(chip,(size(chip))[5]))
   avg=avsigclip(chip,invar)

   ; make bad pixel mask
   if NOT keyword_set(sigma) then sigma = 7
   badmask = bytarr(2048, 4096) +1

   good = where(avg.flux le 1000) ; get rid of really high pixels
   clean = (avg.flux)[good] ; cleaner version of flux values in dark
   rand = randomu(seed, 1e5)*n_elements(clean) ; randomly select for statistics
   djs_iterstat, clean[rand], sigrej = 5, sigma = std
   print, 'stdev of this chip is ', std
   medflux = median(clean[rand])
   hot = where(avg.flux ge (medflux+sigma*std) OR avg.flux le (medflux-sigma/2.*std))
   badmask[hot] = 0

   mwrfits, badmask, hotpixfile
   mwrfits, avg.flux, darkimagefile
endfor



return
end



