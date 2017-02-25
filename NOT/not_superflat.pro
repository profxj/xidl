
; NAME:
;   long_superbias
;
; PURPOSE:
;   Cosntruct a mean bias from a list of images.
;
; CALLING SEQUENCE:
;   long_superbias, filenames, outfile, [ sigrej=, maxiter=, /verbose ]
;                 
;
; INPUTS:
;   filenames  - Input bias file names
;   outfile    - Output file name
;
; OPTIONAL INPUTS:
;   sigrej     - Rejection threshold in call to DJS_AVSIGCLIP() 
;   maxiter    - Number of rejection iterations; default to 3.
;   verbose    - Verbose if set
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   A sigma-clipped mean is computed for each pixel.
;
;   The output header is copied from the first input file, with an
;   additional NEXP keyword (for the number of input files) and FILE*
;   keywords that list each input file name.
;  
; EXAMPLES:
;
; BUGS:
;    
; PROCEDURES CALLED:
;   djs_avsigclip()
;   fileandpath()
;   long_proc
;   mwrfits
;   splog
;   sxaddpar
;
; REVISION HISTORY:
;   11-Mar-2005  Written by J. Hennawi (UCB), D. Schlegel (LBL)
;-
;------------------------------------------------------------------------------
pro not_superflat, filenames, outfile, verbose = verbose $
                   , superbiasfile = superbiasfile
  
   ;----------
   ; Set defaults

   nfile = n_elements(filenames)
   if (nfile EQ 0 OR size(filenames,/tname) NE 'STRING') then $
    message, 'FILENAMES must be set'

   t0 = systime(1)

   ;----------
   ; Loop over all files

   splog, 'Computing bias from ', nfile, ' files'
   for ifile=0L, nfile-1L do begin
      splog, 'Working on file ', filenames[ifile]
      not_proc, filenames[ifile], image1, ivar, superbiasfile = superbiasfile $
                , mask = mask, ihdr = hdr1
      if (NOT keyword_set(image1)) then $
         message, 'Unable to read file '+filenames[ifile]
      dims1 = size(image1, /dimens)
      flat1 = fltarr(dims1[0], dims1[1])
      FOR chip = 1L, 4L DO BEGIN
         ipix = WHERE(mask EQ chip, npix)
         flat1[ipix] = image1[ipix]/median(image1[ipix])
      ENDFOR
      image1 = 0
      ;; If this is the first file, then construct other necessary images
      if (ifile EQ 0) then begin
         hdr = hdr1
         dims = dims1
         flatarr = make_array(dimension = [dims, nfile], /float)
         sxaddpar, hdr, 'NEXP', nfile, 'Number of exposures in this file'
      endif else begin
         if (total(dims NE dims1) NE 0) then $
          message, 'Inconsistent dimensions for input images'
      endelse
      flatarr[*, *, ifile] = flat1
      sxaddpar, hdr, string(ifile+1, format = '("FILE",i2.2)'), $
                fileandpath(filenames[ifile]) $
                , ' File number ' + strtrim(string(ifile)), before = 'EXPTIME'
   endfor                       ;;End loop over files
   
   IF nfile GT 1 THEN flatfinal = djs_median(flatarr, 3) $
   ELSE flatfinal = flatarr

   flatarr = 0
   ;; Write output file
   mwrfits, flatfinal, outfile, hdr, /create
   
   splog, 'Elapsed time = ', systime(1)-t0, ' sec'

   return
end
;------------------------------------------------------------------------------
