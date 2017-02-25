
; NAME:
;   lris_superbias
;
; PURPOSE:
;   Cosntruct a mean bias from a list of images.
;
; CALLING SEQUENCE:
;   niri_superdark, filenames, outfile, [ sigrej=, maxiter=, /verbose ]
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
;   lris_proc
;   mwrfits
;   splog
;   sxaddpar
;
; REVISION HISTORY:
;   06-Jun-2006  Written by J. Hennawi (UCB)
;-
;------------------------------------------------------------------------------
FUNCTION nirspec_superdark, filenames, sigrej = sigrej, maxiter = maxiter $
                            , verbose = verbose
;----------
; Set defaults

nfile = n_elements(filenames)
if (nfile EQ 0 OR size(filenames, /tname) NE 'STRING') then $
  message, 'FILENAMES must be set'

if (NOT keyword_set(sigrej)) then begin
    if (nfile LE 2) then sigrej = 1.0 $ ; Irrelevant for only 1 or 2 files
    else if (nfile EQ 3) then sigrej = 1.1 $
    else if (nfile EQ 4) then sigrej = 1.3 $
    else if (nfile EQ 5) then sigrej = 1.6 $
    else if (nfile EQ 6) then sigrej = 1.9 $
    else sigrej = 2.0
endif
if (NOT keyword_set(maxiter)) then maxiter = 3

t0 = systime(1)
   
;----------
; Loop over all files

splog, 'Computing bias from ', nfile, ' files'
for ifile = 0L, nfile-1L do begin
    splog, 'Working on file ', filenames[ifile]
    
;     Read this file
;    hdr1 = headfits(filenames[ifile])
;    image1 = transpose(mrdfits(filenames[ifile], 1))
    niri_proc, filenames[ifile], image1, hdr = hdr1 $
               , verbose = verbose
    if (NOT keyword_set(image1)) then $
      message, 'Unable to read file '+filenames[ifile]
    dims1 = size(image1, /dimens)
    
;     If this is the first file, then construct other necessary images
    if (ifile EQ 0) then begin
        hdr = hdr1
        dims = dims1
        imgarr = make_array(dimension = [dims, nfile], /float)
;          inmask = make_array(dimension = [dims, nfile], /byte)
        sxaddpar, hdr, 'NEXP', nfile, 'Number of exposures in this file', $
                  before = 'EXPTIME'
    endif else begin
        if (total(dims NE dims1) NE 0) then $
          message, 'Inconsistent dimensions for input images'
    endelse
    
    imgarr[*, *, ifile] = image1
;      inmask[*, *, ifile] = invvar1 LE 0
    
    sxaddpar, hdr, string(ifile+1, format = '("FILE",i2.2)'), $
              fileandpath(filenames[ifile]), $
              ' File number ' + strtrim(string(ifile)), before = 'EXPTIME'
endfor                          ; End loop over files

IF nfile GT 1 THEN $
  imgfinal = djs_avsigclip(imgarr, 3, sigrej = sigrej, maxiter = maxiter) $
;                             inmask = inmask) $
ELSE imgfinal = imgarr0

; Write output file
;  mwrfits, imgfinal, outfile, hdr, /create

splog, 'Elapsed time = ', systime(1)-t0, ' sec'

return, imgfinal
end
;------------------------------------------------------------------------------
