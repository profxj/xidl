; NAME:
;   long_superflat
;
; PURPOSE:
;   Construct a pixel flat for multislit data. 
;
; CALLING SEQUENCE:
;   long_superflat, filenames, superpixflatfile, superillumflatfile, $
;    [ slitfile=, waveile=, biasfile=, $
;    npoly=, use_illum=, use_pixel=, tempdir=, $
;    sigrej=, maxiter=, /verbose, /write_flats ]
;
; INPUTS:
;   filenames    - Input flat-field file names, either dome or twilight flats
;   superpixflatfile  - Output file name for super pixel flat
;   superillumflatfile- Output file name for super illumination flat
;   wavefile     - File with wavelength solution
;
; OPTIONAL INPUTS:
;   slitfile     - File with trace sets describing the slit positions;
;                  if not set, then assume the entire image is a single slit.
;   biasfile     - Bias file to apply to raw images
;   pixflatfile  - Pixel flat to be used if an illumination 
;                  flat is being constructed
;     
;   tempdir      - Directory for temporary files when generating the flats;
;                  default to the current directory
;   npoly        - Number of polynomial terms for B-spline in the spatial
;                  direction; default to 3 for a cubic fit
;                  (But demand at least 10 pixels per row, on average,
;                  per degree of the polynomial.)
;   sigrej       - Rejection threshold in call to DJS_AVSIGCLIP()
;   maxiter      - Number of rejection iterations; default to 3.
;   use_illum    - Array the size of FILENAMES with 1s or 0s 
;                  which specifies which files are to be used for the 
;                  illumination function; default is to use all files
;   use_pixel    - Array the size of FILENAMES with 1s or 0s 
;                  which specifies which files are to be used for the 
;                  illumination function; default is to use all files
;   verbose      - Verbose if set
;   write_flats  - ???
;
; OUTPUTS: 
;
; OPTIONAL OUTPUTS: 
;
; COMMENTS:
;   Each slit of each input image has a model computed.  We divide each
;   input image by its corresponding model, and then take the sigma-clipped
;   mean of all those results.  The model is a B-spline in the wavelength
;   direction, and a low-order polynomial in the X direction.
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
;   bspline_iterfit()
;   bspline_valu()
;   djs_avsiglip()
;   fileandpath()
;   long_proc
;   long_slits2mask()
;   long_slits2x()
;   mrdfits
;   mwrfits
;   splog
;
; REVISION HISTORY:
;   10-Mar-2005 Written by J. Hennawi (UCB), D. Schlegel (LBL), S. Burles (MIT)
;   27-Feb-2012    Eliminated a 'bug' in the slitless flats with a
;   PIXIMG kludge [JXP]
;-
;------------------------------------------------------------------------------
pro long_superflat, filenames, superpixflatfile, superillumflatfile $
                    , slitfile = slitfile $
                    , waveimg = waveimg, piximg = piximg $
                    , wavefile = wavefile, biasfile = biasfile $
                    , darkfiles = darkfiles1 $
                    , npoly = npoly1 $
                    , use_illum = use_illum1, use_pixel = use_pixel1 $
                    ,tempdir = tempdir $
                    , sigrej = sigrej, maxiter = maxiter $
                    , verbose = verbose, write_flats = write_flats $
                    , CHK = CHK, SLITSAMP = slitsamp $
                    , TOL_EDG = TOL_EDG2 $
                    , NIRSPEC = NIRSPEC, SOFI = SOFI, ISAAC = ISAAC $
                    , LUCI = LUCI, SINFONI = SINFONI $
                    , TSPEC=tspec 


  if total(use_illum1) + total(use_pixel1) EQ 0 then return

   nfile = n_elements(filenames)
   if (nfile EQ 0 OR size(filenames,/tname) NE 'STRING') then $
      message, 'FILENAMES must be set'

   ;; For IR datasets allow darkfiles to be an array of files
   IF KEYWORD_SET(darkfiles1) THEN BEGIN
      CASE n_elements(darkfiles1) OF
         1: darkfiles = replicate(darkfiles1, nfile)
         nfile: darkfiles = darkfiles1
         ELSE: message, 'Darkfile must have either 1 element or nfile elements'
      ENDCASE
   ENDIF
   if (n_elements(use_illum1) NE 0) then use_illum = use_illum1 $
    else use_illum = lonarr(nfile)+1
   if (n_elements(use_pixel1) NE 0) then use_pixel = use_pixel1 $
    else use_pixel = lonarr(nfile)+1
   IF n_elements(TOL_EDG2) GT 0 THEN TOL_EDG1 = TOL_EDG2 $
   ELSE TOL_EDG1 = 2.5d ;; Default to masking 2.5 pixels
   n_illum = long(total(use_illum))
   splog, 'Computing illumination flat from ', n_illum, ' files' 
   
   n_pixel = long(total(use_pixel))
   splog, 'Computing pixel flat from ', n_pixel, ' files' 
 
   if (keyword_set(tempdir)) then $
    spawn, '\mkdir -p '+tempdir

   ;----------
   ; Set defaults

   if (NOT keyword_set(npoly1)) then npoly1 = 7
   if (NOT keyword_set(sigrej)) then begin
       if (n_pixel LE 2) then sigrej_pixel = 1.0 $ 
       ; Irrelevant for only 1 or 2 files
       else if (n_pixel EQ 3) then sigrej_pixel = 1.1 $
       else if (n_pixel EQ 4) then sigrej_pixel = 1.3 $
       else if (n_pixel EQ 5) then sigrej_pixel = 1.6 $
       else if (n_pixel EQ 6) then sigrej_pixel = 1.9 $
       else sigrej_pixel = 2.0
       
       if (n_illum LE 2) then sigrej_illum = 1.0 $ 
       ; Irrelevant for only 1 or 2 files
       else if (n_illum EQ 3) then sigrej_illum = 1.1 $
       else if (n_illum EQ 4) then sigrej_illum = 1.3 $
       else if (n_illum EQ 5) then sigrej_illum = 1.6 $
       else if (n_illum EQ 6) then sigrej_illum = 1.9 $
       else sigrej_illum = 2.0
   endif else begin
       sigrej_pixel = sigrej
       sigrej_illum = sigrej
   endelse

   t0 = systime(1)

   ;--------------------
   ; Loop over all files
   ;--------------------

   file_counter = 0L
   for ifile = 0L, nfile-1L do begin
       if ifile EQ 1 then tempillum0 = temp_illum
       fname = fileandpath(filenames[ifile])
       splog, 'Working on file ', fname
       temp_illum = djs_filepath('tempillum-'+fname, root_dir=tempdir)
       temp_spec  = djs_filepath('tempspec-'+fname, root_dir=tempdir)
       temp_pixel = djs_filepath('temppixel-'+fname, root_dir=tempdir)

       splog, 'Storing illumination function in ', temp_illum
       splog, 'Storing spectral dependence in ', temp_spec

;     Read this file
       IF KEYWORD_SET(NIRSPEC) THEN $
          nirspec_proc, filenames[ifile] $
                        , image1, invvar1, hdr = hdr1, darkfile = darkfiles[ifile] $
                        , verbose = verbose $
       ELSE IF KEYWORD_SET(SOFI) THEN $
          sofi_proc, filenames[ifile], image1, invvar1, hdr = hdr1 $
                     , darkfile = darkfiles[ifile], verbose = verbose $
       ELSE IF KEYWORD_SET(LUCI) OR KEYWORD_SET(ISAAC) $
       THEN BEGIN
          image1  = mrdfits(filenames[ifile], 0, hdr1)
          invvar1 = mrdfits(filenames[ifile], 1)
       ;;ENDIF ELSE IF KEYWORD_SET(ISAAC) THEN $
       ;;   isaac_proc, filenames[ifile], image1, invvar1, hdr = hdr1 $
       ;;               , verbose = verbose $
       ENDIF ELSE IF KEYWORD_SET(TSPEC) THEN $
             tspec_proc, filenames[ifile], image1, invvar1, hdr = hdr1 $
                         , verbose = verbose $
       ELSE IF KEYWORD_SET(SINFONI) THEN $
          ;; May want to add dark subtraction here, but none are taken
          sinfoni_proc, filenames[ifile], image1, invvar1, hdr = hdr1 $
          , verbose = verbose, darkfile = darkfiles[ii] $
       ELSE long_proc, filenames[ifile] $
                       , image1, invvar1, hdr = hdr1, biasfile = biasfile $
                       , verbose = verbose 
       telescope = strcompress(sxpar(hdr1, 'TELESCOP'), /rem)
       ;; Check for Kast (needed for long_flatfield_specillum :: JXP
       ;; 21 Aug 2013)
       ;; Also need to turn this on for MMT/Blue Channel ;; KHRR 23 Apr 2014
       ;; ISIS/WHT likes this too MF 2015
       if ifile EQ 0 then begin
          if (stregex(sxpar(hdr1,'INSTRUME'),'.*kast.*',$
                         /boolean,/fold_case) eq 1) or $
             (stregex(sxpar(hdr1,'VERSION'),'kast*', /bool, /fold_case) EQ 1) or $
             (stregex(sxpar(hdr1,'INSTRUME'),'.*mmtbluechan.*',$
                      /boolean,/fold_case) eq 1) or $
             (stregex(sxpar(hdr1,'INSTRUME'),'.*ISIS*.',/boolean,/fold_case) eq 1) then $
                kast = 1 else kast = 0
       endif
       ;;
       if (NOT keyword_set(image1)) then $
          message, 'Unable to read file ' + filenames[ifile]
       dims1 = size(image1, /dimens)
;     If this is the first file, then construct other necessary images
       if (ifile EQ 0) then begin
          hdr = hdr1
          dims = dims1
          nx = dims[0]
          ny = dims[1]
          xarr = findgen(nx)#replicate(1.0, ny)
          yarr = replicate(1.0, nx) # findgen(ny)
          ;; generate slitmask
          IF (keyword_set(slitfile)) THEN $
             tset_slits_orig = xmrdfits(slitfile, 1 $
                                        , silent = (keyword_set(verbose) EQ 0)) $
          ELSE message, 'Did not find slitfile'
          ;; Default to all mask values of 1B, which rejects points
          IF n_pixel GT 0 THEN BEGIN
             imgarr = make_array(dimension = [dims, n_pixel], /float)
             inmask = make_array(dimension = [dims, n_pixel], /byte $
                                 , value = 1B)
          ENDIF
          sxaddpar, hdr, 'NEXP', nfile, 'Number of exposures in this file', $
                    before = 'EXPTIME'
       endif else begin
          if (total(dims NE dims1) NE 0) then $
             message, 'Inconsistent dimensions for input images'
       endelse
       temp_norm_flat = fltarr(nx, ny) + 1.0D 
       tset_slits = tset_slits_orig

       ;; Shift slits to account for flexure. This will improve edge behavior
       IF strcmp(telescope, 'KeckI') OR KEYWORD_SET(NIRSPEC) THEN BEGIN
          ;; Aggressive edge masking for flats
          xshift = long_xcorr_slits(image1, tset_slits, /shift)
          TOL_EDG = TOL_EDG1 >  1.5*abs(XSHIFT)
       ENDIF ELSE IF KEYWORD_SET(SOFI) THEN TOL_EDG = TOL_EDG1 $
       ELSE IF KEYWORD_SET(LUCI) THEN TOL_EDG = TOL_EDG1 $
       ELSE IF KEYWORD_SET(ISAAC) THEN TOL_EDG = TOL_EDG1 $
       ELSE TOL_EDG = TOL_EDG1
       midwidth = tset_slits[1].coeff[0, *] - tset_slits[0].coeff[0, *]
       ximg = long_slits2x(tset_slits, edgmask = edgmask, TOL_EDG = TOL_EDG)
       slitmask = long_slits2mask(tset_slits, nslit = nslit) $
                  * (ximg GT 0. AND ximg LT 1.0)

       ;; If waveimg and piximg didn't get passed in, generate
       ;; them from the wavefile, or create dummy images
       IF KEYWORD_SET(waveimg) EQ 0 AND KEYWORD_SET(piximg) EQ 0 THEN BEGIN 
          IF KEYWORD_SET(wavefile) THEN BEGIN
             pixset  = xmrdfits(wavefile $
                                , silent = (keyword_set(verbose) EQ 0), 1) 
             wavesvfile =  repstr(wavefile, '.fits', '.sav')
             IF file_search(wavesvfile) THEN restore, wavesvfile
             piximg = long_wpix2image(pixset, tset_slits, XFIT = xfit $
                                      , waveimg = waveimg)
          ENDIF ELSE BEGIN
             splog, 'No wavelength image specified.'
             splog, 'Using row number for wavelength solution'
             ;; JXP Kludge  (25 Feb 2012)
             ;;piximg = replicate(1.0, nx)#findgen(ny)
            piximg = replicate(1.0, nx)#findgen(ny) + (findgen(nx)/float(nx-1)/10.)#replicate(1., ny)
            waveimg = piximg + 1000. 
         ENDELSE
       ENDIF
       FOR slitid = 1L, nslit DO BEGIN
          splog, 'Working on slit # ', slitid
          IF KEYWORD_SET(waveimg) THEN BEGIN
             igood = where(slitmask EQ slitid AND waveimg GT 0, ngood) 
             IF ngood EQ 0 THEN BEGIN
                splog,'Slit',slitid,' appears to have problematic wavelengths'
                igood = where(slitmask EQ slitid, ngood) 
             ENDIF
          ENDIF ELSE igood = where(slitmask EQ slitid, ngood)
          if (ngood GT 0) then begin
              ;; Compute the approximate number of pixels per row for this slit
              npercol = floor(float(ngood) / ny) > 1
              ;; Demand at least 10 pixels per row (on average) per degree
              ;; of the polynomial
              npoly = npoly1 < ceil(npercol / 10.)
              npoly = npoly > 1
              
              ;; Sort the pixels in wavelength 
              ii = igood[sort(piximg[igood])]
              ;; Always use the same breakpoints 
              ybkpt = 0
              if keyword_set(temppillum0) then begin
                 illum_set = xmrdfits(tempillum0, slitid, /silent)
                  ybkpt = illum_set.fullbkpt
              endif
              
              tempivar = invvar1[igood]
              ;; remove extreme edges as they seem to ruin the fits...
              ;xsort = sort(ximg[igood])
              ;tempivar[xsort[0:dims[1]*2]] = 0
              ;tempivar[xsort[ngood-dims[1]*2:*]] = 0
              long_flatfield_specillum, piximg[igood], ximg[igood] $
                                        , image1[igood], spec_set $
                                        , illum_set $
                                        , SLITSAMP=slitsamp $
                                        , invvar = tempivar $ 
                                        , slitwidth = midwidth[slitid-1] $
                                        , modelfit = modelfit $
                                        , FINEBKPT=kast $
                                        , ybkpt = ybkpt $
                                        , npoly = npoly, CHK = CHK $
                                        , PIXFIT = (use_pixel[ifile] EQ 1);, DEBUG=IGOOD
            ; Divide this flat by the B-spline model (for this slit)
            IF use_pixel[ifile] EQ 1 THEN BEGIN
                qgood = modelfit GT 0
                temp_norm_flat[igood] = qgood*image1[igood]/ $
                  (modelfit + (qgood EQ 0))
                index = long(total(use_pixel[0:ifile]))-1L
                imgarr[nx*ny*index + igood] = $
                  qgood * image1[igood] / (modelfit + (qgood EQ 0))
                inmask[nx*ny*index + igood] = $
                  (invvar1[igood] LE 0 OR qgood EQ 0) ; =0 for good
                ;CHK=1
                IF KEYWORD_SET(CHK) THEN BEGIN
                    xmin = min(xarr[igood])
                    xmax = max(xarr[igood])
                    ymin = min(yarr[igood])
                    ymax = max(yarr[igood])
                    edginds = WHERE(edgmask EQ 1, nedg)
                    qa_temp = temp_norm_flat
                    IF nedg GT 0 THEN qa_temp[edginds] = 1.0
                    qa_temp = qa_temp*float(slitmask EQ slitid)
                    qaimg = qa_temp[xmin:xmax, ymin:ymax]
                    xatv, qaimg, wv = waveimg, min = 0.9, max = 1.1, /block
                    qa_temp = 0
                ENDIF
            ENDIF
            mwrfits, illum_set, temp_illum, create = (slitid EQ 1)
            mwrfits, spec_set, temp_spec, create = (slitid EQ 1)
         endif
    endfor ;; End loop over slits for this image
    ;; write out the temporary flats
    IF KEYWORD_SET(use_pixel[ifile]) THEN BEGIN
        splog, 'Storing pixel dependence in ', temp_pixel
        mwrfits, temp_norm_flat, temp_pixel, /create
    ENDIF
    illum_now = long_slitillum(temp_illum, slitmask, ximg, edgmask)
    spawn, '/bin/cp -f ' + temp_illum +  '  ' +  temp_illum + '.tmp'
    mwrfits, illum_now, temp_illum, /create
    FOR slitid = 1L, nslit DO BEGIN
        illum_now = xmrdfits(temp_illum + '.tmp', slitid, /silent)
        mwrfits, illum_now, temp_illum
    ENDFOR
    spawn, '/bin/rm -f ' + temp_illum + '.tmp'
    sxaddpar, hdr, string(ifile+1, format = '("FILE",i2.2)'), $
              fileandpath(filenames[ifile]), $
              ' File number ' + strtrim(string(ifile)), before = 'EXPTIME'
endfor   ;; End loop over files

   ;----------
   ; Average the stack of images, with outlier-rejection
   if (n_pixel GT 0) then begin
       if size(imgarr,/n_dimensions) EQ 3 then $
         imgfinal = djs_avsigclip(imgarr, 3, sigrej=sigrej_pixel, $
                                  maxiter=maxiter, inmask=inmask) $
       else imgfinal=imgarr
       ;; Now be careful to assign the pixels on the slit edges and 
       ;; and the pixels in between the slits to flat field values of 1.0
       unitinds = WHERE(slitmask LE 0.0 OR EDGMASK, nunit)
       IF nunit GT 0 THEN imgfinal[unitinds] = 1.0
      ; Write output file
      mwrfits, imgfinal, superpixflatfile, hdr, /create
  endif

   ;----------
   ; Construct the illumination flat

   if (n_illum GT 0) then begin
       temp_ifinal = djs_filepath('tempifinal-'+superillumflatfile $
                                  , root_dir = tempdir)
       for slitid = 1L, nslit do begin
         xspace = 0.1 /midwidth[slitid-1]
         ntmp = long(1./xspace)+1
         xtmp = findgen(ntmp)*xspace
         illum_array = fltarr(ntmp, n_illum)
         illum_inds = WHERE(use_illum)
         for ij = 0L, n_illum-1L do begin
            fname = fileandpath(filenames[illum_inds[ij]])
            temp_illum = djs_filepath('tempillum-'+fname, $
             root_dir=tempdir)
            illum_set = mrdfits(temp_illum, slitid, /silent)
            if (keyword_set(illum_set)) then $
              illum_array[*, ij] = bspline_valu(xtmp, illum_set) $
            else $
              splog, 'Warning: Unable to read illumination for file ', $
              fname, ' slit #', slitid
         endfor

         if (n_illum EQ 1) then $
          illum_comb = illum_array $
         else $
          illum_comb = djs_avsigclip(illum_array, 2, sigrej = sigrej_illum  $
                                     , maxiter = maxiter)
         illum_comb_set = illum_set
         ec = bspline_fit(xtmp, illum_comb, illum_comb*0+1, illum_comb_set)
         mwrfits, illum_comb_set, temp_ifinal, create = (slitid EQ 1L) $
                  , /silent
     endfor
     illum = long_slitillum(temp_ifinal, slitmask, ximg, edgmask)
     mwrfits, illum, superillumflatfile, /create
     FOR slitid = 1L, nslit do begin
         illum_set = xmrdfits(temp_ifinal, slitid, /silent)
         mwrfits, illum_set, superillumflatfile
     ENDFOR
 endif
 
   
   splog, 'Elapsed time = ', systime(1)-t0, ' sec'

   return
end
;------------------------------------------------------------------------------
