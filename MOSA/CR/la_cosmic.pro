;------------------------------------------------------------------------------
;+
; NAME:
;    reckon_statsec
;
; PURPOSE:
;   parse the statsec string in la_cosmic and return the indices
;   which bracket the section to be exaimed.
;
;
; CALLING SEQUENCE:
;     reckon_statsec,  statstring,arrsize
;
; INPUTS:
;   statstring:  The string such as "[23:34,*]"
;   arrsize   :  Array with the x and y sizes of the image
;
; OUTPUTS:
;   returns the indices of the statsec in the form [x1,x2,y1,y2]
;
;
; PROCEDURES CALLED:
;
; COMMENTS:
;   A good deal of error checking is done to ensure that
;   a statsection will be valid.
;
; NOTES:
; BUGS:
;
; REVISION HISTORY:
;   20-May-2001  Written by Joshua Bloom, Caltech (jsb@astro.caltech.edu)
;-
;------------------------------------------------------------------------------
function reckon_statsec, statstring,arrsize
;; must be of the form
;; X1:X2,Y1:Y2

tmp = size(statstring)

if (tmp[1] ne 7) then begin
    print, 'RECKON_STATSEC: Warning, statsec not valid, using full image'
    return, [0,arrsize[0]-1,0,arrsize[1]-1]
endif
tmp = strsplit(statstring,'[',/extract)
tmp = strsplit(tmp[0],']',/extract)

;; break up the string by the comma and evaluate
str = strsplit(tmp[0],',',/extract)
nstr = n_elements(str)
if (nstr ne 2) then begin
    print, 'RECKON_STATSEC: Warning, statsec not valid, using full image'
    return, [0,arrsize[0]-1,0,arrsize[1]-1]
endif
retarr = lonarr(4)
for i=0,1 do begin
    ; now look at each string and parse
    str1 = strsplit(str[i],':',/extract)
    nstr1 = n_elements(str1)
    if (nstr1 gt 2) then begin
        ;; malformed strsep
        retarr[i*2] = 0
        retarr[i*2 + 1] = arrsize[i] - 1
    endif
    if (nstr1 eq 1) then begin
        
        if (stregex(str1[0],'[*]',/boolean)) then begin
            ;; the user wants the entire image
            retarr[i*2] = 0
            retarr[i*2 + 1] = arrsize[i] - 1
        endif else begin
            ;; it's a number, so convert it 
            retarr[i*2] = long(str1[0])
            retarr[i*2 + 1] = long(str1[0])
        endelse
    endif else begin
        retarr[i*2] = long(str1[0])
        retarr[i*2 + 1] = long(str1[1])
    endelse
endfor

return, retarr
end

;------------------------------------------------------------------------------
;+
; NAME:
;   lacos_replace, arr, repval, low, high 
;
; PURPOSE:
;  replace pixels whose value are between low and high with value = repval
;  modelled after IRAF (!) IMREPLACE
;
; CALLING SEQUENCE:
;     lacos_replace, arr, repval, low, high
;
; INPUTS:
;   arr:       number array of any size or dimension.
;   repval     valid replacment value
;   low,high   bracket values to replace (can be 'INDEF') to replace all
;
; OUTPUTS:
;   returns the array = arr but with the replaced pixels
;
; PROCEDURES CALLED:
;
; COMMENTS:

; NOTES:
; BUGS:
;
; REVISION HISTORY:
;   20-May-2001  Written by Joshua Bloom, Caltech (jsb@astro.caltech.edu)
;-
;------------------------------------------------------------------------------
function lacos_replace, arr, repval, low, high
extr = 1d40
slow  = size(low)
shigh = size(high)

if (shigh[1] eq 7) then begin
    if (high ne 'INDEF') then begin
        print, 'LACOS_REPLACE: Sorry must call with INDEF not'
        print, high
        return, arr
    endif
    high = extr
endif

if (slow[1] eq 7) then begin
    if (low ne 'INDEF') then begin
        print, 'LACOS_REPLACE: Sorry must call with INDEF not'
        print, slow
        return, arr
    endif
    low = -1d0 * extr
endif

high = double(high)
low  = double(low)
bads = where((arr le high) and (arr ge low), n)
if (n ge 1) then arr[bads] = repval
return, arr
end

;------------------------------------------------------------------------------

;------------------------------------------------------------------------------
;+
; NAME:
;   la_cosmic
;
; PURPOSE:
;   Remove cosmic rays from imaging data.  Images must be debiased for
;     gain estimation to work properly.
;
; CALLING SEQUENCE:
;   la_cosmic, imlist, [outlist=, masklist=, sigclip=, gain=, readn=, $
;               skyval=,objlim=, niter=,sigfrac=,verbose=,statsec=, $
;               zeroindexed=,masksuffix=, outsuff=,isbig=,blocksize=]
;
; INPUTS:
;   imlist:     List (strarr) of images to be cleaned or string with
;                regexp of files      
;
; OPTIONAL INPUTS:
;   outlist:    List (string array) of output cleaned images
;   masklist:   List of mask files associated with the cleaned images
;   sigclip     Level of cr clipping
;   gain        Gain of the CCD (< 0 to estimate) (default=-1.0)
;   readn       Readnoise of CCD (default=0.0)
;   skyval      Avg. of sky level already subtracted off the images (array)
;   objlim      Detection threshold for objects (not crs)
;   niter       How many passes to take
;   sigfrac     Sigfrac: see PASP paper
;   verbose     Verbose printing of steps? 
;   statsec     image region for statistics
;                string which can be of the form "[23:45,100:300]" or
;                                                 [*,100:400], etc.
;   zeroindexed Is the image region zeroindexed?
;   masksuffix  Suffix to automatically determine mask output name
;   outsuff     Suffix to automatically determine output image name
;   blocksize   Size of working blocks.  Keep in integer multiples of 512
;   isbig       Tell the routine to chop up the images in to more
;                manageable sections of size blocksize x blocksize
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES CALLED:
;   DJS_ITERSTAT
;   DJS_MEDIAN   
;   reckon_statsec
;   lacos_replace
;   astro-library 
;
; COMMENTS:
;  This routine is based on Pieter Van Dokkum's "LACOSMIC" Cosmic
;  ray rejection routine for IRAF.
;
;  If you find that after ~4 iterations that the routine is still
;  finding up cosmic rays chances are that you are digging into the
;  noise and or objects.  Try setting the sigclip number a bit higher
;
; DEFAULTS:
;   outlist:    This will be set to the input list with the suffix
;                outsuffix + '.fits' 
;   masklist:   This will be set to the input list with the suffix
;                masksuffix + '.fits' 
;   sigclip     4.5
;   gain        -1.0 (e-/DN) (forces routine to estimate for gain)
;   readn       0.0
;   skyval      0.0 
;   objlim      4.0
;   niter       4
;   sigfrac     0.5
;   verbose     1 (yes)
;   statsec     "[*,*]" (use the entire image to estimate the gain)
;   zeroindexed 0 (no)
;   masksuffix  "-mask"
;   outsuff     "-out"
;   isbig       0
;   blocksize   1024
;
; NOTES:
;    (1) This routine will only work on .fits images 
;    (2) Haven't checked how useful it is on spectroscopic images.
;    (3) Speed is clearly an issue.  It takes about 21 seconds per iteration
;        on a fairly zippy (650 PIII Linux 2.4 kernel) laptop on a 1k x 1k image.
;        Using the isbig=1 the speed should scale as the number of
;        pixels. So the scaling for speed is,
; 
;                t = 21 s * (n/1024)^2 * (cpu/650 MHz)^(-1)
;
;        So that a 2k x 2k image will take ~80s per interation.


; EXAMPLES:
;   1. Remove the cosmic rays from a debiased flat-fielded HST STIS/Clear image
;   hst1.fits and created a replaced image called hst1-cleaned.fits
;   with a mask file called mask.fits. Set the read noise to be 4.46 e-
;   and the gain to be = 1.0 (e-/DN)
;
;      IDL> la_cosmic, ['hst1.fits'], masklist=['mask.fits'], $
;                      outsuffix="-cleaned", readn=4.46, gain=1.0
;
;   2. Remove the cosmic rays from all images starting with the name
;   hst and create masks "hst*-mask.fits" and output images
;   "hst*-out.fits". Set sigclip = 4. Let la_cosmic determine the gain
;   from a region of [50:100,*] (indexed exactly as in IRAF, ie. unity
;   indexed).
;
;      IDL> la_cosmic, 'hst*.fits', outsuffix="-out",
;                masksuffix="-mask",statsec="[50:100,*]",zeroindexed=0, 
;                gain = -1, sigclip = 4.0
;
; BUGS:
;
;  1. If the image has not been debiased, then the gain estimation
;     will go horribly wrong.  Could add a "biassec" parameter to
;     allow the bias to be estimated.
;  2. Speed scaling of the routine works well until the point that
;     the images are too large to store entirely in memory.  Could write
;     a little section which chops up the image in to manageable
;     chunks and pastes those chuncks together...
;
; REVISION HISTORY:
;   20-May-2001  Written by Joshua Bloom, Caltech (jsb@astro.caltech.edu)
;-
;------------------------------------------------------------------------------
pro la_cosmic, imlist, outlist=outlist, masklist=masklist, sigclip=sigclip, $
               gain=gain, $
               readn=readn, indx=indx, $
               skyval=skyval,objlim=objlim, niter=niter,sigfrac=sigfrac, $
               verbose=verbose,statsec=statsec,zeroindexed=zeroindexed,$
               masksuffix=masksuff, outsuff=outsuff, isbig=isbig,blocksize=blocksize

;; set some sensible defaults
if ( (not keyword_set(outlist)) and (not keyword_set(outsuff))) then begin
    outsuff = "-out"
endif

if ( (not keyword_set(masklist)) and (not keyword_set(masksuff))) then begin
    masksuff = "-mask"
endif
nfiles = n_elements(imlist)
if (nfiles eq 1) then begin
    ;; nice to allow user to find a whole bunch of files
    ;; which match a particular describtion
    imlist = findfile(imlist[0])
endif
nfiles = n_elements(imlist)
if not keyword_set(indx) then indx = 0
if (not keyword_set(outlist)) then outlist = imlist
if (not keyword_set(masklist)) then masklist = imlist
if (not keyword_set(gain)) then gain       = -1.0
if (not keyword_set(readn)) then readn     = 0.0
if (not keyword_set(skyval)) then begin    
    skyval   = fltarr(nfiles)
endif
if (not keyword_set(sigclip)) then sigclip = 4.5
if (not keyword_set(sigfrac)) then sigfrac = 0.5
if (not keyword_set(objlim)) then objlim   = 3.0
if (not keyword_set(niter)) then niter     = 4
if (not keyword_set(verbose)) then verbose = 0
if ((not keyword_set(isbig)) and (not keyword_set(blocksize))) then begin
    isbig = 0
endif
if keyword_set(blocksize) then begin
    isbig = 1
endif
if (not keyword_set(blocksize)) then begin
    blocksize = 512l
endif
if (not keyword_set(statsec)) then statsec = "[*,*]"
if (not keyword_set(zeroindexed)) then zeroindexed = 1

gain = double(gain)
readn = double(readn) > 0d0

if (nfiles le 0) then begin
    print, "LACOSMIC: Sorry there's not much to do here. Returning."
    return
endif
if ((nfiles ne n_elements(outlist) or (nfiles ne n_elements(masklist)))) then $
  begin
    print, "LACOSMIC: Sorry the number of elements in outlist, masklist, &" + $
           " imlist must"
    print, "      match up.  Leave masklist and outlist blank to " + $
           "derive the outputted"
    print,        "      names from the input names"
    return
endif

if (verbose) then begin
    print, ""
    print, "----------------------------------------------------"
    print, ""
    print, " L.A. Cosmic: Laplacian cosmic ray removal"
    print, ""
    print, "      by Pieter van Dokkum"
    print, "    IDL version by Josh Bloom"
    print, ""
    print, "   Imaging version 1.0 (April 2001) "
    print, "----------------------------------------------------"
endif
;; make the kernel as an array of 3x3
lakernel = [[0.0, -1.0, 0.0],[-1.0,4.0,-1.0],[0.0,-1.0,0.0]]
gkernel  = [[1,1,1],[1,1,1],[1,1,1]]

for i=0, nfiles -1 do begin
    filename = imlist[i]
    if (verbose) then begin
        print, "LACOSMIC: Working on image " + filename
    endif
;    fxread, filename, oldoutput, oldheader
    ;; JXP
    oldoutput = xmrdfits(filename, indx, oldheader, /silent)
    xsize = long(n_elements(oldoutput[*,0]))
    ysize = long(n_elements(oldoutput[0,*]))
    outmask = bytarr(xsize,ysize)

    usegain = gain
    sstop = 0
    iter = 1
    ;if (skyval[i] gt 0.0) then oldoutput = temporary(oldoutput) + skyval[i]
    xblock = blocksize ; 512l
    yblock = blocksize ; 512l

    if (isbig) then begin
        ;; we may have to do some chopping here.  Chop the array into
        ;; blocks
        if ((yblock GT ysize) and (xsize GT xblock)) then begin
            ;;; there's really no need to chop up here.
            isbig = 0
            potx = 1
            poty = 1
            regx = intarr(1,2)
            regy = intarr(1,2)
            regx[0,0] = 0l
            regx[0,1] = xsize-1l
            regy[0,0] = 0l
            regy[0,1] = ysize-1l
        endif else begin
            potx = long(xsize) / long(xblock)
            poty = long(ysize) / long(yblock)
            ;; create an array of indices which mark off the regions
            ;; of smaller blocks
            regx = intarr(potx,2)
            regy = intarr(poty,2)
            regx[0,0] = 0
            regx[0,1] = xblock - 1
            regy[0,0] = 0
            regy[0,1] = yblock - 1

            for k=1,potx-2 do begin
                regx[k,0] = regx[k-1,1] + 1
                regx[k,1] = regx[k,0] + xblock - 1
            endfor
            ;; the last block may be bigger than the others
            regx[potx-1,0] = regx[potx-2,1] + 1
            regx[potx-1,1] = xsize - 1
            for k=1,poty-2 do begin
                regy[k,0] = regy[k-1,1] + 1
                regy[k,1] = regy[k,0] + xblock - 1
            endfor
            regy[poty-1,0] = regy[poty-2,1] + 1
            regy[poty-1,1] = ysize - 1
        endelse
    endif else begin
                    ;;; there's really no need to chop up here.
        isbig = 0
        potx = 1
        poty = 1
        regx = intarr(1,2)
        regy = intarr(1,2)
        regx[0,0] = 0l
        regx[0,1] = xsize-1l
        regy[0,0] = 0l
        regy[0,1] = ysize-1l

    endelse

    while(not sstop) do begin
        if (verbose) then begin
            print, "-------------------Iteration" + $
                   strtrim(string(iter),2) + "------------------------"
        endif
        ;; add back in the background if the user so desires
        ;oldoutput = oldoutput + skyval[i]
        if (skyval[i] gt 0.0) then oldoutput = temporary(oldoutput) + skyval[i]
        if (gain le 0.0) then begin
            if (verbose and (iter eq 1)) then print, $
              "Trying to determine gain automatically: "
            if (verbose and (iter gt 1)) then print, $
              "Improving gain estimate: "
            ;; figure out what statsection to use from the statsec
            ;; string
            arr=reckon_statsec(statsec,[xsize,ysize])
            if (zeroindexed eq 0) then begin
                ;; user gave a statsec region like IRAF.  Unity indexed..
                arr = arr - 1
            endif
            arr[1] = (arr[0] > arr[1]) < (xsize - 1)
            arr[3] = (arr[2] > arr[3]) < (ysize - 1)
            djs_iterstat, oldoutput[arr[0]:arr[1],arr[2]:arr[3]], $
                          sigrej=5.0, maxiter=10.0,mean=estmean,$
                          median=estmedian, sigma=estsigma
            skylev = estmedian
            ;sigima = abs(oldoutput[arr[0]:arr[1],arr[2]:arr[3]] - $
            ;  median(oldoutput[arr[0]:arr[1],arr[2]:arr[3]],7,/even))
            sigima = abs(oldoutput[arr[0]:arr[1],arr[2]:arr[3]] - $
                         djs_median(oldoutput[arr[0]:arr[1],arr[2]:arr[3]],$
                                    width=7,boundary='reflect'))
            djs_iterstat, sigima, $
                          sigrej=5.0, maxiter=10.0,mean=estmean,$
                          median=estmedian, sigma=estsigma
            sig = estmedian * 1.48
            usegain = skylev/sig^2
            if (verbose) then begin
                print, "  Approximate sky level = " + $
                       strtrim(string(skylev),2) + ' ADU'
                print, "  Sigma of sky = " + strtrim(string(sig),2)
                print, "  Estimated gain = " + strtrim(string(usegain),2)
            endif
            if (usegain le 0) then begin
                print, 'LACOSMIC: error.  Gain was found to be less than zero'
                print, '  is it possible you forgot to give a "skyval"?'
                return
            endif
        endif
        if (verbose) then begin
            print, 'Convolving image with Laplacian kernel'
            print, ' '
        endif
        ;; we've got to chop this image up
        nchop = 1
        ncosmicray = 0
        for xchop = 0, potx - 1 do begin
        for ychop = 0, poty - 1 do begin
            oldoutputwork = oldoutput[regx[xchop,0]:regx[xchop,1], $
                                      regy[ychop,0]:regy[ychop,1]]
            outmaskwork = outmask[regx[xchop,0]:regx[xchop,1], $
                                      regy[ychop,0]:regy[ychop,1]]
            if (verbose) then begin
                print, 'Working on block #' + strtrim(string(nchop),2) + $
                       ' out of ' +  strtrim(string(potx*poty),2)
            endif

        ;; rebin the image and convolve with the kernel then rebin
            tmpxsize = regx[xchop,1] - regx[xchop,0] + 1
            tmpysize = regy[ychop,1] - regy[ychop,0] + 1

        im2 = rebin( (convol(rebin(oldoutputwork,2*tmpxsize, $
                           2*tmpysize,/sample), $
                  lakernel,1,/center,/edge_truncate) > 0.0), $
                    tmpxsize, tmpysize)

        ;; compute the 5x5 median to remove all the cosmic rays
        med5 = djs_median(oldoutputwork,width=5,boundary='reflect')
        ;med5 = median(oldoutput,5,/even)
        ;; NOTE: the boundary is not handled the same as in IRAF
        ;; This affects the outer 2 pixel boundary...could change
        ;; with some kludges...

;        SAME but slower: med5 = lacos_replace(med5,0.0001,'INDEF', 0.0)
        bad = med5 LE 0.0
        med5 = 0.0001 * bad + temporary(med5) * (1.0 - bad)
        bad = [0]

        ;; create a noise model based on this median image knowing the
        ;; gain, and readnoise of the image 
        ;; note that this step supposes that no sky background subtration
        ;; has been done
        if (verbose) then begin
            print, 'Creating noise model using:'
            print, '  gain = ' + strtrim(string(usegain),2)
            print, '  readnoise = ' + strtrim(string(readn),2)
        endif
        noise  = sqrt(med5*usegain + readn^2)/usegain
        med5 = [0]
        sigmap =  im2/noise/2.d 
        ;; free up some memory
        im2 = [0]
        sigmap =  -1.0 * (djs_median(sigmap,width=5,boundary='reflect') - $
                          temporary(sigmap))
        
        ;; do a replacement setting all the high values to 1
        firstsel = sigmap
        firstsel = firstsel * (firstsel GT sigclip)
        firstsel = firstsel * (firstsel LT 0.1) + (firstsel GE 0.1)
        ;med3 = median(oldoutput,3,/even)
        ;med7 = median(med3,7,/even)
        med3 = djs_median(oldoutputwork,width=3,boundary='reflect')
        med7 = djs_median(med3,width=7,boundary='reflect')

        med3 = (temporary(med3) - med7)/noise
        noise = [0]
        med7 = [0]
        med3 = temporary(med3) > 0.01

        starreject = firstsel*sigmap/med3
        med3 = [0]
        firstsel = temporary(firstsel) * (starreject GT objlim)
        starreject = [0]

        ;; grow CRs by one pixel and check in original sigma map
        gfirstsel = convol(firstsel,gkernel,/center,/edge_truncate)
        firstsel = [0]
        gfirstsel = sigmap * ((gfirstsel GT 0.5) + $
                             gfirstsel*(gfirstsel LE 0.5))
        gfirstsel = (temporary(gfirstsel) GT sigclip)
        sigcliplow = sigfrac * sigclip
        
        finalsel = convol(gfirstsel,gkernel,/center,/edge_truncate)
        finalsel = sigmap * ((finalsel GT 0.5) + $
                             finalsel*(finalsel LE 0.5))
        sigmap = [0]
        finalsel = (temporary(finalsel) GT sigcliplow)
        
        ;; how many cosmic rays were found here
        gfirstsel = (1.0 - outmaskwork)*finalsel
        ttt = where(gfirstsel ge 0.5,npix)
        ttt = [0]
        gfirstsel = [0]
        outmaskwork = (temporary(outmaskwork) + finalsel) < 1
        finalsel = [0]

        ncosmicray = ncosmicray + npix
        inputmask = (1.0 - 10000.0*outmaskwork)*oldoutputwork
        inputmask = lacos_replace(inputmask, !VALUES.F_NAN,'INDEF',-9999)
        ;med5 = median(inputmask,5,/even)*outmask
        med5 = djs_median(inputmask,width=5,boundary='reflect')*outmaskwork
        inputmask = [0]
        output = (1.0 - outmaskwork)*oldoutputwork + med5
        med5 = [0]
        oldoutput[regx[xchop,0]:regx[xchop,1],regy[ychop,0]:regy[ychop,1]] = $
                  output[0:(tmpxsize-1),0:(tmpysize-1)]
        outmask[regx[xchop,0]:regx[xchop,1],regy[ychop,0]:regy[ychop,1]] = $
                  outmaskwork[0:(tmpxsize-1),0:(tmpysize-1)]
        nchop = nchop + 1
        endfor
        endfor

        print, 'Found ' + strtrim(string(ncosmicray)) + ' cosmic rays in iteration ' + strtrim(string(iter))
        if (ncosmicray eq 0) then sstop = 1
        iter = iter + 1
        if (iter gt niter) then sstop = 1
        if (skyval[i] > 0) then oldoutput = temporary(oldoutput) - skyval[i]
        ;; groovy ---------------------------------------

    endwhile
    ;; save the output image and the outmask
    troot = strsplit(imlist[i],".fits",/extract,/regex)

    tmp = strsplit(outlist[i],".fits",/extract,/regex)
    if (outlist[i] eq imlist[i]) then begin
        outname = tmp[0] + outsuff + ".fits"
    endif else begin
        outname = tmp[0] + '.fits'
    endelse
    tmp = strsplit(masklist[i],".fits",/extract,/regex)
    if (masklist[i] eq imlist[i]) then begin
        maskname = tmp[0] + masksuff + ".fits"
    endif else begin
        maskname = tmp[0] + '.fits'
    endelse
    print, 'Writing mask   image = ' + maskname
;    print, '        output image = ' + outname
    fxhmake,h,outmask
    fxwrite,maskname,h,outmask

;    fxhmake,h,oldoutput
;    fxwrite,outname,h,oldoutput
endfor
return
end






