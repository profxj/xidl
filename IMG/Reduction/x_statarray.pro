;+ 
; NAME:
; x_statarray
;     Version 1.1
;
; PURPOSE:
;    This set of routines takes a series of flats observed through the
;    diffuser and creates a normalized Flat used to correct
;    pixel-to-pixel response variations.  A principal challenge with
;    MIKE in its current state (April 2004) is that it is difficult to
;    get sufficient counts on the blue side.  
;
;    mike_mkmflat :: The main routine simply does some basic
;    accounting, organizes bias subtraction and performs I/O.
;
;    mike_mkflat_work :: Accepts the name(s) of a series of OV
;    subtracted milky flats.  It then:
;      1.  Opens the file and creates an inverse variance array
;      2.  Takes out the low order variation in the image (lo_interp)
;      3.  Performs a series of medians along the columns rejecting bad
;          pixels (replace by local mean). 
;      4.  If multiple images were input, they are stacked with
;          rejection and the final image is returned
;
; CALLING SEQUENCE:
;   
;  mike_mkmflat, mike, setup, [side]
;
; INPUTS:
;   mike     -  MIKE structure
;   setup    -  Setup identifier 
;   [side]   -  Blue (1), Red (2), or both [1,2L]    (Default: [1,2L])
;
; RETURNS:
;
; OUTPUTS:
;  One normalized flat per setup per side with names like
;  'Flats/Flat_B_01_M.fits.gz' 
;
; OPTIONAL KEYWORDS:
;   /CLOBBER - Overwrite Output MilkyFlat
;   /OVCLOB  - Overwrite OV files if they exist for the flats
;   /SVOV    - Save the OV files created during this step
;   /USEBIAS - Use bias frame in OV subtraction
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   9-May-2005 Adapted by JXP from existing programs by SB
;-
;------------------------------------------------------------------------------
pro x_statarray, image, iaxis, mean=mean, stddev=stddev, $
       sigrej=sigrej, niter=niter, OVRIDE=ovride, NOREJCOLM=norejcolm

     if N_PARAMS() LT 2 then begin
        print, 'Please give an axis to calculate'
        print, 'Use statarray, statarray, image, 0, mean=mn.... for full array'
        return
     endif


     ndim = (size(image))[0] 
     if ndim NE 3 then stop ;; Requires an array of 2D images
     if ndim LT iaxis then begin
        print, "iaxis is larger than ndim of image"
        return
     endif
     if NOT keyword_set(sigrej) then sigrej=3.0
     if NOT keyword_set(niter) then niter=1

     ;; Size
     sz = size(image, /dimensions)
     if iaxis NE ndim then stop  ;; Necessary for memory issues

     nbefore = 1L
     nafter  = 1L
     if (iaxis EQ 0) OR (ndim EQ 1) then nsum = n_elements(image)  $
     else begin
       nsum = (size(image))[iaxis]
       for i=1,iaxis-1 do    nbefore = nbefore * (size(image))[i]
       for i=iaxis+1,ndim do nafter  = nafter * (size(image))[i]
     endelse

     ; first iteration
     mn1 = djs_median(image,iaxis)
     
     mnfull = transpose(reform((mn1[*]) # $
                replicate(1,nsum),nbefore,nafter,nsum),[0,2,1])
     stddev = sqrt(djs_median((image - mnfull)^2,iaxis))

     
     ;; Deal with zero stddev points (rare)
     a = where(stddev EQ 0., na)
     if na NE 0 then begin
         for ii=0L,na-1 do $
           stddev[a[ii]] = mean( stddev[a[ii]-lindgen(5)]+ $
                                 stddev[a[ii]+lindgen(5)] )
     endif

     if not keyword_set( NOREJCOLM ) then begin
         sm = sqrt(median(djs_avsigclip(stddev,2), 5) # $
                   median(djs_avsigclip(stddev,1), 5)  > 0)
         badmask = stddev GT 5.0*sm
         totc = total(badmask,2)
         badc = where(float(totc)/sz[1] GT 0.3, complement=goodc)
         badmask[goodc,*] = 0
     endif

     ;; More aggressively..
     bad = where(stddev EQ 0., nbad)
     if nbad NE 0 then begin
         if keyword_set(OVRIDE) then stddev[bad] = median(stddev) else stop
     endif

     for iter=1,niter do begin
         maskfull = bytarr(sz)
;        sdfull = transpose(reform((stddev[*]) # $
;        replicate(1,nsum),nbefore,nafter,nsum),[0,2,1])
         ;; Loop on images
         for jj=0L,nsum-1 do begin
             maskfull[*,*,jj] = (image[*,*,jj] GT mnfull[*,jj] $
                                 - sigrej*stddev) AND $
               (image[*,*,jj] LT mnfull[*,jj] + sigrej *stddev)
         endfor
         
         norm = total(maskfull, iaxis)

         mean = total(image*maskfull,iaxis) /(norm + (norm EQ 0)) * (norm GT 0)
         mnfull = transpose(reform((mean[*]) # $
                 replicate(1,nsum),nbefore,nafter,nsum),[0,2,1])
         stddev = sqrt(total((image - mnfull)^2*maskfull,iaxis) / $
                       (norm - (norm NE 1))) * (norm GT 1 AND badmask EQ 0)
    endfor

return
end
