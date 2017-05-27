;+
; NAME:
;    quicklook
;
; PURPOSE:
;    Examine one frame of DEIMOS science data in 'quicklook' manner
;    for analysis at the telescope.   subset of slitlets for
;    SN, studies alignment boxes for A-band test
;
; CALLING SEQUENCE:
;    quicklook, mask, inputfile
;
; INPUTS:
;    mask -- integer name of mask (e.g. 1145)
;    inputfile -- name of science frame to analyse
;
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;    outdir -- directory for output of data (current if not specified)
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
;    generates several plots ??.
;    generate calibSlit files first, then call this routine! 
;    Input data comes from current directory.
;    Output data goes into current directory.
;
; REVISION HISTORY:
;   MD, 6Aug2002
;
;----------------------------------------------------------------------
pro quicklook,  mask,  inputfile,  outdir=outdir, flat=flat


  if n_elements(flat) eq 0 then flat = 0

  if NOT keyword_set(outdir) then outdir = ''
  chiplist = [1, 2, 3, 4, 5, 6, 7, 8] ;use all chips

; -------- Get environment variables, set paths
  deimos_data = getenv('DEIMOS_DATA')+'/'
  calib_data = getenv('CALIB_DATA')+'/'

;  if deimos_data eq '/' then message, 'You need to set $DEIMOS_DATA!'
  if deimos_data eq '/' then deimos_data = ''
  if calib_data eq '/' then message, 'You need to set $CALIB_DATA!'

  deep_dir = getenv('DEEP_DIR')+'/'
  if deep_dir eq '/' then message, 'You need to set $DEEP_DIR!'
  badmaskname = calib_data+'deimos_badmask.fits.Z'
;  pixflatname = calib_data+'deimos_pix_mult_flat.fits'


  outdatadir = outdir
  masks = string(mask, format='(i4)')
  calib_flist = findfile(outdatadir+'calibSlit.'+masks+'*', count=filecount)
  
  if filecount EQ 0 then $
     message, 'No calibSlit files for this mask: Run ql_calibrate first!'

; get chip list for slitlets
  chip_calib = intarr(filecount)
  for i=0, filecount -1 do $
    chip_calib[i] = sxpar(headfits(calib_flist[i],  ext=1), 'CHIPNO')  


; process this frame
    readnoise = 2.32
    objectcat = deimos_tables(inputfile, bluslits=slitcoords, $
                            slitobjmap=slitobjmap) 
    sci_header = headfits(inputfile) ; primary FITS header

  for i=0, n_elements(chiplist) -1  do begin 
     chipno = chiplist[i]
  

     jj = where(chip_calib eq chipno, nff)
     if nff gt 0 then begin 
       print
       print, 'Chipno: ', chipno
       print, 'Number of slitlets: ', nff

       fflist = calib_flist[jj]
; -------- read subimage
;       print, 'Reading bad mask: ', badmaskname
       badchip = mrdfits(badmaskname, chipno,  /silent) 
;       deimos_pixflat_mul = deimos_badchip*0.+1.

; process the slitlets on this chip
;      deimos_spslit, chipno, maskno, fflist, $
;         deimos_badchip, sciencenames, $
;         outdatadir=outdatadir, pixmap=deimos_pixflat_mul, /flat 

      specimage = deimos_read_chip(inputfile, chipno)  
;     specimage = specimage*pixmap
;
; get ivar
      medspec = djs_median(specimage, width=5, boundary='reflect') > 0.01 ;get median so as to flag CR's
;  use median for invvar unless CR or other anomaly makes it locally high
      specivar = (badchip eq 0)*((specimage/medspec gt 3.)/$
                         ((specimage > 5.) +readnoise^2) + $
       (specimage/medspec le 3.)/((medspec > 5.)+readnoise^2))
 ;     specivar = specivar/pixmap^2
      if chipno EQ 5 then begin 
         med = median(specimage[*, 4000:4040])
        print, 'Chip 5 -- fixing up image top with median ', med
         for k=4020, 4095 do specimage[*, k] = $
            specimage[*, k]-median(specimage[*, k])+med
         specimage[*, 4090:4095] = med
         specivar[*, 4090:4095] = 0
      endif 

     for islit=0, nff-1 do begin ;loop over all slitlets
     
        calib = mrdfits(fflist[islit], 1, header, /silent)
        slitn = sxpar(header, 'SLITNO')
        slitwid = sxpar(header, 'SLITWID')
        x0 = calib.x0 
        x1 = calib.x1
        rect_spec = deimos_rectify_slit(specimage, specivar, x0, x1, $
                 /interp,  npad=0, mask=specmask)
        rect_spec = rect_spec*specmask
        rect_specivar = deimos_rectify_slit(specivar, specivar, $
                                            x0, x1, /interp, npad=0)
        rect_specivar = rect_specivar*specmask

       interp_mask = deimos_rectify_slit(float(badchip EQ 1B), specivar*0+1., x0, x1, /interp, npad=0) NE 0

        vig_mask = deimos_rectify_slit(float(badchip AND 2B), specivar*0+1., x0, x1, /interp, npad=0) EQ 0

    if keyword_set(flat)  then begin 
           rect_spec = deimos_applyflat(rect_spec, calib.flat, $
              calib.slitfn, invvar=rect_specivar)
;deimos_flatfield(rect_spec, calib.flat, flat2d, $
;                 flat1d=slitfn, invvar=rect_specivar, /twod, mask=goodflat)
    endif


; remove this?

        sizey = (size(rect_spec, /dimens))[1]
        ncol  = (size(rect_spec, /dimens))[0]

; no flat applied??
        wave = lambda_eval(calib) ;extract lambda from structure
;
;
; data is now read in; compute SN, extract spectrum in case of stars
;  then run abandtest
;

        bitmask = (1B-calib.mask) + interp_mask*2B + (1B-specmask)*4B

; determine skyind, should be fast; later routines might want this.

        yobj = objpos_on_slit(slitcoords, slitobjmap, slitn, $
                              nrow=sizey, nobj=nobj)
        skyok = bytarr(sizey)
        skyok[5:sizey-6] = 1B
        ; remove rows near objects (within wid pixels)
        wid = 14
        
        if nobj eq 0 then print, 'DEIMOS_SPSLIT:   Sky only slit?'
        for k=0, nobj-1 do $
          skyok[(yobj[k]-wid/2) > 0:(yobj[k]+wid/2) < (sizey-1)] = 0B
        skyind = where(skyok, ct)
        if ct LE 6 then begin 
           print, 'DEIMOS_SPSLIT: not enough sky rows - trying narrow mask.  Slit ', slitn
           wid = 11 ; try again
           skyok[4:sizey-5] = 1B
           for k=0, nobj-1 do $
             skyok[(yobj[k]-wid/2) > 0:(yobj[k]+wid/2) < (sizey-1)] = 0B
           skyind = where(skyok, ct)
           
           if ct LE 6 then begin 
              skyind = lindgen(sizey)
              print, 'DEIMOS_SPSLIT: Not enough sky rows, object probably lost!!!   Slit ', slitn
           endif 
        endif



        spSlit = {flux: rect_spec, $
                  ivar: rect_specivar, $
;                  func: calib.func, xmin: calib.xmin, xmax: calib.xmax, $
;                  coeff: calib.coeff, $
                  mask: bitmask , $
                  slitfn: calib.slitfn, $
                  lambdax: calib.lambdax, tiltx:calib.tiltx, $
                  yobj: yobj, $
                  skyrow: skyok, $
                  skyind: skyind, $
                  slitwidth: slitwid, $
                  dlam:calib.dlam}

; look for object, while rejecting cosmic rays, using the badmask,
; using the tilts of the lines, and subtracting out the sky with a mode
 
; alternatively, we could do the following stuff at the end, after
; files are written out.  Would probably make a fair amount of sense
      
;        lightprofile = find_object(spslit, profivar=ivarout, npix=npixout, $
;                        /CR, /BPM, /MMM, /USETILT)

;        peakinfo, lightprofile, pkcol, fwhm, pk_cent=pkcent, $
;          profivar=ivarout, window=9, s2n=sig2noise

; may want to dump the lightprofile to a PS file, along with the
; predicted location of the peak (yobj, above)

        
                                ; take only the most significant peak.
                                ; The quantities we have measured are:


;        pkcol = pkcol[0]  ; location of maximum in the profile along the
                             ;slit (spatial direction), in pixels
;        fwhm = fwhm[0]    ; peak fwhm, in pixels
;        pkcent = pkcent[0] ; centroid of the peak profile
;        sig2noise = sig2noise[0] ; S/N for detecting the object, given 
                                ;an ivar-weighted integration over the
                                ;9-spatial-pixel window
                                ;defined above
;        sig2noiseperpix = sig2noise[0]/sqrt(npixout[pkcol]) ; typical signal-to-noise in an extracted spectrum, off sky lines, per (wavelength) pixel.  


; make new FITS header with selected keywords 
        keylist = ['FRAMENO', 'OUTFILE', 'EXPTIME', 'DARKTIME', 'OBSERVER', $
            'OBJECT', 'OBSTYPE', 'ROTATVAL', 'DATE-OBS', 'AIRMASS', $
            'TARGNAME', 'EQUINOX', 'DEC',  'RA', 'AZ', 'EL', 'HA', $ 
            'ST', 'UTC', 'MJD', 'MJD-ANG', 'PARANG', 'SYNOPSIS', 'DWFILNAM', $
            'SLMSKNAM', 'GRATEPOS',  'GRATING'] 

        hdr = copy_keywords(sci_header, keylist, 'spSlit')
        vers = spec2d_version() ;get version of code
        sxaddpar, hdr, 'SP2DVERS', vers, 'Version of spec2d'
        sxaddpar, hdr, 'AUTHOR', 'Berkeley DEEP team'
        sxaddpar, hdr, 'CHIPNO', chipno, 'chip number for this slitlet'
        sxaddpar, hdr, 'SLITNO', slitn, 'slit number'
        sxaddpar, hdr, 'SLITX0', calib.x0[2048], 'center of trace'
        sxaddpar, hdr, 'SLITX1', calib.x1[2048], 'center of trace'
 
        slitstr = string(slitn, format='(I3.3)')
        maskstr = size(masks, /tname) EQ 'STRING' ? $
          masks : string(masks, format='(I4.4)')
        color = (chipno gt 4) ? 'R':'B' ;red or blue side??

        inputstem = strmid(inputfile,strpos(inputfile,'.fits')-10, 10)

        fname = outdatadir+'quickslit.'+maskstr+'.'+inputstem+'.'+slitstr+color+'.fits'

        hdr2 = hdr
; Write spSlit file; use /create if exposure eq 0. 
        create = 1
        mwrfits, spSlit, fname, hdr, create=create

; NOTE: THIS DIFFERS FROM AN SPSLIT FILE IN NOT HAVING THE BSPLINE
; SAVESET ATTACHED.

      endfor
    endif  
  endfor

quickcalc, inputfile


  return
end
