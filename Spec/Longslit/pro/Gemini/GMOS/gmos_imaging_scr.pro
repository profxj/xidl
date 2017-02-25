

path = '/Users/joe/GMOS_redux/'
rawpath = path + 'Raw/'
irafpath = path + 'IRAF/'

datafiles = ['S20130905S0044.fits' $
             , 'S20130905S0046.fits' $
             , 'S20130905S0048.fits' $
             , 'S20130905S0051.fits' $
             , 'S20130905S0053.fits' $
             , 'S20130905S0055.fits']
procfiles = 'proc_imag-' + datafiles
ivarfiles = 'proc_ivar-' + datafiles

flatfiles = ['S20130903S0032.fits' $
             , 'S20130903S0033.fits' $
             , 'S20130903S0034.fits' $
             , 'S20130903S0035.fits' $
             , 'S20130903S0036.fits' $
             , 'S20130903S0037.fits' $
             , 'S20130903S0038.fits']


biasfiles = ['S20130904S0108.fits' $
             , 'S20130904S0109.fits' $
             , 'S20130904S0110.fits' $
             , 'S20130904S0111.fits' $
             , 'S20130904S0112.fits']

darkfiles =  ['S20130905S0038.fits' $
              , 'S20130905S0039.fits' $
              , 'S20130905S0040.fits' $
              , 'S20130905S0070.fits' $
              , 'S20130905S0071.fits' $
              , 'S20130905S0072.fits' $
              , 'S20130905S0073.fits' $
              , 'S20130905S0074.fits'] 

;; Create a superflat
superflatfile = path+ 'superflat-' + flatfiles[0] 
gmos_imag_superflat, rawpath + flatfiles,  superflatfile 

;;; create the superdark and the darkmask
superbiasfile = path + 'superbias-' + biasfiles[0] 
superdarkfile = path + 'superdark-' + darkfiles[0] 
gmos_superdark, rawpath + biasfiles, rawpath + darkfiles, superbiasfile, superdarkfile $
                , /IMAGING

ndata = n_elements(datafiles)
darkmask = mrdfits(superdarkfile, 1)
FOR ii = 0L, ndata-1L DO BEGIN 
;;; For a single file, run this
   long_proc, rawpath + datafiles[ii], imag1, invvar1 $
              , biasfile = superdarkfile, pixflatfile = superflatfile
   imag_arr = gmos_unmosaic(imag1*darkmask)
   ivar_arr = gmos_unmosaic(invvar1*darkmask)

   datafile_now = rawpath  + datafiles[ii]
   hdr0 = headfits(datafile_now, exten = 0)
   hdr1 = headfits(datafile_now, exten = 1)
   hdr2 = headfits(datafile_now, exten = 2)
   hdr3 = headfits(datafile_now, exten = 3)
   dims = size(imag_arr, /dim)
   nx = dims[0]
   ny = dims[1]
   biassec1 = strcompress(sxpar(hdr1, 'BIASSEC'), /rem)
   bias_arr1 = long(strsplit(biassec1, '[*:*,*:*]', /extract))
   oscan_pix1 = bias_arr1[1]-bias_arr1[0] + 1L
   
   IF ii EQ 0 THEN BEGIN
      ;; write out the darkmask in IRAF format
      dark_arr = gmos_unmosaic(darkmask)
      dark_dum1 = fltarr(nx+oscan_pix1, ny)
      dark_dum2 = dark_dum1
      dark_dum3 = dark_dum1
      ;; Chip1 and Chip2 have oscan on right side, Chip 3 has oscan on left side
      dark_dum1[0:nx-1L, *] = dark_arr[*, *, 0]
      dark_dum2[0:nx-1L, *] = dark_arr[*, *, 1]
      dark_dum3[oscan_pix1:(oscan_pix1 + nx-1L), *] = dark_arr[*, *, 2]

      darkfile_out = irafpath + 'proc_darkmask-' + darkfiles[0]
      mwrfits, 0, darkfile_out, hdr0, /create
      mwrfits, dark_dum1, darkfile_out, hdr1, /iscale
      mwrfits, dark_dum2, darkfile_out, hdr2, /iscale
      mwrfits, dark_dum3, darkfile_out, hdr3, /iscale

      ;; write out the flat field image in IRAF format
      flat = mrdfits(superflatfile, 0)
      flat_arr = gmos_unmosaic(flat)
      flat_dum1 = fltarr(nx+oscan_pix1, ny)
      flat_dum2 = flat_dum1
      flat_dum3 = flat_dum1
      ;; Chip1 and Chip2 have oscan on right side, Chip 3 has oscan on left side
      flat_dum1[0:nx-1L, *] = flat_arr[*, *, 0]
      flat_dum2[0:nx-1L, *] = flat_arr[*, *, 1]
      flat_dum3[oscan_pix1:(oscan_pix1 + nx-1L), *] = flat_arr[*, *, 2]
      
      flatfile_out = irafpath + 'proc_flat-' + flatfiles[0]
      mwrfits, 0, flatfile_out, hdr0, /create
      mwrfits, flat_dum1, flatfile_out, hdr1, /iscale
      mwrfits, flat_dum2, flatfile_out, hdr2, /iscale
      mwrfits, flat_dum3, flatfile_out, hdr3, /iscale
   ENDIF

   ;; write out the imag in IRAF format
   imag_dum1 = fltarr(nx+oscan_pix1, ny)
   imag_dum2 = imag_dum1
   imag_dum3 = imag_dum1
   ;; Chip1 and Chip2 have oscan on right side, Chip 3 has oscan on left side
   imag_dum1[0:nx-1L, *] = imag_arr[*, *, 0]
   imag_dum2[0:nx-1L, *] = imag_arr[*, *, 1]
   imag_dum3[oscan_pix1:(oscan_pix1 + nx-1L), *] = imag_arr[*, *, 2]

   procfile_out = irafpath + procfiles[ii]
   mwrfits, 0, procfile_out, hdr0, /create
   mwrfits, imag_dum1, procfile_out, hdr1, /iscale
   mwrfits, imag_dum2, procfile_out, hdr2, /iscale
   mwrfits, imag_dum3, procfile_out, hdr3, /iscale
   
   ;; Now write out the ivar in IRAF format
   ivar_dum1 = fltarr(nx+oscan_pix1, ny)
   ivar_dum2 = ivar_dum1
   ivar_dum3 = ivar_dum1
   ;; Chip1 and Chip2 have oscan on right side, Chip 3 has oscan on left side
   ivar_dum1[0:nx-1L, *] = ivar_arr[*, *, 0]
   ivar_dum2[0:nx-1L, *] = ivar_arr[*, *, 1]
   ivar_dum3[oscan_pix1:(oscan_pix1 + nx-1L), *] = ivar_arr[*, *, 2]

   ivarfile_out = irafpath + ivarfiles[ii]
   mwrfits, 0, ivarfile_out, hdr0, /create
   mwrfits, ivar_dum1, ivarfile_out, hdr1, /iscale
   mwrfits, ivar_dum2, ivarfile_out, hdr2, /iscale
   mwrfits, ivar_dum3, ivarfile_out, hdr3, /iscale

ENDFOR

   

END






