pro mage_makepix, arcfile, pixfile, tset_slits, chk=chk

   mage_proc, arcfile, arcimg
   pixset = long_wavepix(arcimg, tset_slits, fwhm=3.0 $
                         , box_radius = 3.0 $
                         , sig_thresh = 3.0 $
                         , pkwdth = 5.0 $
                         , TOLER = 2.0 $
                         , CHK = CHK); removing this for
                                ;now...doesn't do anything

   piximg = long_wpix2image(pixset, tset_slits)
   mwrfits, piximg, pixfile, /create
   mwrfits, pixset, pixfile

end
