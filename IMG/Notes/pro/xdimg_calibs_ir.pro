pro nirc2_24apr05_calibs_ir, dimg
  xdimg_badpix, dimg, CCD='NIRC2W'
  xdimg_badpix, dimg, CCD='NIRC2N'
  xdimg_dark, dimg
  xdimg_dflat, dimg, /IR, CCD='NIRC2W', BPM='PixMask/BadPixMask_w.fits'
  xdimg_dflat, dimg, /IR, CCD='NIRC2N', BPM='PixMask/BadPixMask_n.fits'
  xdimg_skyflt, dimg, /IR, CCD='NIRC2W', /NOMSK
  xdimg_skyflt, dimg, /IR, CCD='NIRC2N', /NOMSK
end
