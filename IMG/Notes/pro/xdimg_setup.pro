pro lris_13jan99_setup, dimg

; Create the structure
  xdimg_strct, dimg, 'LRISR', 'Keck', /mkdir
; Edit
  dimg[15:16].flg_anly = 0
  dimg[15:16].type = 'TWI'
  dimg[85].flg_anly = 0
  dimg[86].flg_anly = 0
  dimg[97].flg_anly = 0
  dimg[99].flg_anly = 0
  dimg[111].flg_anly = 0
  dimg[128].flg_anly = 0
; Write
  write_dimgstr, dimg, FITS='lris_13jan99.fits'
end
