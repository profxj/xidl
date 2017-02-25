pro z1p_to_fits

  ;; Read CSV
  ;; Need to unzip the file first
  readcol, 'z1p5_4.csv', QSO, num, wave, trans, zabs, zsig, bval, bsig, logN, sig_logN, $
           ref, comment, FORMAT='A,L,D,A,D,D,F,F,F,F,A,A', delimi=','

  nline = n_elements(num)

  ;; Create binary FITS table
  line = { $
         qso: '', $
         wave: 0.d, $
         trans: '', $
         zabs: 0.d, $
         zsig: 0.d, $
         bval: 0., $
         bsig: 0., $
         logN: 0., $
         sigN: 0., $
         ref: '', $
         comment: '' $
         }
  all_lines = replicate(line, nline)
  all_lines.qso = qso
  all_lines.wave = wave
  all_lines.trans = strtrim(trans,2)
  all_lines.zabs = zabs
  all_lines.zsig = zsig
  all_lines.bval = bval
  all_lines.bsig = bsig
  all_lines.logN = logN
  all_lines.sigN = sig_logN
  all_lines.ref = strtrim(ref,2)
  all_lines.comment = comment

  outfil = 'z1p5_4.fits'
  mwrfits, all_lines, outfil, /create
  spawn, 'gzip -f '+outfil
  return
end
