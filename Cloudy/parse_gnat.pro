pro parse_gnat

  ;; file
  fil = 'gnat_CIE.dat'

  ;; 
  close, /all
  openr, 1, fil

  ;; 
  nlin = 327L-127+1

  lin = ''

  ;; Dummy lines
  dumlin = 126
  for ss=0L,dumlin-1 do readf, 1, lin

  tmp = { $
        T: 0., $
        H: fltarr(2), $
        He: fltarr(3), $
        C: fltarr(7), $
        N: fltarr(8), $
        O: fltarr(9), $
        Neon: fltarr(11), $
        Mg: fltarr(13), $
        Si: fltarr(15), $
        S: fltarr(17),$
        Fe: fltarr(27) $
        }
  gnat_strct = replicate(tmp, nlin)

  ;; Good stuff
  for ss=0L,nlin-1 do begin
      readf, 1, lin
      prs = strsplit(lin, ' ', /extract)
      ind = prs[1:*]
      ;; Parse
      gnat_strct[ss].T  = prs[0]
      gnat_strct[ss].H  = ind[0:1]
      gnat_strct[ss].He = ind[2:4]
      gnat_strct[ss].C  = ind[5:11]
      gnat_strct[ss].N  = ind[12:19]
      gnat_strct[ss].O  = ind[20:28]
      gnat_strct[ss].Neon = ind[29:39]
      gnat_strct[ss].Mg = ind[40:52]
      gnat_strct[ss].Si = ind[53:67]
      gnat_strct[ss].S  = ind[68:84]
      gnat_strct[ss].Fe = ind[85:111]
  endfor
  close, /all
  
  outfil = 'gnat_CIE.fits'
  mwrfits, gnat_strct, outfil, /create
  spawn, 'gzip -f '+outfil


  return
end
