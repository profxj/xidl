; fig_taux, psfile='Figures/fig_taux.ps', timbin=[5e8,2e8],
; OUTFIL='lx.dat'
pro sdss_dlalox, GZFIL=gzfil, DLALST=dlalst, $
              RSDSS=rsdss, RTOT=rtot, TIMBIN=timbin, STRCT=strct, $
              BINS=bins, TIFFF=tifff, OUTFIL=outfil

  if not keyword_set( GZFIL ) then gzfil = '~/SDSS/DR5_QSO/dr5_dlagz_s2n4.fits'
  if not keyword_set(BINS) then $
    bins = [ [2.2, 2.4], [2.4,2.7], [2.7, 3.0], $
             [3.0, 3.5], [3.5, 4.], [4.,5.3] ]
  sz = size(bins,/dimensions)
  nbin = sz[1]

  ;; Read gz
  gzstr = xmrdfits(gzfil, 1, /silent)
  dgzz1 = gzstr.gzz - shift(gzstr.gzz,1)

  sdss_dlastrct, alldr, /all, GZFIL=GZFIL, /stat

  ;; DLA Samples
  ndr = n_elements(alldr)

  ;; Output
  dndx = fltarr(nbin)
  medz = fltarr(nbin)
  errz = fltarr(nbin)
  errnx1 = fltarr(nbin)
  errnx2 = fltarr(nbin)
  svnd = lonarr(nbin)

  ;; Calculate dn/dz
  for ii=0L,nbin - 1 do begin
      ;; SDSS
      gd = where(alldr.zabs GE bins[0,ii] AND alldr.zabs LT bins[1,ii], nsdss)
      ;; Parse out further
      if nsdss NE 0 then begin
          msk = bytarr(nsdss)
          for jj=0L,nsdss-1 do begin
              x_radec, alldr[gd[jj]].qso_ra, alldr[gd[jj]].qso_dec, rad, decd
              idx = where(abs(gzstr.ra-rad) LT 0.001 AND $
                          abs(gzstr.dec-decd) LT 0.001, nidx)
              if nidx NE 1 then continue
              if alldr[gd[jj]].zabs GE gzstr.z1[idx] AND $
                alldr[gd[jj]].zabs LE gzstr.z2[idx] AND $
                gzstr.flg_bal[idx] NE 2 then msk[jj] = 1B
          endfor
      endif else msk = bytarr(10)
      allgd = where(msk EQ 1B, nsdss)
      if nsdss NE 0 then gd = gd[allgd] else gd = gd[0]  ; simple kludge
      medz[ii] = median([alldr[gd].zabs])

      ;; Calculate
      ndla = nsdss
      dXs = sdss_dladx( bins[0,ii], bins[1,ii], GZSTR=gzstr)
      dX = dXs
      dndx[ii] = ndla/dX
      val = x_poisscl(double(ndla), 0.8413)
      errnx1[ii] = val[0]/dX - dndx[ii]
      errnx2[ii] = dndx[ii] - val[1]/dX

      svnd[ii] = ndla
      errz[ii] = medz[ii]-bins[0,ii]
  endfor
      
  strct = { $
          medz: medz, $
          errz: errz, $
          dndx: dndx, $
          bins: bins, $
          ndla: svnd, $
          errnx1: errnx1, $
          errnx2: errnx2 $
          }
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print, 'sdss_dlalox:  All done!'

  return
end
      
      
