pro crremove, list, msklist, iraffil

  ;; Read list
  readcol, list, filnm, FORMAT='A'
  readcol, msklist, mskfil, FORMAT='A'

  ;; IRAFFIL
  if keyword_set(iraffil) then begin
      close, /all
      openw, 1, iraffil
  endif
  ;; Loop on list
  nfil = n_elements(filnm)

  for qq=0L,nfil-1 do begin
      ;; Loop on sub-images
      for i=1L,8 do begin
          ;; Read image
;          img = xmrdfits(strtrim(filnm[qq],2)+'.fits', i, /silent)

          ;; Run la_cosmic
          msknm = strtrim(mskfil[qq],2)+'_'+strtrim(i,2)
          la_cosmic, [strtrim(filnm[qq],2)+'.fits'], $
            masklist=[msknm+'.fits'], $
            gain=3.1, $
            readn=5.6, sigclip=4.5, sigfrac=0.5, objlim=2., niter=3, $
            blocksize=1024L, indx=i
          ;; IRAFFIL (convert to pl)
          printf, 1, 'imcopy '+msknm+'.fits '+msknm+'.pl'
          printf, 1, 'imdel '+msknm+'.fits'
      endfor
      ;; More IRAF
      printf, 1, 'crplusbpmask "', strmid(filnm[qq],3), '"'
      for i=1L,8 do $
        printf, 1, '#fixpix ', strtrim(filnm[qq],2)+'['+strtrim(i,2)+'] BPM'
  endfor

  ;;
  close, /all

  return
end
