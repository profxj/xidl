pro hires_qsubbias, fil, outfil, IMG=ovimg, GAIN=gain, HEAD=main_head, $
                   DEBUG=debug

  ;;
  if not keyword_set( GAIN ) then gain = 1.14
  if not keyword_set( OUTFIL ) then begin
      pos = strpos(fil, '/')
      pos2 = strpos(fil, 'fits')
      pos = pos > 0
      len = strlen(strtrim(fil,2))
      outfil = 'OV/'+strmid(fil,pos+1,pos2-pos+7)
  endif

  if x_chkfil(outfil+'*',/silent) NE 0 $
    and not keyword_set( CLOBBER ) then begin
      print, 'hires_subbias: File exists, returning'
      if arg_present(OVIMG) then $
        ovimg = xmrdfits(outfil, 0,main_head,/silent)
      return
  endif

  ;; Grab blue frame
  raw = xmrdfits(fil, 1, head, /fscale, /silent)
  sz = size(raw, /dimensions)

  ;; OV subtract
  ov = djs_median(raw[683:*,*],1)
  fitstr = x_setfitstrct()
  f1 = x_fitrej(findgen(sz[1]), ov,fitstr=fitstr)
  ovimg = raw[0:681,*] - replicate(1.,682L)#f1

  if keyword_set(DEBUG) then stop
  ovimg = ovimg*gain

  ;;
  mkhdr, main_head, ovimg
  sxdelpar, main_head, 'END'
  old_head = xheadfits(fil, /silent)
  sxdelpar, old_head, 'BITPIX'
  sxdelpar, old_head, 'SIMPLE'
  sxdelpar, old_head, 'NAXIS'
  sxdelpar, old_head, 'BZERO'
  sxdelpar, old_head, 'DATSUM'
;  sxdelpar, old_head, 'END'

  nhd = n_elements(old_head)

  for qq=0L,nhd-1 do begin
      if strlen(strtrim(old_head[qq],2)) GT 0 then $
        main_head = [main_head, old_head[qq]]
  endfor
  
;  mjd = sxpar(main_head, 'MJD')
;  exptim = sxpar(main_head, 'EXPTIME')
;  sz = size(ovimg, /dimension)
;  sxaddpar, main_head, 'NAXIS', 2
;  sxaddpar, main_head, 'NAXIS1', sz[0]
;  sxaddpar, main_head, 'NAXIS2', sz[1]
;  sxaddpar, main_head, 'BITPIX', -32
;  sxdelpar, main_head, 'BZERO'
;  sxdelpar, main_head, 'DATASUM'
  mwrfits, ovimg, outfil, main_head, /create
  spawn, 'gzip -f '+outfil

  return

end
