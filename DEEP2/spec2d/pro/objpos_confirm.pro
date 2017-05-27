;;; a routine to compare the object positions derived from the
;;; observed spatial profile to the positions given by the design
;;; specifications of the slitmask. 
PRO objpos_confirm, n, noplot=noplot 

;;; get all of the slitfiles in directory.
  files = findfile('slit*', count=nfiles)
  midpt = FLOOR(nfiles/2.)
  bar = fltarr(8192,2) - 999.
  nrows_old = 0
  if n eq 1 then begin
    i = 0
    limit = midpt-1 
  endif else begin
    i = midpt
    limit = nfiles-1
  endelse
;;; iterate thru slit files.
  while i lt limit do begin
;;; read-in slit file.
    slit1 = mrdfits(files[i], 1, hdr1, /silent)
    slit2 = mrdfits(files[i+1], 1, hdr2, /silent)
;;; parse header for slit number.
    slitno = SXPAR(hdr1, 'SLITNO')

    if slitno EQ SXPAR(hdr2, 'SLITNO') then begin
      if i eq 0 then slitnames = slitno $
      else slitnames = [slitnames,slitno]
;;; determine number of pixels in spatial direction.
      nrows1 = n_elements(slit1.flux[0,*])
      nrows2 = n_elements(slit2.flux[0,*])
      nrows = nrows1 < nrows2
;;; get design position.
      cat_objpos = find_objpos(slitno, nrows, fwhm=cat_fwhm)      
;      print, 'Slit', slitno,  cat_objpos
;;; get position from spatial profile.
      sprof1 = find_object(slit1)
      peakinfo, sprof1, pkcol1, fwhm1, pk_cent=pkc1
      if pkc1 gt 0 then objpos1 = pkc1 else objpos1 = pkcol1
      sprof2 = find_object(slit2)
      peakinfo, sprof2, pkcol2, fwhm2, pk_cent=pkc2
      if pkc2 gt 0 then objpos2 = pkc2 else objpos2 = pkcol2
;;; build an image containing all of the 2-d slit images.
      if i eq 0 then slit_img = [slit1.flux[*,0:nrows-1], slit2.flux[*,0:nrows-1]] $
      else slit_img = [ [slit_img], [bar], [slit1.flux[*,0:nrows-1], slit2.flux[*,0:nrows-1]] ]
;;; build array to hold the object positions.
;;; for the blue slit files...
      if i eq 0 then yobj1 = objpos1 $
      else yobj1 = [yobj1, objpos1+nrows_old]
      if i eq 0 then wid1 = fwhm1 $
      else wid1 = [wid1, fwhm1]
;;; for the red slit files...
      if i eq 0 then yobj2 = objpos2 $
      else yobj2 = [yobj1, objpos2+nrows_old]
      if i eq 0 then wid2 = fwhm2 $
      else wid2 = [wid2, fwhm2]
;;; for the desing specs...
      if i eq 0 then cat_yobj = cat_objpos[0] $
      else if cat_objpos[0] gt 0. then $
        cat_yobj = [cat_yobj, cat_objpos[0]+nrows_old] $
        else cat_yobj = [cat_yobj, 0.]
      if i eq 0 then cat_wid = cat_fwhm[0] $
      else cat_wid = [cat_wid, cat_fwhm[0]]

      if i eq 0 then nrows_old = nrows else nrows_old = nrows_old + nrows + 2
      i = i+2
    endif else begin
      if i eq 0 then slitnames = slitno $
      else slitnames = [slitnames,slitno]  
;;; determine number of pixels in spatial direction.
      nrows = n_elements(slit.flux[0,*])
;;; get design position.
      cat_objpos = find_objpos(slitno, nrows, fwhm=cat_fwhm)
;      print, 'Slit ' + slitno, cat_objpos
;;; get position from spatial profile.
      sprof1 = find_object(slit)
      peakinfo, sprof1, pkcol1, fwhm1, pk_cent=pkc1
      if pkc1 gt 0 then objpos = pkc1 else objpos1 = pkcol1
;;; build an image containing all of the 2-d slit images.
      if i eq 0 then slit_img = [slit1.flux, fltarr(4096,nrows)] $
      else slit_img = [[slit_img], [bar], [slit1.flux, flatarr(4096,nrows)]]
;;; build array to hold the object positions.
;;; for the blue slit files...
      if i eq 0 then yobj1 = objpos1 $
      else yobj1 = [yobj1, objpos1+nrows_old]
      if i eq 0 then wid1 = fwhm1 $
      else wid1 = [wid1, fwhm1]
;;; for the red slit files...
      if i eq 0 then yobj2 =objpos1 $
      else yobj2 = [yobj1, objpos1+nrows_old]
      if i eq 0 then wid2 = fwhm1 $
      else wid2 = [wid2, fwhm1]
;;; for the desing specs...
      if i eq 0 then cat_yobj = cat_objpos[0] $
      else if cat_objpos[0] gt 0. then $
        cat_yobj = [cat_yobj, cat_objpos[0]+nrows_old] $
        else cat_yobj = [cat_yobj, 0.]
      if i eq 0 then cat_wid = cat_fwhm[0] $
      else cat_wid = [cat_wid, cat_fwhm[0]]

      if i eq 0 then nrows_old = nrows else nrows_old = nrows_old + nrows + 2

      i = i+1
    endelse
  endwhile
  
;;; atv the image.
  if NOT(KEYWORD_SET(noplot)) then begin
    atv, slit_img,  max=30, min=-30
;;; atv the blue objpos
    for i=0, n_elements(yobj1)-1 do $
      atvplot, [0, 4096], [1, 1]*yobj1[i]+wid1/2, color=2
    for i=0, n_elements(yobj1)-1 do $
      atvplot, [0, 4096], [1, 1]*yobj1[i]-wid1/2, color=2
;;; atv the red objpos
    for i=0, n_elements(yobj2)-1 do $
      atvplot, [4097, 8192], [1, 1]*yobj2[i]+wid2/2, color=2
    for i=0, n_elements(yobj2)-1 do $
      atvplot, [4097, 8192], [1, 1]*yobj2[i]-wid2/2, color=2
;;; atv the design objpos
    for i=0, n_elements(cat_yobj)-1 do $
      if cat_yobj[i] gt 0. then $
      atvplot, [0, 8192], [1, 1]*cat_yobj[i]+cat_wid/2, color=1, linesty=2
    for i=0, n_elements(cat_yobj)-1 do $
      if cat_yobj[i] gt 0. then $
      atvplot, [0, 8192], [1, 1]*cat_yobj[i]-cat_wid/2, color=1, linesty=2
    
    for i=0,n_elements(yobj1)-1 do atvxyouts, 0, yobj1[i], $
      string(slitnames[i], format='(I3)'), charsize=2, color=2, align=1.05
    dex = where(cat_yobj eq 0., cnt)
    if cnt gt 0 then $
      for i=0,cnt-1 do atvxyouts, 0, yobj1[dex[i]], $
      'sky', charsize=2, color=5, align=1.05

  endif

;;; make output structure (after looping thru the second time)
  img1 = slit_img
  struc1 = {objposB:yobj1, fwhmB:wid1, objposR:yobj2, fwhmR:wid2, $
            cat_objpos:cat_yobj, cat_fwhm:cat_wid}


end





