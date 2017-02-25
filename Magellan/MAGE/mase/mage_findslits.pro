function mage_findslits, tracefiles, orderfile = orderfile, no_chk_trc=no_chk_trc

  nfiles = n_elements(tracefiles)
  mage_proc, tracefiles[0], trcimg, trcivar, hdr = hdr
  trcimg = trcimg*(trcivar GT 0) 
  IF nfiles GT 1 THEN BEGIN
     dims = size(trcimg, /dim)
     imgarr = make_array(dimension = [dims, nfiles], /float)
     imgarr[*, *, 0] = trcimg
     FOR ifile = 1L, nfiles-1L DO BEGIN
        mage_proc, tracefiles[ifile], tmp, tmpivar, hdr = hdr
        imgarr[*, *, ifile] = tmp*(tmpivar GT 0)
     ENDFOR
     trcimg = djs_avsigclip(imgarr, 3, sigrej = sigrej $
                            , maxiter = maxiter, inmask = inmask) 
  ENDIF
  chk_trc = keyword_set(no_chk_trc) ? 0:1 ; use old IDL convention
  tset_slits = mage_traceorders(trcimg, chk=chk_trc)
  ordermask = mage_ordermask(tset_slits) 
  ;;slitmask = long_slits2mask(tset_slits)
  ;;ordermask = 0 * slitmask
  ;;ordermask[WHERE(slitmask GT 0)] = -slitmask[WHERE(slitmask GT 0)] + 21L
   if keyword_set(orderfile) then begin
      mwrfits, ordermask, orderfile, /create
      mwrfits, tset_slits, orderfile
   endif else begin
      file = "Orders.fits"
      mwrfits, ordermask, file, /create
      mwrfits, tset_slits, file
   endelse

   return, tset_slits

   print, "mage_findslits: all done!"

end
