FUNCTION tspec_findslits, traceflat, orderfile, ksize = ksize

  tspec_proc, traceflat, trcimg, hdr = hdr
  tset_slits = tspec_traceorders(trcimg, ksize = ksize, /CHK, /CRUDE)
  ordermask = tspec_ordermask(tset_slits) 
  
  mwrfits, ordermask, orderfile, /create
  mwrfits, tset_slits, orderfile
  
END


