; top level control script to reduce a mask of DEIMOS data
; intended to be run in batch mode
; formerly known as reduce_it.pro

; you must CD to the directory you want to run this in. 

; SET NLSKY to run nonlocal sky subtraction in addition to local

; MD & DPF

pro make_2d_file, planfile,filestring

  if NOT keyword_set(planfile) then message, 'You must specify a plan file!'

  read_planfile, planfile, maskname, rawdatadir, outdatadir, flatnames, arcnames, sciencenames

  head_flat = headfits(flatnames[0])
  deimos_grating, head_flat, grating, grangle, lambda_c

  minlambda=lambda_c - 1300*(1200./grating) 
  bluelim=5500.-1500.*(minlambda lt 5500)

  slitfiles = findfile(filestring, count=nfiles)

  slitfiles = slitfiles[sort(slitfiles)]
  isred = (strpos(slitfiles, 'R.fits') ne -1)

; Make simple 2d images, 2 per mask
  epos = strpos(slitfiles[0], '.fits')
  masknumber = strmid(slitfiles[0], 4, epos-4-4) ;get mask name '.xxxx.'
  image = slit_merge_lambda(slitfiles[0:nfiles-1], hdr, blue=bluelim)
  ;image = slit_merge_lambda(slitfiles[0:nfiles-1], hdr,blue=bluelim)
  writefits, 'Highz'+masknumber+'fits', image, hdr

end







