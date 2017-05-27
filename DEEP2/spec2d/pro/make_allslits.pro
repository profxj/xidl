; top level control script to reduce a mask of DEIMOS data
; intended to be run in batch mode
; formerly known as reduce_it.pro

; you must CD to the directory you want to run this in. 

; SET NLSKY to run nonlocal sky subtraction in addition to local

; MD & DPF

pro make_allslits, planfile,chiplist=chiplist, nlsky=nlsky

  if n_elements(nlsky) eq 0 then nlsky = 0

  if n_params() lt 3 then begin
    planfile = findfile('*.plan')
    planfile = planfile[0]
  endif

  read_planfile, planfile, maskname, rawdatadir, outdatadir, flatnames, arcnames, sciencenames

 head_flat = headfits(flatnames[0])
  deimos_grating, head_flat, grating, grangle, lambda_c

  minlambda=lambda_c - 1300*(1200./grating) 
  bluelim=5500.-1500.*(minlambda lt 5500)

  slitfiles = findfile('slit*.fits', count=nfiles)
  slitfiles = slitfiles[sort(slitfiles)]
  isred = (strpos(slitfiles, 'R.fits') ne -1)
  trans = max(where(isred and lindgen(nfiles) le nfiles/2))
  if trans eq -1 then trans = floor(nfiles/2)

; Make simple 2d images, 2 per mask
  epos = strpos(slitfiles[0], '.fits')
  masknumber = strmid(slitfiles[0], 4, epos-4-4) ;get mask name '.xxxx.'
  image = slit_merge_lambda(slitfiles[0:trans], hdr,blue=bluelim)
  writefits, 'Allslits0'+masknumber+'fits', image, hdr
  if trans gt 1 then begin
    image = slit_merge_lambda(slitfiles[trans+1:nfiles-1], $
                              hdr, blue = bluelim)
    writefits, 'Allslits1'+masknumber+'fits', image, hdr
  endif



end







