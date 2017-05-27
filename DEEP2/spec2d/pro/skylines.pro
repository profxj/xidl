;   slitwid    - slit width [pix]
pro sky_sanity, b, slitwid=slitwid

  spec = b.flux
  specivar = b.ivar*(b.mask EQ 0)



  if NOT keyword_set(slitwid) then slitwid = 3 ; pixels
  npix = (size(spec, /dimens))[0]
  nrow = (size(spec, /dimens))[1]
  ymid = nrow/2
;  lam0 = b.lambda0+b.lambda[*, ymid]
;  lambda = b.lambda+ b.lambda0#(fltarr(nrow)+1)
  lambda = lambda_eval(b)

  dir = getenv('DEEP_DIR')+'/'
  if dir eq '/' then message, 'you must set $DEEP_DIR'
  fname = dir+'spec2d/etc/sdss_skylines.dat'

  skyline = read_lampfile(fname)
  skyline = skyline[where(skyline.good)]


;  readcol, fname, linelam, junk, linegood, element, comments, $
;    format='(F,I,A,A,A)'

  lam0 = lambda[*, ymid]
  xstart = interpol(findgen(npix), lam0, skyline.lambda)
  w = where((xstart GT 0) and (xstart LT npix-1), nw) 
  if nw lt 1 then message, 'No lines on chip!'
  skyline = skyline[w]
  xstart = xstart[w]

  xcen1 = trace_crude(spec, specivar, xstart=xstart, ystart=ymid, $
            radius=4+slitwid/2., yset=ycen, nave=1, nmed=1, $
            maxerr=0.5, maxshifte=.9, maxshift0=2., xerr=xerr)

  lamobs = interpolate(lambda, xcen1, ycen)
  lamdiff = lamobs-transpose(skyline.lambda#(fltarr(nrow)+1))

  
  plot, lamobs, ycen, ps=6, syms=.25, xrange=minmax(lam0)+[-200, 200], /xst, /nodata, $
    title='Sky sanity check', xtitle='wavelength [Ang]', ytit='row', $
    chars=1.5

  for i=0, nw-1 do begin 
     oplot, lamobs[*, i], ycen[*, i], ps=6, syms=.25
     oplot, lamobs[*, i]+lamdiff[*, i]*250, ycen[*, i], ps=-1

     mn = djs_mean(lamdiff[*, i])
     sig = djsig(lamdiff[*, i])
     print, 'Line:', i, '  lambda:', skyline[i].lambda, '  offset:', mn, $
       ' +/-', sig, ' Ang', format='(A,I2,A,F8.2,A,F7.2,A,F6.2,A)'

  endfor 

  return
end


pro check_sanity

  b=mrdfits('spSlit.2101.063R.fits.gz',1)
  sky_sanity, b, slitwid=slitwid

  return
end



pro foo
  b=mrdfits('spSlit.2101.063R.fits.gz',1)
  sky_sanity, b, slitwid=slitwid

  
  spec = b.skymodel[*,0]
  lam = b.lambda0
  p=findpeaks(spec, mask,nsig=8)


  wp = where(p, np)
  if np lt 5 then message, 'not enough lines!'

  spec3 = [[spec], [spec], [spec]]
  wp3 = float([[wp], [wp], [wp]])

  xcen = trace_fweight(spec3, wp3, wp3-wp3+1., rad=5)
  xcen = xcen[*, 0]

  splot, spec
  soplot, xcen, xcen*0+100, ps=7

  lp = interpolate(lam, xcen)
  print, lp
  return
end


