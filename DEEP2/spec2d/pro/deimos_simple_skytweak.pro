; Finkbeiner 14 July - get the sky tweak and call it dlamsky. 

;   slitwid    - slit width [pix]
pro deimos_simple_skytweak, slit, slitwid=slitwid, dlamsky=dlamsky, skyline=skyline, plot=plot

  spec     = slit.flux
  specivar = slit.ivar*(slit.mask EQ 0)

  if NOT keyword_set(slitwid) then slitwid = 3 ; pixels
  npix   = (size(spec, /dimens))[0]
  nrow   = (size(spec, /dimens))[1]
  ymid   = nrow/2
  lambda = lambda_eval(slit)

  if NOT keyword_set(skyline) then begin 
     dir = getenv('DEEP_DIR')+'/'
     if dir eq '/' then message, 'you must set $DEEP_DIR'
     fname = dir+'spec2d/etc/deep_skylines.dat'
     if findfile(fname) ne fname then fname = dir+'spec2d/etc/sdss_skylines.dat'
     skyline = read_lampfile(fname)
     skyline = skyline[where(skyline.good)]
;     stop
  endif 

  lamb = lambda[*, 0]
  lamt = lambda[*, nrow-1]
  xstartb = interpol(findgen(npix), lamb, skyline.lambda)
  xstartt = interpol(findgen(npix), lamt, skyline.lambda)
  vigpix = total(slit.mask NE 0 OR slit.ivar EQ 0, 2) gt 0.2*nrow
  wheregood = where(vigpix eq 0)
  goodrange = minmax(wheregood)

  w = where( ((xstartb < xstartt) GT (goodrange[0] > 5)) and ((xstartb > xstartt) LT (goodrange[1] < (npix-6))), nw) 

  if nw lt 1 then message, 'No lines on chip!'
  skylinew = skyline[w]
  xstart = (xstartb[w]+xstartt[w])/2.

; trace sky lines
  xcen1 = trace_crude(spec, specivar, xstart=xstart, ystart=ymid, $
            radius=fix(2+slitwid/2.), yset=ycen, nave=1, nmed=1, $
            maxerr=0.5, maxshifte=.9, maxshift0=2., xerr=xerr)

  lamobs  = interpolate(lambda, xcen1, ycen)
  dlamsky = lamobs-transpose(skylinew.lambda#(fltarr(nrow)+1))

  if keyword_set(plot) then begin 
     plot, lamobs, ycen, ps=6, syms=.25, xrange=minmax(lamb)+[-200, 200], $
       /xst, /nodata, title='Sky sanity check', $
       xtitle='wavelength [Ang]', ytit='row', chars=1.5
     
     for i=0, nw-1 do begin 
        oplot, lamobs[*, i], ycen[*, i], ps=6, syms=.25
        oplot, lamobs[*, i]+dlamsky[*, i]*250, ycen[*, i], ps=-1
        
        mn = djs_mean(dlamsky[*, i])
        sig = djsig(dlamsky[*, i])
        print, 'Line:', i, '  lambda:', skylinew[i].lambda, '  offset:', mn, $
          ' +/-', sig, ' Ang', format='(A,I2,A,F8.2,A,F7.2,A,F6.2,A)'
        
     endfor 
  endif 
  
  return
end

