; DEIMOS ID test 2001-Aug-22
; DEIMOS readiness review plots 2001-Sep-18
; D. Finkbeiner
; try 2d- bsplines! MD, 27sep

; plate scales
; 1- 8.56  5- 8.56
; 2- 8.47  6- 8.47
; 3- 8.47  7- 8.47
; 4- 8.56  8- 8.56
;pro readiness, chipno, ext1, ext2,  flat=flat

  device, pseudo=8

  ybin = 8
  chipno = 2
  flat = 1

; read subimage
  fname = '/deepscr2/marc/DEIMOS/sci1_rot45fast_g3_2001.fits'   ; flat
  arcname = '/deepscr2/marc/DEIMOS/sci1_rot45fast_g3_2000.fits'
  im = deimos_read(fname, chipno) -1600 ;remove overscan, roughly
  
  invvar = (im LT 30000)
;  im = im*invvar + 3000.*(1-invvar)

; Find some traces - ncoeff=4 means third order fit 
  deimos_trace_crude, im, tset1, tset2, imgbin, ybin=ybin, ncoeff=4, $
    invvar=invvar

; weed out the black sheep
  clean_tset, tset1
  clean_tset, tset2
  match_tset, tset1, tset2

  delvarx, ypos1
  delvarx, ypos2
  traceset2xy, tset1, ypos1, xpos1
  traceset2xy, tset2, ypos2, xpos2

  ny = 4096
  yind = lindgen(ny/ybin)*ybin
;  atv, imgbin, min=0, max=60000
;  atvplot, xpos1[yind, *], ypos1[yind, *]/ybin, ps=3,color=4
;  atvplot, xpos2[yind, *], ypos2[yind, *]/ybin, ps=3,color=6
  
; ------------------------------------------------------------------------
; Plot 1 - Flats
;  dfpsplot, 'flat1.ps', /sq, /color, bits=8
;  imps = imgbin[0:511, *]
;  display,bytscl(imps,min=-100,max=65000), xtit='spatial [pix]', ytit=$
;    'lambda [pix]', chars=1.5, xmargin=[7, 2]
;  nline = (size(xpos1))[2]
;  for i=0, nline-1 do oplot, xpos1[yind, i], ypos1[yind, i]/ybin,color=4
;  for i=0, nline-1 do oplot, xpos2[yind, i], ypos2[yind, i]/ybin,color=6
;  dfpsclose

; read in mask file
  maskfile = '/deep0/marc/deep/targetselection/testmask.11.fits'
  m = mrdfits(maskfile, 2)

; slit tops:
  xfound = xpos1[2048, *]
  deimos_slit_id, m, xfound, indtop, xtop, /top, xoffs = xoffs, pltscl=pltscl
;  atvplot, xtop, fltarr(n_elements(xtop))+2048/ybin, ps=7, thick=2, syms=1.5
  info = string('Chipno:', chipno, ' Top Scale:', pltscl, $
                ' (pix/asec)  Xoffs:', xoffs, format='(A,I2,A,F8.4,A,F10.3)')
  print, info

;stop
; ------------------------------------------------------------------------
; Plot 2 - slit ID
;  dfpsplot, 'slitid.ps', /sq
;  plot, xfound-xtop, ps=1, xtit='Slit number', ytit='Position resid. [pix]', $
;    chars=1.2, title=info
;  dfpsclose

; slit bottoms
  xfound = xpos2[2048, *]
  deimos_slit_id, m, xfound, indbot, xbot, /bot, xoffs = xoffs, pltscl=pltscl
;  atvplot, xbot, fltarr(n_elements(xbot))+2048/ybin, ps=6, thick=2, syms=1.5
  print, 'Chipno:', chipno, ' Bot Scale:', pltscl, ' (pix/asec)  Xoffs:', xoffs, $
    format='(A,I2,A,F8.4,A,F10.3)'


  arc = deimos_read(arcname, chipno) -1600
  arcivar = (arc*0+1./900*(arc LT 60000))
;  atv, arc
;  atvplot, xpos1, ypos1, ps=3, color=4
;  atvplot, xpos2, ypos2, ps=3, color=6


; ------------------------------------------------------------------------
; Plot 3 - Arcs
;  dfpsplot, 'arc1.ps', /sq, /color, bits=8
;  imps = arc[0:511, 0:511]
;  display,bytscl(imps,min=-10000,max=65000), xtit='spatial [pix]', ytit=$
;    'lambda [pix]', chars=1.5, xmargin=[7, 2]
;  nline = (size(xpos1))[2]
;  for i=0, nline-1 do oplot, xpos1[*, i], ypos1[*, i],color=4
;  for i=0, nline-1 do oplot, xpos2[*, i], ypos2[*, i],color=6
;  dfpsclose



; Read lamp file into a structure
  lampfilename = '~dfink/deep/spec2d/etc/lamphgcdne.dat'
  lamps = read_lampfile(lampfilename)

  lampwave = lamps.lambda
  intensity = lamps.intensity


  for k=0, 27 do begin ;loop over all slitlets
;  k = 3

	  x0 = xpos1[*, k]
 	 x1 = xpos2[*, k]
; extract arc
  	extvert = deimos_rectify_slit(arc, arcivar, x0, x1, /interp, $
                                /recen, xshift=xshift, npad=0)

; extract corresponding flat 

  	if keyword_set(flat)  then begin
     	extflatvert = deimos_rectify_slit(im, invvar, x0, x1, /interp, $
                                       /recen, xshift=xshift, npad=0)
     	extvert = deimos_flatfield(extvert,  extflatvert)
  	endif

  	ext = transpose(extvert)
  
; just as a test
  	ymid = (size(ext, /dim))[1]/2
  	spec = ext[*, ymid]
  
  	spec = reform(spec, n_elements(spec))
  
;  lamscale = 0.4314 ; Angstrom per pixel
  
  	t1 = systime(1)
  
  	wset1 = lin_arcfit_guess(spec, lampwave, intensity, color=color, $
                           func=func, bestcorr=bestcorr, acoeff=acoeff,  $
                           dcoeff=dcoeff, nsteps=nsteps)
 	print, 'Elapsed time', systime(1)-t1
  	print, 'shift in pix ', wset1.coeff[0]
  	print, 'masky: ', m[indtop[k]].masky
  	traceset2xy, wset1, xtemp, lambda
;  splot, lambda, spec
;  soplot, lampwave, intensity*4, ps=7
  
  	print, 'Wavelength range on chip in this trace: ', lambda[0], $
    	  'A to ',  lambda[ny-1], 'A'
  	print, 'pixel size (bottom) ', lambda[1]-lambda[0], 'A'
 	 print, '              (top) ', lambda[ny-1]-lambda[ny-2], 'A'
  	print
  	print
  
;stop  
; ------------------------------------------------------------------------
; Plot 4 - wavematch
  
;  dfpsplot, 'wavematch.ps', /sq
;  plot, lambda, spec, thick=1, xtit='lamgda [Ang]', ytit='counts', $
;    title='Wavelength matching', chars=1.5, /xst
;  oplot, lampwave, intensity*4, ps=7, syms=1.5, thick=2
;  dfpsclose
  
  	arcline_x = traceset2pix(wset1, lampwave)
  
  	deimos_arcfit, ext, arcline_x, lamps, wset, lamdif=lamdif
; lamdif is in Ang
  	lno = findgen(n_elements(lamdif))/(size(lamdif))[2]
;  splot, lno, transpose(lamdif)


;stop  
; ------------------------------------------------------------------------
; Plot 5 - arcres

;  dfpsplot, 'arcres.ps', /square
;  plot, lno, transpose(lamdif), ps=1, xtit='Line number', $
;    ytitle='Arc fit residual [Ang]', title='Wavelength residual', $
;    yrange=[-0.04, 0.04], chars=1.5
;  dfpsclose
  
  	traceset2xy, wset, pixnorm, wave

  	nrow = (size(ext, /dimens))[1]

  	r0 = 10
  	r1 = nrow-11
        if(r1 lt r0) then r1=r0  ;short rows -- just a test
  	arcwave = wave[*, r0:r1]
  	arcflux = ext[*, r0:r1]
  	arcivar = (ext*0+1./900*(ext LT 60000))[*, r0:r1]
 
     ystart=(xpos1[0,k]+xpos2[0,k])/2.
     ypos=ystart-(r1-r0)/2 +findgen(r1-r0+1)
     if(k eq 0) then begin
        all_arcwave=transpose(arcwave)
        all_arcflux=transpose(arcflux)
        all_arcivar=transpose(arcivar)
        all_ypos=ypos
     endif else begin
        all_arcwave=[all_arcwave, transpose(arcwave)] ;concatenate arrays
        all_arcflux=[all_arcflux, transpose(arcflux)]
        all_arcivar=[all_arcivar, transpose(arcivar)]
        all_ypos   =[all_ypos, ypos]

     endelse 
endfor
       all_arcwave=transpose(all_arcwave)  ;orient arrays so lambda is first
       all_arcflux=transpose(all_arcflux)
       all_arcivar=transpose(all_arcivar) 
       ncols=(size(all_arcwave,/dimen))[0]
       yarray=replicate(1,ncols)#all_ypos


delvarx, im,  arc, pixnorm,  imgbin 
save,  filename='readiness.2.sav'

; need to generate an array of same size as others, listing x2 positions
 ;ind = sort(arcwave)
; try a model of a high resolution spline in spectral direction, low
; order polynomial function in spatial direction


 minlambda = min(all_arcwave, max=maxlambda)
 bkpt = findgen(4096)*(maxlambda-minlambda)/4095. +minlambda 
;
;set breakpoints to be uniform in this space, but why are we in log
;lambda??
;
 tstart = systime(1)
 sset = bspline_iterfit(all_arcwave, all_arcflux, invvar=all_arcivar, $
         x2=yarray, npoly=4, xmin=min(all_ypos), xmax=max(all_ypos), $         
         upper=20, lower=20, maxiter=3, bkpt=bkpt, yfit=skyfit)
 print, 'time for bspline estimation:', systime(1)-tstart 
 ; wsky = (findgen(100000)+0.5)/100000*(max(wave)-min(wave))+min(wave)
 ; sky = bspline_valu(wsky, sset, x2=yarray)

;  splot, 10^arcwave, arcflux, ps=3
;  soplot,10^wsky,sky, color=2

;  dfpsplot, 'arcfit.ps', /sq
;  plot, 10^arcwave, arcflux, ps=3, xr=[8120, 8170], yr=[1000, 2500], $
;    xtit='wavelength [Ang]', ytit='counts', title='Arcfit', chars=1.5
;  oplot,10^wsky,sky, color=3, thick=3
;  xyouts, 8142, 1400, 'Ghost', /data, chars=1.5
;  dfpsclose

   atv,  all_arcflux - skyfit,  max=1000.,  min=-1000.

  ;arcmodel = bspline_valu(wave, sset)

  ;mask = ext LT 60000
  ;block = [[(ext-arcmodel)*mask],[arcmodel-1800],[ext-1800]]
  ;img = block[1010:1521, *]

; ------------------------------------------------------------------------
; Plot 6 - Arc subtraction
;  dfpsplot, 'arcsub.ps', /sq, /color, bits=8
;  bim = bytscl(img,min=-1000,max=5000, top=245)+8B
;  bim2 = bytscl((ext-arcmodel)[1010:1521, *],min=-250,max=250, top=245)+8B
;  display,[[bim2], [bim]], title='arc; arcmodel; residual; res stretched', $
;    xtit='lambda [pix]', ytit=$
;    '[pix]', chars=1.5, xmargin=[7, 2]
;  dfpsclose


  ;atv,block, max=5000, min=-1000



end
;  7293  -31.6 arcsec
;  7489   46.0
;  7224  -59
;  7408,  12.3
