; DEIMOS ID test 2001-Aug-22
; DEIMOS readiness review plots 2001-Sep-18
; D. Finkbeiner
; try 2d- bsplines! MD, 27sep
; use new data from end of 2001


; plate scales
; 1- 8.56  5- 8.56
; 2- 8.47  6- 8.47
; 3- 8.47  7- 8.47
; 4- 8.56  8- 8.56
;pro 2dreduce, chipno, ext1, ext2,  flat=flat

  device, pseudo=8

  ybin = 8
  chipno = 1
  flat = 1
  anamorph = 1.6 ; anamorphic factor, which is a function of grating and tilt
  spline =  1 ; do spline fit

; read subimage

  maskno = 3110

  deimos_data = getenv('DEIMOS_DATA')+'/'
  if deimos_data eq '/' then message, 'You need to set $DEIMOS_DATA!'

  deep_dir = getenv('DEEP_DIR')+'/'
  if deep_dir eq '/' then message, 'You need to set $DEEP_DIR!'

  maskstr = string(maskno, format='(I4.4)')
  maskdir = deimos_data+maskstr+'/'

  fname = maskdir+'davis.1000.fits' ; flat
  arcname = maskdir+'davis.1001.fits'

; perhaps things like this should be versioned???
  maskname = '/deep0/marc/deep/spectral/deimos_badmask.fits.Z'

  im = deimos_read(fname, chipno) -1600 ;remove overscan, roughly
  deimos_badchip = mrdfits(maskname, chipno) ;get mask file of bad regions
  
  invvar = (im LT 60000)*deimos_badchip

; what are these???  ptr should be reserved for IDL pointers
  cum_flat1d = fltarr(2048) ; resulting 1d flat showing slit variations
  cum_slitlet_ptr = 0
  slitlet_ptr = intarr(100) -1  ;array of initial points of each slitlet

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
  nslitlets = (size(xpos1, /dimens))[1]

  ny = 4096
  yind = lindgen(ny/ybin)*ybin
  atv, imgbin, min=0, max=60000
  atvplot, xpos1[yind, *], ypos1[yind, *]/ybin, ps=3,color=4
  atvplot, xpos2[yind, *], ypos2[yind, *]/ybin, ps=3,color=6
  
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

; read in mask file-- hardwired for the moment.  Note that the FITS
;                     table information will in future be included in
;                     data file. This routine does not yet handle
;                     multiple objects/slitlet.
  maskfile = '/deep0/marc/deep/targetselection/masks/add_xdistortion/' + $
    'testmask.3110.fits'
  mask = mrdfits(maskfile, 2)
  maskfitsfile = maskdir+'masl3110.fits'
  objectcat = deimos_tables(maskfitsfile, bluslits=slitcoords) 
; 
; end of hardwired region
;


; slit tops:
; note that deimos_slit_id should actually refer to objectcat and
; bluslits, not the testmask.xxxx.fits file used to construct the
; mask! TBD !!

  xfound = xpos1[2048, *]

  deimos_slit_id, mask, xfound, indtop, xtop, /top, xoffs = xoffs, $
      pltscl=pltscl
  atvplot, xtop, fltarr(n_elements(xtop))+2048/ybin, ps=7, thick=2, syms=1.5
  info = string('Chipno:', chipno, ' Top Scale:', pltscl, $
                ' (pix/asec)  Xoffs:', xoffs, format='(A,I2,A,F8.4,A,F10.3)')
  print, info

;TBD all this needs to be cleaned up!
  objno_selected = mask[indtop].objno 
  object_mask = long( objectcat.object)
;unique object numbers of selected objects in mask

  nobj = n_elements(objno_selected)
  blu_ptr = intarr(nobj)  ;list of objects to examine on this chip
  for i=0, nobj-1 do begin
    ii = where(object_mask eq objno_selected[i], jjj)
    if jjj eq 0 then print,  'problem in matching lists!'
    blu_ptr[i] = ii[0] ;get first in list
  endfor
  slitcoords = slitcoords[blu_ptr] ;keep only appropriate slitlets


; ------------------------------------------------------------------------
; Plot 2 - slit ID
  dfpsplot, 'slitid.ps', /sq
 plot, xfound-xtop, ps=1, xtit='Slit number', ytit='Position resid. [pix]', $
   chars=1.2, title=info
  dfpsclose

; slit bottoms
  xfound = xpos2[2048, *]
  deimos_slit_id, mask, xfound, indbot, xbot, /bot, xoffs = xoffs, $
     pltscl=pltscl
  atvplot, xbot, fltarr(n_elements(xbot))+2048/ybin, ps=6, thick=2, syms=1.5
  print, 'Chipno:', chipno, ' Bot Scale:', pltscl, ' (pix/asec)  Xoffs:', xoffs, $
    format='(A,I2,A,F8.4,A,F10.3)'


  arc = deimos_read(arcname, chipno, header=arc_header) -1000
;  arc_ivar = (arc*0+1./900*(arc LT 60000))*deimos_badchip
  arc_ivar = float(arc*0.+1./(abs(arc)+1))*deimos_badchip 
;give saturated points nonzero wt.
  atv, arc
  atvplot, xpos1, ypos1, ps=3, color=4
  atvplot, xpos2, ypos2, ps=3, color=6


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
;  lampfilename = '/deep0/marc/deep/spectral/lamphgcdne.dat'
  lampfilename = '/deep0/marc/deep/spectral/lamplist.dat'
  lamps = read_lampfile(lampfilename)


  lampwave = lamps.lambda
  intensity = lamps.intensity

;
;TBD -- replace all this with optical model calls
;


;get optical model estimates for all slitlets in question
  model_lambda =  deimos_omodel(chipno, slitcoords, arc_header)
  xtemp = findgen(4096)


  for slitno=0, nslitlets-1 do begin ;loop over all slitlets
     tstart = systime(1) 
     
     
     x0 = xpos1[*, slitno] 
     x1 = xpos2[*, slitno] 
; extract arc
     extarcvert = deimos_rectify_slit(arc, arc_ivar, x0, x1, /interp, $
                                      xshift=xshift, npad=0) ;don't /recen
     extivar = deimos_rectify_slit(arc_ivar, arc_ivar, x0, x1, /interp, $
                                   xshift=xshift, npad=0)
     
     sizex = (size(extarcvert, /dimens))[0]
     
     slitlet_ptr[slitno] = cum_slitlet_ptr
     
; extract corresponding flat 
     
     if keyword_set(flat)  then begin
        extflatvert = deimos_rectify_slit(im, invvar, x0, x1, /interp, $
                                          xshift=xshift, npad=0)
        extarcvert = deimos_flatfield(extarcvert,  extflatvert, $
                                      flat1d=outflat, invvar=extivar,  /twod)
        cum_flat1d[cum_slitlet_ptr:cum_slitlet_ptr + sizex-1] = outflat
        
     endif
     cum_slitlet_ptr = cum_slitlet_ptr + sizex ;increment for next slitlet
     
     
     ext = transpose(extarcvert)
     extivar = transpose(extivar)
     
     
     ext = convol(ext,[1.,2.,1.],3.,/center) 
;slightly smooth arc in spectral direction 
     
     ymid = (size(ext, /dim))[1]/2
     spec = ext[*, ymid]
     
     spec = reform(spec, n_elements(spec))
     
;since we have an optical model, use it instead of guessing initial solution!
     if model_lambda[slitno].lambda_y[0] eq 0 then begin 
; use this only in case of no optical model!
        print, 'no optical model available for slitlet ', slitno
        
        t1 = systime(1)
        
        wset1 = lin_arcfit_guess(spec, lampwave, intensity, color=color, $
                                 func=func, bestcorr=bestcorr, acoeff=acoeff,  $
                                 dcoeff=dcoeff, nsteps=nsteps)
        print, 'Elapsed time for guess ', systime(1)-t1
;  	print, 'shift in pix ', wset1.coeff[0]
;  	print, 'masky: ', mask[indtop[slitno]].masky
        traceset2xy, wset1, xtemp, lambda
        splot, lambda, spec
        soplot, lampwave, intensity*4, ps=7
        arcline_x = traceset2pix(wset1, lampwave)
;  	print, 'Wavelength range on chip in this trace: ', lambda[0], $
;    	  'A to ',  lambda[ny-1], 'A'
;  	print, 'pixel size (bottom) ', lambda[1]-lambda[0], 'A'
; 	 print, '              (top) ', lambda[ny-1]-lambda[ny-2], 'A'
;  	print
;  	print
        
; ------------------------------------------------------------------------
; Plot 4 - wavematch
        
;  dfpsplot, 'wavematch.ps', /sq
;  plot, lambda, spec, thick=1, xtit='lamgda [Ang]', ytit='counts', $
;    title='Wavelength matching', chars=1.5, /xst
;  oplot, lampwave, intensity*4, ps=7, syms=1.5, thick=2
;  dfpsclose
        
        
; use model to get initial lambda
     endif else begin
        lambda =  poly(xtemp, model_lambda[slitno].lambda_y ) ;lambda(y)
        
        splot, lambda, spec, thick=1, xtit='lamgda [Ang]', ytit='counts', $
          title='Wavelength matching', chars=1.5, /xst
        soplot, lampwave, intensity*4, ps=7, syms=1.5, thick=2
        
        arcline_x = poly(lampwave,  model_lambda[slitno].y_lambda )
        
     endelse
     
     slitpa_in = (atan(tan(mask[indtop[slitno]].slitpa/!radeg)/anamorph) - $
                  atan(2.*mask[indtop[slitno]].maskx/480.*(40./480.)))*!radeg 
;tilt of slit, including effect of bowing, assuming bow of 40 pixels
;over range of 4096 pixels in quadratic behavior. Anamorphic factor
; must be included!
     wave =  deimos_arcfit(ext, arcline_x, lamps, wset, lamdif=lamdif,  $
                           ncoeff=4, slitpa=slitpa_in,  arcivar=extivar,  spline=spline)
     
; lamdif is in Ang
     if spline eq 0 then begin
        lno = findgen(n_elements(lamdif))/(size(lamdif))[2]
        splot, lno, transpose(lamdif), xtit='line number', ytit='delta lambda'
; ------------------------------------------------------------------------
; Plot 5 - arcres
        
        dfpsplot, 'arcres.ps', /square
        plot, lno, transpose(lamdif), ps=1, xtit='Line number', $
          ytitle='Arc fit residual [Ang]', title='Wavelength residual', $
          yrange=[-0.04, 0.04], chars=1.5
        dfpsclose
     endif
     
     
     
     ncol = (size(ext, /dimens))[0]
     nrow = (size(ext, /dimens))[1]
     bitmask = bytarr(ncol, nrow)
     
; Generate spSlit structure for output
     spSlit = {flux: ext, $
               ivar: extivar, $
               lambda: wave, $
               mask: bitmask}
     


     slitstr = string(mask[indtop[slitno]].slitn, format='(I3.3)')
     fname = maskdir+'spSlit.'+maskstr+'.'+slitstr+'.fits'
     print, ' Writing: ', fname
     mwrfits, spSlit, fname
     
     r0 = 0
     r1 = nrow-1                ;keep everything for saving
     
     if (r1 lt r0) then r1=r0    ;short rows -- just a test
     arcwave = wave[*, r0:r1]
     arcflux = ext[*, r0:r1]
     arcivar = extivar[*, r0:r1]
     
     ystart=(xpos1[0,slitno]+xpos2[0,slitno])/2.
     ypos=ystart-(r1-r0)/2 +findgen(r1-r0+1)
     if(slitno eq 0) then begin
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

     print, systime(1)-tstart,  ' seconds to process slitlet ', slitno

  endfor

  all_arcwave=transpose(all_arcwave) ;orient arrays so lambda is first
  all_arcflux=transpose(all_arcflux)
  all_arcivar=transpose(all_arcivar) 
  ncols=(size(all_arcwave,/dimen))[0]
  yarray=replicate(1,ncols)#all_ypos
  
  
  delvarx, im,  arc,  imgbin,  deimos_badchip, invvar, arc_ivar
  save,  all_arcwave, all_arcflux,  all_arcivar,  cum_flat1d, yarray, ncols, $
    slitlet_ptr, filename='2dtest.sav'
  
; need to generate an array of same size as others, listing x2 positions
                                ;ind = sort(arcwave)
; try a model of a high resolution spline in spectral direction, low
; order polynomial function in spatial direction
  
  
end

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
 print, 'time for bspline estimation: ', systime(1)-tstart 
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
