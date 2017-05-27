;+
; NAME:
;   deimos_reduce_chip
;
; PURPOSE:
;   rectify, remove sky for every slitlet on a given chip/exposure
;
; CALLING SEQUENCE:
; 
; INPUTS:
;
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
;   Data arrays are passed to this routine; it knows nothing about 
;    where files are located.
;
; REVISION HISTORY:
;
;       Wed Feb 20 17:52:02 2002, Douglas Finkbeiner (dfink)
;             Split of from 2dtest, from 2001-Aug-22
;
;
;----------------------------------------------------------------------

pro deimos_reduce_chip, chipno, maskno, arc_header, $
          flatimage, flativar, arcimage, arcivar, specimage, specivar, $
          pixmask, maskdef, objectcat, slitcoords, lamps, $
          ybin=ybin, anamorph=anamorph, flat=flat, spline=spline

;, ext1, ext2,  flat=flat

; -------- Set defaults:
  if NOT keyword_set(ybin) then ybin = 8
  if NOT keyword_set(anamorph) then anamorph = 1.6
  plot = 1B

; -------- Find some traces - ncoeff=4 means third order fit 
  deimos_trace_crude, flatimage, tset1, tset2, flatbin, ybin=ybin, ncoeff=4, $
    invvar=ivar

; -------- weed out the unpaired traces
  clean_tset, tset1
  clean_tset, tset2
  match_tset, tset1, tset2

; -------- compute tracesets for top and bottom of slits
  delvarx, ypos1
  delvarx, ypos2
  traceset2xy, tset1, ypos1, xpos1
  traceset2xy, tset2, ypos2, xpos2
  nslitlets = (size(xpos1, /dimens))[1]

; -------- debugging plots
  ny = 4096
  if keyword_set(plot) then begin 
     yind = lindgen(ny/ybin)*ybin
     atv, flatbin, min=0, max=60000
     atvplot, xpos1[yind, *], ypos1[yind, *]/ybin, ps=3,color=4
     atvplot, xpos2[yind, *], ypos2[yind, *]/ybin, ps=3,color=6
  endif 

; -------- postscript plots (for paper)
  if keyword_set(pspath) then begin 
     dfpsplot, pspath+'flat1.ps', /sq, /color, bits=8
     imps = flatbin[0:ny/ybin-1, *]
     display,bytscl(imps,min=-100,max=65000), xtit='spatial [pix]', ytit=$
       'lambda [pix]', chars=1.5, xmargin=[7, 2]
     nline = (size(xpos1))[2]
     for i=0, nline-1 do oplot, xpos1[yind, i], ypos1[yind, i]/ybin,color=4
     for i=0, nline-1 do oplot, xpos2[yind, i], ypos2[yind, i]/ybin,color=6
     dfpsclose
  endif 


; slit tops:
; note that deimos_slit_id should actually refer to objectcat and
; bluslits, not the testmask.xxxx.fits file used to construct the
; mask! TBD !!

  xfound = xpos1[2048, *]

; this is simple-minded; probably want to cross-correlate???
  deimos_slit_id, maskdef, xfound, indtop, xtop, /top, xoffs = xoffs, $
    pltscl=pltscl

  if keyword_set(plot) then begin 
     atvplot, xtop, fltarr(n_elements(xtop))+2048/ybin, ps=7, thick=2, syms=1.5
  endif 

  info = string('Chipno:', chipno, ' Top Scale:', pltscl, $
                ' (pix/asec)  Xoffs:', xoffs, format='(A,I2,A,F8.4,A,F10.3)')
  print, info

;TBD all this needs to be cleaned up!
  objno_selected = maskdef[indtop].objno 
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
  deimos_slit_id, maskdef, xfound, indbot, xbot, /bot, xoffs = xoffs, $
     pltscl=pltscl
  print, 'Chipno:', chipno,' Bot Scale:', pltscl,' (pix/asec)  Xoffs:',xoffs, $
    format='(A,I2,A,F8.4,A,F10.3)'

  if keyword_set(plot) then begin 
     atvplot, xbot, fltarr(n_elements(xbot))+2048/ybin, ps=6, thick=2, syms=1.5
     atv, arcimage
     atvplot, xpos1, ypos1, ps=3, color=4
     atvplot, xpos2, ypos2, ps=3, color=6
  endif 

; Postscript Plot 3 - Arcs
  if keyword_set(pspath) then begin 
     dfpsplot, pspath+'arc1.ps', /square, /color, bits=8
     imps = arcimage[0:511, 0:511]
     display,bytscl(imps,min=-10000,max=65000), xtit='spatial [pix]', ytit=$
       'lambda [pix]', chars=1.5, xmargin=[7, 2]
     nline = (size(xpos1))[2]
     for i=0, nline-1 do oplot, xpos1[*, i], ypos1[*, i],color=4
     for i=0, nline-1 do oplot, xpos2[*, i], ypos2[*, i],color=6
     dfpsclose
  endif


  lampwave  = lamps.lambda
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
     extarcvert = deimos_rectify_slit(arcimage, arcivar, x0, x1, /interp, $
                                      xshift=xshift, npad=0) ;don't /recen
     extivar = deimos_rectify_slit(arcivar, arcivar, x0, x1, /interp, $
                                   xshift=xshift, npad=0)
     
     sizex = (size(extarcvert, /dimens))[0]
     
; extract corresponding flat 
     
     if keyword_set(flat)  then begin
        extflatvert = deimos_rectify_slit(im, invvar, x0, x1, /interp, $
                                          xshift=xshift, npad=0)
        extarcvert = deimos_flatfield(extarcvert,  extflatvert, $
                                      flat1d=outflat, invvar=extivar,  /twod)
     endif
     
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
     
     slitpa_in = (atan(tan(maskdef[indtop[slitno]].slitpa/!radeg)/anamorph) - $
                  atan(2.*maskdef[indtop[slitno]].maskx/480.*(40./480.)))*!radeg 
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
     


     slitstr = string(maskdef[indtop[slitno]].slitn, format='(I3.3)')
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

     print, systime(1)-tstart,  ' seconds to process slitlet ', slitno

  endfor
  
  
  
end

