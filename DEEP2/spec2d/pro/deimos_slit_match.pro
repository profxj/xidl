;+
; NAME:
;   deimos_slit_match
;
; PURPOSE:
;   Associate slits found in flats with mask blueprint
;
; CALLING SEQUENCE:
;   deimos_slit_match, chipno, slitcoords, tset1, tset2, indblu, $
;    xoffs = xoffs, pltscl=pltscl, plot=plot,badmatch=badmatch
; 
; INPUTS:
;   chipno   - (1-8) which chip is this?
;   xblu1,2  - x [mm] from mask blueprint (start=1 or end=2)
;   tset1,2  - tracesets for slit boundaries in 2-D image. 
;
; OUTPUTS:
;   indblu   - index of given tset in bluprint
;
; KEYWORDS:
;   xoffs
;   pltscl  - plate scale (hardwired ???)
;   plot    - make plots?
;   badmatch - set to 1 on output if the slit tracing was problematic
;   model    - optical model structure (input)
;
; EXAMPLES:
;
; COMMENTS:
;   

;   Use discrete cross-correlation to find correct spatial offset. 
;   This won't work if the platescale is very far off!
;
; REVISION HISTORY:
;
;       Sat Feb 23 13:02:26 2002, Douglas Finkbeiner (dfink)
;		Now we delete duplicate matches
;
;       Thu Feb 21 09:17:23 2002, Douglas Finkbeiner (dfink)
;		Major restructuring, commenting. 
;
;       Thu Aug 23 09:17:41 2001, Douglas Finkbeiner (dfink)
;		First draft
;
;  2002-Oct-12  major changes, now fit for found X position [pixels] as a
;  function of xblu and yblu (in mm).  xblu and yblu are passed to
;  subroutines.  deimos_slit_coeff_eval contains everything needed to
;  evaluate the coeff array.  This has reduced the RMS error in trace
;  edge fitting from a few pixels to a few tenths.  - D. Finkbeiner
;
;  2002-Oct-14  added synth array so we can track which slit edges
;                were synthesized by deimos_slit_match_fix - DPF
;----------------------------------------------------------------------
function deimos_slit_coeff_eval, x, y, coeff

; used to be
;   a = coeff[0]+coeff[1]*xx

  a = coeff[0]+coeff[1]*x+coeff[2]*y


  return, a
end



pro deimos_slit_match_sort, ind, tset, synth

  sind = sort(ind)

  sset = {func:tset.func, xmin:tset.xmin, xmax:tset.xmax, $
          coeff: tset.coeff[*, sind]}
  tset = sset
  ind = ind[sind]
  synth = synth[sind]

  return
end



pro deimos_slit_match_sub, xblu, yblu, tset, ind, coeff, sigres=sigres
;----------------------------------------------------------------------
;  lower_break = [-5000., -3000., -1000., 1000.] ;pixel limits per quadrant
;  upper_break = [-1000.,  1000., 3000., 5000.]
; -------- compute tracesets for top (or bottom) of slits
  delvarx, xpos, ypos
  traceset2xy, tset, ypos, xpos
  npix = (size(xpos, /dimens))[0]
  ntrace = (size(xpos, /dimens))[1]


; I get tired of remembering what pltscl means!
  scale_pix_per_asec = 8.52
  scale_mm_per_asec = 0.73
  scale_mm_per_pix = scale_mm_per_asec / scale_pix_per_asec 
  print, 'using ', scale_mm_per_pix, ' mm per pixel'

; Convert blueprint values to pixel values
  xblupix = xblu/scale_mm_per_pix

; where this should appear in the image according to blu, in pixels
;  xblupix = xblu/scale_mm_per_pix
; N.B. change sign!!
;  ilow =  where(xblupix gt lower_break[quadrant])
;  ihigh = where(xblupix lt upper_break[quadrant])
;  i1 = ilow[0]                    ;first 
;  i2 = ihigh[n_elements(ihigh)-1] ; and last slitlets to search over

i1 = 0
i2 = n_elements(xblupix)-1 

;test
;print, i1, i2, ' limits for search over slitlets'

  found_x = reform(xpos[npix/2, *], ntrace)
; get matches between found slits and blueprint


  ind = discrete_correlate_match(found_x, xblupix[i1:i2], step=5)+i1
; ind has same dimension as found_X. (so does yfit)
; One element of xblupix is assigned via ind (repeats possible) for
; EACH found_x. 


; Linear fit to tsets
;  coeff = linfit(xpos[npix/2, *], xblupix[ind])
; 3-sigma clip
  sigrej = 3
;  poly_iter, found_x, xblupix[ind], 1, sigrej, yfit, coeff=coeff
; out of spec - throw away
; try multilinear regression - DPF Oct 9 02
;  now that we have the index array (ind) we should fit as a function
;  of xmm and ymm mask design params. 
; xy = transpose([[found_x], [slitcoords[ind].ymm]])
;  result = regress(xy, xblupix[ind], yfit=yfit, const=const)


; do this with multilinear regression with outlier rejection

  weights = ind*0.+1.

  if djsig(yblu[ind]) gt 5. then begin
     onexy = transpose([[fltarr(n_elements(ind))+1], [xblu[ind]], [yblu[ind]]])
     hogg_iter_linfit, onexy, found_x, weights, coeff
  endif else begin
     onexy = transpose([[fltarr(n_elements(ind))+1], [xblu[ind]]])
     hogg_iter_linfit, onexy, found_x, weights, coeff
	coeff=[reform(coeff),0.]
  endelse
  yfit = deimos_slit_coeff_eval(xblu[ind], yblu[ind], coeff)

;  coeff = [const, reform(result, n_elements(result)) ]

; compute residual
  res = yfit-found_x
  sigres = djsig(res, sigrej = sigrej)
  out = abs(res) GT sigrej*sigres

; -------- check for duplicate indices
  bad = bytarr(n_elements(ind))+1B
  bad[uniq(ind)] = 0B
  wbad = where(bad, nbad)
  if nbad GT 0 then begin 
     print, 'Trimming ', nbad, ' bad matches in deimos_slit_match'
     for i=0, nbad-1 do begin 
        badind = ind[wbad[i]]
        w = where(ind EQ badind)
        bad[w] = 1B ; set all with (ind eq badind) bad for a moment
;        dif = blupred[ind[w]] - xpos[npix/2, w]
        foo = min(abs(res[w]), wdif)
        bad[w[wdif]] = 0B  ; set this one not bad
     endfor 
     bad = bad OR out
     good = where(bad EQ 0)
     ind = ind[good]
     sset = {func:tset.func, xmin:tset.xmin, xmax:tset.xmax, $
             coeff: tset.coeff[*, good]}
     tset = sset
  endif 

  return
end

pro deimos_slit_match_model, xblu, tset, ind, coeff, sigres=sigres
;----------------------------------------------------------------------
;  lower_break = [-5000., -3000., -1000., 1000.] ;pixel limits per quadrant
;  upper_break = [-1000.,  1000., 3000., 5000.]
; -------- compute tracesets for top (or bottom) of slits
  delvarx, xpos, ypos
  traceset2xy, tset, ypos, xpos
  npix = (size(xpos, /dimens))[0]
  ntrace = (size(xpos, /dimens))[1]


; I get tired of remembering what pltscl means!
;  scale_pix_per_asec = 8.52
;  scale_mm_per_asec = 0.73
;  scale_mm_per_pix = scale_mm_per_asec / scale_pix_per_asec 
;  print, 'using ', scale_mm_per_pix, ' mm per pixel'

; Convert blueprint values to pixel values
;  xblupix = xblu/scale_mm_per_pix

  xblupix=xblu

; where this should appear in the image according to blu, in pixels
;  xblupix = xblu/scale_mm_per_pix
; N.B. change sign!!
;  ilow =  where(xblupix gt lower_break[quadrant])
;  ihigh = where(xblupix lt upper_break[quadrant])
;  i1 = ilow[0]                    ;first 
;  i2 = ihigh[n_elements(ihigh)-1] ; and last slitlets to search over

;i1 = 0
;i2 = n_elements(xblupix)-1 

i1=min(where(xblu ne -1.))
i2=max(where(xblu ne -1.))


;test
;print, i1, i2, ' limits for search over slitlets'

  found_x = reform(xpos[npix/2, *], ntrace)
; get matches between found slits and blueprint


  ind = discrete_correlate_match(found_x, xblupix[i1:i2], step=1,lag=[-50,50])+i1
; ind has same dimension as found_X. (so does yfit)
; One element of xblupix is assigned via ind (repeats possible) for
; EACH found_x. 


; Linear fit to tsets
;  coeff = linfit(xpos[npix/2, *], xblupix[ind])
; 3-sigma clip
  sigrej = 3
;  poly_iter, found_x, xblupix[ind], 1, sigrej, yfit, coeff=coeff
; out of spec - throw away
; try multilinear regression - DPF Oct 9 02
;  now that we have the index array (ind) we should fit as a function
;  of xmm and ymm mask design params. 
; xy = transpose([[found_x], [slitcoords[ind].ymm]])
;  result = regress(xy, xblupix[ind], yfit=yfit, const=const)


; do this with multilinear regression with outlier rejection

  weights = ind*0.+1.

;  if djsig(yblu[ind]) gt 5. then begin
;     onexy = transpose([[fltarr(n_elements(ind))+1], [xblu[ind]], [yblu[ind]]])
;     hogg_iter_linfit, onexy, found_x, weights, coeff
;  endif else begin
     onexy = transpose([[fltarr(n_elements(ind))+1], [xblu[ind]]])
     hogg_iter_linfit, onexy, found_x, weights, coeff
	coeff=[reform(coeff),0.]
;  endelse
  yfit = deimos_slit_coeff_eval(xblu[ind], xblu[ind]*0., coeff)

;  coeff = [const, reform(result, n_elements(result)) ]

; compute residual
  res = yfit-found_x
  sigres = djsig(res, sigrej = sigrej)
  out = abs(res) GT sigrej*sigres



; -------- check for duplicate indices
  bad = bytarr(n_elements(ind))+1B
  bad[uniq(ind)] = 0B
  wbad = where(bad, nbad)
  if nbad GT 0 then begin 
     print, 'Trimming ', nbad, ' bad matches in deimos_slit_match'
     for i=0, nbad-1 do begin 
        badind = ind[wbad[i]]
        w = where(ind EQ badind)
        bad[w] = 1B ; set all with (ind eq badind) bad for a moment
;        dif = blupred[ind[w]] - xpos[npix/2, w]
        foo = min(abs(res[w]), wdif)
        bad[w[wdif]] = 0B  ; set this one not bad
     endfor 
     bad = bad OR out
     good = where(bad EQ 0)
     ind = ind[good]
     sset = {func:tset.func, xmin:tset.xmin, xmax:tset.xmax, $
             coeff: tset.coeff[*, good]}
     tset = sset
  endif 

  return
end


;+
; NAME:
;   deimos_slit_match_fix
;
; PURPOSE:
;   Repair traceset for blueprint slits not found in traceset
;
; CALLING SEQUENCE:
;   deimos_slit_match_fix, need, ind, tset, xblu, yblu, coeff, synth
;
; INPUTS:
;   need    - boolean: for each slit in blueprint, is it on this chip?
;   ind     - blueprint index corresponding to each found slit
;   tset    - traceset for found slits (same dim as ind)
;   xblu    - blueprint X position for appropriate end of slit
;   yblu    - blueprint Y position
;   coeff   - coefficients of bilinear slit edge position model
; 
; OUTPUTS:
;   ind     - same as input, but augmented with interpolated slits
;   tset    - fixed up tset (corresponding to ind)
;   synth   - boolean: 0=slit was found, 1=slit synthesized by this
;             routine 
;
; COMMENTS:
;   BOT         - Written by D. Finkbeiner
;   2002-Oct-14 - added synth array - DPF
;
pro deimos_slit_match_fix, need, ind, tset, xblu, yblu, coeff, synth

; -------- Remove anything we have and don't need
  good = where(need[ind], ngood)
  ind = ind[good]
  sset = {func:tset.func, xmin:tset.xmin, xmax:tset.xmax, $
          coeff: tset.coeff[*, good]}
  tset = sset
  synth = bytarr(n_elements(ind)) 

; -------- mark off the ones we already have
  needadd = need
  needadd[ind] = 0B
  
  needind = where(needadd, nneed) ; slits we still need to put in tset
  if nneed GT 0 then begin 
     xpred = deimos_slit_coeff_eval(xblu, yblu, coeff)
     newcoeff = tset_coeff_interpol(tset, xpred[ind], xpred[needind])

     sset = {func:tset.func, xmin:tset.xmin, xmax:tset.xmax, $
              coeff:([[tset.coeff],[newcoeff]])}
     ind = [ind, needind]
     tset = sset
     synth = [synth, bytarr(nneed)+1B]
  endif

  return
end



; main routine:
pro deimos_slit_match, chipno, tset1, tset2, indblu, synth1, synth2, $
          xoffs = xoffs, plot=plot, slitcoords=slitcoords, $
                       badmatch=badmatch, model=model

; remove duplicates and anything too far away from expected place
; maybe pass yblu instead of slitcoords???
; maybe yblu1 and yblu2 ???

; -------- use information from slitcoords (from blueprint) to match
;           found slit positions with slit numbers.             

; convert these to x-positions outside of deimos_slit_match
  xblu1 = slitcoords.xmm - slitcoords.slitlength/2
  xblu2 = slitcoords.xmm + slitcoords.slitlength/2
  yblu1 = slitcoords.ymm - slitcoords.slitlength*tan(slitcoords.pa/!radeg)/2
  yblu2 = slitcoords.ymm + slitcoords.slitlength*tan(slitcoords.pa/!radeg)/2

  usemodel=0
  nmod=total(model.xt gt 0 OR model.xb gt 0)
  if total(nmod) lt 7 then begin
      usemodel=1
      deimos_slit_match_model,model.xb,tset1, ind1, coeff1, sigres=sigres1
      deimos_slit_match_model,model.xt,tset2, ind2, coeff2, sigres=sigres2
  endif else begin
      deimos_slit_match_sub, xblu1, yblu1, tset1, ind1, coeff1, sigres=sigres1
      deimos_slit_match_sub, xblu2, yblu2, tset2, ind2, coeff2, sigres=sigres2
  endelse
  badmatch=0



; Check sigres (sigma of residual of fit)
  print
  print, 'SLIT_MATCH: RMS residuals for xpos1, xpos2: ', sigres1, ',', $
    sigres2, ' pixels', format='(A,2(F6.2,A))'
  if usemodel eq 0 then $
    if (sigres1 > sigres2) GT mean(abs(xblu2-xblu1)*coeff1[1])/5 then begin
      message, 'Warning:  Slit matching may be totally wrong!!', /info
      badmatch=1
    endif

  if keyword_set(plot) then begin 
     delvarx, ypos1
     delvarx, ypos2
     traceset2xy, tset1, ypos1, xpos1
     traceset2xy, tset2, ypos2, xpos2
     npix = (size(xpos1, /dimens))[0]
     splot,xpos1[npix/2, *], xblu1[ind1],ps=-4,/yno
     soplot,ps=-4,xpos2[npix/2, *], xblu2[ind2]
     xx = findgen(2000)
;     soplot, xx, deimos_slit_coeff_eval(xx, 0, coeff1)   ; no longer relevent
;     soplot, xx, deimos_slit_coeff_eval(xx, 0, coeff2), line=2

  endif 

; -------- Fix up bad/missing traces. 
;  blulim1 = deimos_slit_coeff_eval(0, 0, coeff1) < deimos_slit_coeff_eval(0, 0, coeff2) 
;  blulim2 = deimos_slit_coeff_eval(2047., 0, coeff1) < deimos_slit_coeff_eval(2047., 0, coeff2) 

  if usemodel then begin
      xblu1=model.xb
      yblu1=xblu1*0+1.
      xblu2=model.xt
      yblu2=xblu2*0+1.
  endif



  xpred1 = deimos_slit_coeff_eval(xblu1, yblu1, coeff1)
  xpred2 = deimos_slit_coeff_eval(xblu2, yblu2, coeff2)



;  blulim2 = (coeff1[0]+coeff1[1]*2047.) > (coeff2[0]+coeff2[1]*2047.)

; need exactly one occurance of each index in "need"
;  need = (xblupix1 GT blulim1) AND (xblupix2 LT blulim2)
  buffer = 20.   ;set this to a bigger number if you want bigger slits
;  need = (xpred1 GT (0-buffer) ) AND (xpred2 LT (2047+buffer)) ; if central row on chip, we keep
  need = (xpred2 gt buffer  and xpred1 lt 2047.-buffer)

; at least one edge must be on chip, and at least buffer pixels must
; be on chip.

  deimos_slit_match_fix, need, ind1, tset1, xblu1, yblu1, coeff1, synth1
  deimos_slit_match_fix, need, ind2, tset2, xblu2, yblu2, coeff2, synth2

; -------- match tset1 and tset2 by sorting
  deimos_slit_match_sort, ind1, tset1, synth1
  deimos_slit_match_sort, ind2, tset2, synth2

  indblu = ind1


  return
end






