;+
; NAME:
;   deimos_slit_id
;
; PURPOSE:
;   Associate slits found in flats with mask blueprint
;
; CALLING SEQUENCE:
; 
; INPUTS:
;   mask     - mask structure from mask definition file
;   xfound   - x (spatial coord) [pix] where an edge was found
;   
; OUTPUTS:
;   ind
;   xslit
;
; KEYWORDS:
;   top     - set of xfound refers to top of slitlet
;   bot     - not used ???
;   xoffs
;   pltscl  - plate scale (hardwired ???)
;   plot    - make plots?
;
; EXAMPLES:
;
; COMMENTS:
;   Use discrete cross-correlation to find correct spatial offset. 
;   This won't work if the platescale is very far off!
;
; REVISION HISTORY:
;
;       Thu Feb 21 09:17:23 2002, Douglas Finkbeiner (dfink)
;		Major restructuring, commenting. 
;
;       Thu Aug 23 09:17:41 2001, Douglas Finkbeiner (dfink)
;		First draft
;
;----------------------------------------------------------------------
pro deimos_slit_id, maskdef, xfound, ind, xslit, top=top, bot=bot, $
  xoffs = xoffs, pltscl=pltscl, plot=plot

  nfound = n_elements(xfound)

; hardwire the plate scale guess (pix / asec)
  if NOT keyword_set(pltscl) then pltscl = 8.52

; top/bottom of each slitlet
  if keyword_set(top) then begin 
     xedge = maskdef.maskx+mask.topdist
  endif else begin
     xedge = maskdef.maskx-maskdef.botdist
  endelse

; make a prediction where this should appear in the image, in pixels
  xpred = -xedge*pltscl    ; fit for offset later

; get offset in pixels between xfound and xpred
  xoffs = deimos_slitpos_offset(xfound, xpred)
  xpred = xpred-xoffs

; for each FOUND slit, get the nearest predicted slit position
  ind = lonarr(nfound)
  for i=0, nfound-1 do ind[i] = where(min(abs(xfound[i]-xpred)) EQ $
                                      abs(xfound[i]-xpred))
  dx = xfound-xpred[ind]
  if keyword_set(plot) then begin 
     splot, xfound, dx, ps=1, xtit='x found', ytit='delta x, found-predicted'
  endif 

; Pass 2 remove linear trend - i.e. adjust plate scale
  poly_iter, xfound, dx, 1, 3, xfit, coeff=coeff
;  lf    = linfit(xfound, dx)
  pltscl = pltscl*(1+coeff[1])
  xpred = -xedge*pltscl
  xoffs = deimos_slitpos_offset(xfound, xpred)
  xpred = xpred-xoffs
;  print, 'Plate Scale ', pltscl, ' pixels per arcsecond'
  for i=0, nfound-1 do ind[i] = where(min(abs(xfound[i]-xpred)) EQ $
                                      abs(xfound[i]-xpred))
  dx = xfound-xpred[ind]
  if keyword_set(plot) then begin 
     splot, xfound, dx, ps=1, xtit='x found', ytit='delta x, found-predicted' 
  endif 

; Pass 3 tweak offset 
  print, 'Tweak offset', median(dx), ' pixels'
  xpred = xpred+median(dx)
  for i=0, nfound-1 do ind[i] = where(min(abs(xfound[i]-xpred)) EQ $
                                      abs(xfound[i]-xpred))

  xslit = xpred[ind]
  dx = xfound-xslit

  poly_iter, maskdef[ind].masky, dx, 1, 3, xfit, coeff=lf
;  lf = linfit(maskdef[ind].masky, dx, yfit=xfit)

;  xslit = xslit+xfit
;  dx = xfound-xslit
  if keyword_set(plot) then begin 
     splot, xslit, dx, ps=1,  xtit='slit x position',  ytit='delta, found-predicted'
  endif 
  print, 'linear fit for y displacement', lf[1], ' sig', djsig(dx)

  return
end
