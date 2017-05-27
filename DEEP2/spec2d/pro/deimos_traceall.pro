;+
; NAME:
;   deimos_traceall
;
; PURPOSE:
;   find traces of slitlets in a chip
;
; CALLING SEQUENCE:
;   deimos_traceall, chipno, flatimage, flativar, $
;         slitcoords, nslitlets, xpos1, xpos2, ypos1, ypos2, indblu, $
;         ybin=ybin, plot=plot, pspath=pspath  ,badmatch=badmatch, $
;         model=model
;
; INPUTS:
;    chipno -- which of the 8 chips (1-8)
;    flatimage -- spectral flat for this mask
;    flativar  -- invvar of flat
;    slitcoords -- structure detailing bluslit information
;    
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;   ybin -- (default 8) binning of raw spectrum for traceset
;   plot  -- set to generate test plots on screen
;   pspath -- if set, path for ps files 
;   model -- optical model structure
;
; OUTPUTS:
;   nslitlets -- how many slitlets on this chip
;   xpos1 -- xposition of start edge of traces
;   xpos2 -- xposition of end edge of traces
;   ypos1 -- yposition of start edge of traces
;   ypos2 -- yposition of end edge of traces
;   synth1 -- array[nslit] =1 if xpos1 was synthesized from neighbors
;   synth2 -- array[nslit] =1 if xpos2 was synthesized from neighbors
;   indblu -- index of structure detailing slitlets, from bluprint FITS header
;   badmatch - set to 1 on output if the slit tracing was problematic
; OPTIONAL OUTPUTS:
;
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
;   split off from deimos_spslit because this code will be called from both
;    quick-look and final spec2d analysis
;
; REVISION HISTORY:
;
;       Wed Feb 20 17:52:02 2002, Douglas Finkbeiner (dfink)
;             Split of from 2dtest, from 2001-Aug-22
;       2002-apr-13  MD - split from deimos_spslit 
;       2002-apr-25  DPF - cleaned up a bit more
;       2002-oct-14  DPF - removed flipsign and other deadwood
;
;----------------------------------------------------------------------
pro deimos_traceall, chipno, flatimage, flativar, $
          slitcoords, nslitlets, xpos1, xpos2, ypos1, ypos2, $
          synth1, synth2, indblu, $
          ybin=ybin, plot=plot, pspath=pspath, badmatch=badmatch  , $
                     model=model
  
; -------- Set defaults:
  if NOT keyword_set(ybin) then ybin = 8
  if NOT keyword_set(flativar) then begin 
     print, 'DEIMOS_TRACEALL:  You should really set flativar!'
  endif else ivar = flativar NE 0  ; treat as mask (1=good) for now

; -------- Find some traces - ncoeff=4 means third order fit 
  deimos_trace_crude, flatimage, tset1, tset2, flatbin, ybin=ybin, ncoeff=4, $
    invvar=ivar

; -------- reject traces with strange (don't match neighbors) coefficients
  ntset=(size(tset1.coeff, /dimens))[1] < (size(tset2.coeff, /dimens))[1]

  if ntset gt 5 then begin
      clean_tset, tset1
      clean_tset, tset2
  endif

; -------- compute tracesets for top and bottom of slits
  delvarx, ypos1
  delvarx, ypos2
  traceset2xy, tset1, ypos1, xpos1
  traceset2xy, tset2, ypos2, xpos2
  nslitlets = (size(xpos1, /dimens))[1]

; -------- debugging plots
  ny = 4096
  yind = lindgen(ny/ybin)*ybin
  if keyword_set(plot) then begin 
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
     nline1 = (size(xpos1))[2]
     nline2 = (size(xpos2))[2]
     for i=0, nline1-1 do oplot, xpos1[yind, i], ypos1[yind, i]/ybin,color=4
     for i=0, nline2-1 do oplot, xpos2[yind, i], ypos2[yind, i]/ybin,color=6
     dfpsclose
  endif 

; -------- use information from slitcoords (from blueprint) to match
;           found slit positions with slit numbers.             
;          This overwrites tset1, tset2 to include anything we missed. 

  deimos_slit_match, chipno, tset1, tset2, indblu, $
    synth1, synth2, slitcoords=slitcoords, badmatch=badmatch, model=model

; indblu is index in bluslit (or slitcoords) for each slit found

; -------- compute tracesets for top and bottom of slits (again)
  delvarx, ypos1
  delvarx, ypos2
  traceset2xy, tset1, ypos1, xpos1
  traceset2xy, tset2, ypos2, xpos2
  nslitlets = (size(xpos1, /dimens))[1]

  if keyword_set(plot) then begin 
     yind = lindgen(ny/ybin)*ybin
     atvplot, xpos1[yind, *], ypos1[yind, *]/ybin, ps=3,color=5
     atvplot, xpos2[yind, *], ypos2[yind, *]/ybin, ps=3,color=7
     ws1 = where(synth1, nws1)
     if nws1 NE 0 then atvplot, xpos1[*, ws1], ypos1[*, ws1]/ybin, ps=3,color=1
     ws2 = where(synth2, nws2)
     if nws2 NE 0 then atvplot, xpos2[*, ws2], ypos2[*, ws2]/ybin, ps=3,color=2
     abandy = 300-(slitcoords[indblu].ymm*.87)
     for i=0, nslitlets-1 do begin 
        ii = i*12+100
        yband = abandy[i]*ybin
        atvplot, [xpos1[yband, i], xpos2[yband, i]], $
          [1, 1]*abandy[i], thick=5, color=3
        atvxyouts, 0.5*(xpos1[yband, i]+xpos2[yband, i]), abandy[i]+8, $
          string(indblu[i], format='(I3.3)'), align=0.5, charsize=1.8, $
          color=2
     endfor
  endif 

  return
end
