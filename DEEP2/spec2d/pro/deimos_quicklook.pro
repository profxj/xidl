 ;+
; NAME:
;   deimos_quicklook
;
; PURPOSE:
;   rectify, find wavelength solution for flat and arc, save solutions 
;           for limited list of slitlets
;
; CALLING SEQUENCE:
;   deimos_spslit, chipno, maskno, arc_header, flatimage, flativar, $
;         arcimage, arcivar, slitcoords, $
;        lamps,  slitlet_list, outdir=outdir, [plot=plot], $
;          [ybin=ybin, anamorph=anamorph, flat=flat, spline=spline]
;
; INPUTS:
;    chipno -- which of the 8 chips (1-8)
;    maskno --
;    flatimage -- spectral flat for this mask
;    flativar  -- invvar of flat
;    arc_header -- FITS header for arc spectrum
;    arcimage -- 2d arc for chip
;    arcivar -- invvar of arc image
;    slitcoords -- structure detailing bluslit information
;    lamps  -- structure detailing spectral linelist
;    slitlet_list -- which slitlets shall be analyzed
;
;    
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;   ybin -- (default 8) binning of raw spectrum for traceset
;   anamorph -- (default 1.6) anamorphic factor
;   flat   -- if set, normalize by flat
;   spline -- if set, do lambda fit by splines, rather than polynomial
;   outdir -- full path for output directory
;   plot  -- set to generate test plots on screen
;   pspath -- if set, path for ps files 
;
; OUTPUTS:
;   generates one save-set for each call, with information of all slitlets 
;   chosen for analysis.  Useful to prepare for quicklook analysis.
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
;       17Apr02 -- variation of deimos_spslit (MD)
;
;
;----------------------------------------------------------------------

pro deimos_quicklook, chipno, maskno, arc_header, $
          flatimage, flativar, arcimage, arcivar,  $
          slitcoords, lamps, slitlet_list, $
          outdir=outdir, plot=plot, pspath=pspath,  $
          ybin=ybin, anamorph=anamorph, flat=flat, spline=spline

 if NOT keyword_set(anamorph) then anamorph=1.6

;----------------------------------------------------------------------
; first deal with the tracing of the slitlets
  deimos_traceall, chipno, maskno, arc_header, flatimage, flativar, $
          slitcoords,  nslitlets, xpos1, xpos2, ypos1, ypos2,indblu, $
          ybin=ybin, plot=plot, pspath=pspath  
         
;----------------------------------------------------------------------

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



; get optical model estimates for all slitlets in question
  model_lambda =  deimos_omodel(chipno, slitcoords, arc_header)

  vprint, 2, 'Beginning loop over slits'
;TBD - look over only selected slitlets

nslits = n_elements(slitlet_list)

  for ii=0, nslits-1 do begin ;loop over all slitlets
     slitno = slitcoords[indblu[ii]].slitno  ; needs work~!!
     tstart = systime(1) 

     x0 = xpos1[*, slitno] 
     x1 = xpos2[*, slitno] 

     if median(x1-x0) LT 2 then message, 'Bad slit!!'

; rectify arc
     rect_arc = deimos_rectify_slit(arcimage, arcivar, x0, x1, /interp, $
                 xshift=xshift, npad=0, mask=arcmask) ;don't /recen

     rect_arcivar = deimos_rectify_slit(arcivar, arcivar, x0, x1, /interp, $
                                        xshift=xshift, npad=0)

     
;     rect_spec = deimos_rectify_slit(specimage, arcivar, x0, x1, /interp, $
;                 xshift=xshift, npad=0, mask=arcmask) ;don't /recen
;
;     rect_specivar = deimos_rectify_slit(specivar, arcivar, x0, x1, /interp, $
                                        xshift=xshift, npad=0)
     
     sizey = (size(rect_arc, /dimens))[1]
     

     if keyword_set(flat)  then begin 
; extract corresponding flat 
        rect_flat = deimos_rectify_slit(flatimage, flativar, x0, x1,/interp,$
                             xshift=xshift, npad=0)

;        deimos_slit_function, rect_flat, slitfn, slitfnsig, noflat
        rect_spec =  deimos_flatfield(rect_spec,  rect_flat, flat2d, $
               flat1d=slitfn, invvar=rect_specivar,  /twod, mask=noflat)
     endif

;--------- Fit for wavelengths across slitlet--------------

     deimos_fitwavelength, rect_arc,rect_arcivar,lamps,chipno, $
          slitcoords,wave, $
          spline=spline, flat=flat



; Generate spSlit structure for output
     spSlit = {flux: rect_arc, $
               ivar: rect_arcivar, $
               lambda: wave, $
               skymodel: sky, $
               mask: bitmask, $
               slitfn: slitfn}
     
; add version to header
     vers = spec2d_version()
     fxbhmake, hdr, 1, 'spSlit'
     sxaddpar, hdr, 'SP2DVERS', vers, 'Version of spec2d'
     sxaddpar, hdr, 'AUTHOR', 'Finkbeiner & Davis'



; use slit number from blueprint structure (slitcoords) for filename
     bluslitno = slitcoords[indblu[slitno]].slitno
     slitstr = string(bluslitno, format='(I3.3)')
     maskstr = string(maskno, format='(I4.4)')

     fname = outdir+'spSlit.'+maskstr+'.'+slitstr+'.fits'
     vprint, 3, 'Writing: ', fname
     mwrfits, spSlit, fname, hdr, /create
     
     vprint, 3, systime(1)-tstart,  ' s -- slitlet ', bluslitno

  endfor
  
  
  
end

