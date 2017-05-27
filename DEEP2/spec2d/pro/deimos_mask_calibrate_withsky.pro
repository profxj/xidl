;+
; NAME:
;   deimos_mask_calibrate
;
; PURPOSE:
;   top level routine to generate flats and arcs to a given mask,
;   with separate FITS files for each slitlet.
;77
; CALLING SEQUENCE:
;   deimos_mask_calibrate, planfile, chiplist=chiplist, $
;                          onlyalign=onlyalign, slitlist=slitlist
; 
; INPUTS:
;   planfile -- an ascii file detailing what flats, arcs, science
;               frames to use for this mask, as well as directory(s)
;               from which to take raw data and where to put results
; 
; KEYWORDS: 
;   chiplist -- if set, what subset of slitlets to use, otherwise all
;   slitlist -- if set, what subset of slitlets to calibrate
;               (+alignment stars), otherwise all
;            are reduced on selected chips
;   onlyalign -- if set, only alignment star slitlets are analyzed
;
; OUTPUTS:
;   separate FITS files for each slitlet, consisting of one HDU,
;   separate files for R and B sides of spectrum
;
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
;   working code, still evolving
;
; REVISION HISTORY:
;  Doug Finkbeiner, Marc Davis 3Jun2002
;  JAN 07jul02 - pix-flat normalization included
;  2002-Oct-09 - changed slitno value in header to bluslitno - DPF
;----------------------------------------------------------------------
FUNCTION get_anamorph, gratenum, lambda_c

;;; READ-IN THE DEIMOS SPECS FROM FILE.
  READCOL, djs_filepath('deimos_anamorph.dat', ROOT=GETENV('DEEP_DIR'), $
                        SUBDIR='spec2d/etc'), grating, lambda, anamor, $
    FORMAT='F,F,F', COMMENT='#', /silent

;;; NOW INTERPOLATE TO GET THE CORRECT ANAMORPH VALUE.
  dex = WHERE(grating EQ gratenum, cnt)
  IF cnt GT 0 THEN anamorph = INTERPOL(anamor[dex], lambda[dex], lambda_c) $
    ELSE MESSAGE, 'ERROR: Data for ' + STRCOMPRESS(STRING(gratenum), /REMOVE_ALL) + $
    ' grating not in file deimos_anamorph.dat!'

;;; RETURN VALUE.
  RETURN, anamorph
END

; a procedure to find the slit pa and mask pa for a given slit by
; sifting thru the bintab file for the needed information. slitno is
; assumed to be of type LONG and slitpa and maskpa will be returned as
; type FLOAT.
pro get_pa, slitno, slitpa=slitpa, maskpa=maskpa, $
            desitab=desitab, masktab=masktab

  if keyword_set(desitab) then begin
;  if keyword_set(masktab) and keyword_set(desitab) then begin
; cross-reference the two tables and get the correct slitpa and maskpa
; values.
      slitdex = where(long(desitab.slitname) eq slitno, cnt)
      if cnt eq 0 then begin
          vprint, 3,string('No entry found in DesiSlit table for slit ' + $
            string(slitno, format='(I3.3)'))
          slitpa = 999.
          maskpa = 999.
          return
      endif
      slitpa = desitab[slitdex[0]].slitlpa[0]
      if keyword_set(masktab) then begin
          desno = desitab[slitdex[0]].desid ;type LONG
          slitdex = where(long(masktab.desid) eq desno, cnt)
          if cnt eq 0 then begin
              vprint,3, $
                string('No entry found in MaskDesign table for slit ' + $
                   string(slitno, format='(I3.3)'))
              maskpa = 999.
              return
          endif 
          maskpa = float(masktab[slitdex[0]].pa_pnt[0])
      endif else maskpa = 999.0
  endif else begin
      slitpa = 999.
      maskpa = 999.
  endelse

end



pro deimos_mask_calibrate_withsky,  planfile,  chiplist=chiplist, noplot=noplot, $
     slitlist=slitlist, onlyalign=onlyalign,quick=quick

  if n_elements(quick) eq 0 then quick = 0
  if n_elements(planfile) eq 0 then planfile = '*.plan'
  if quick then noplot=1
  verbset, 4-3*(quick gt 0)

; -------- Get environment variables, set paths
  deimos_data = getenv('DEIMOS_DATA')+'/'
  calib_data = getenv('CALIB_DATA')+'/'

;  if deimos_data eq '/' then message, 'You need to set $DEIMOS_DATA!'
  if deimos_data eq '/' then deimos_data = ''
  if calib_data  eq '/' then message, 'You need to set $CALIB_DATA!'

  deep_dir = getenv('DEEP_DIR')+'/'
  if deep_dir eq '/' then message, 'You need to set $DEEP_DIR!'

; input the plan file for this mask
  read_planfile, planfile, maskname, rawdatadir, outdatadir, $
             flatnames, arcnames, sciencenames, slitnums, $
             nbright, ndim, nsimp, nlong, $
             polyflag=polyflag,$
	     chips=chiptext, bluearc = bluearc, redarc = redarc, $
                 linelist = linelist

  if quick gt 0 then quickslits, planfile, slitnums

  if n_elements(slitlist) eq 0 AND n_elements(slitnums) ne 0 then $
    slitlist=slitnums

   maskstr = size(maskname, /tname) EQ 'STRING' ? $
          maskname : string(maskname, format='(I4.4)')


; check if this is a deep mask...set isdeep=1:0
  deimos_isdeep, isdeep

; extract the design table and the mask table from the bintab file.
  bfile = findfile(strmid(maskstr,0,4)+'*bintab*.fits*', count=nfiles)
  if nfiles eq 0 then vprint, 1, 'ERROR: No bintab file found!!!' $
  else begin
      bfile = bfile[0]
      fits_open, bfile, fcb
      extnames = fcb.extname
      fits_close, fcb
      desidex = where(extnames eq 'DesiSlits', desicnt)
      maskdex = where(extnames eq 'MaskDesign', maskcnt)
      if desicnt eq 0 then $
        vprint, 1, $
           string('No DesiSlits table found in file ' + bfile + '!!!') $
      else desitab = mrdfits(bfile, desidex[0], /silent)
      if maskcnt eq 0 then $
        vprint, 1, $
        string('No MaskDesign table found in file ' + bfile + '!!!') $
      else masktab = mrdfits(bfile, maskdex[0], /silent)
  endelse

; set the read-noise value.
  readnoise = 2.32
  noisesq=readnoise^2+5. ; floor for counts in IVAR calculations
  
; call astrolib if needed
  defsysv, '!TEXTOUT', exist=exist
  if NOT exist then astrolib

;  if NOT keyword_set(chiplist) then chiplist = [1, 2, 3, 4, 5, 6, 7, 8]

  device, pseudo=8
  verbset, 4
  ybin = 8
  flat = 1
  pspath = './ps/' ;where to put ps files
  spawn, 'mkdir -p '+pspath ; make sure directory exists
  if NOT keyword_set(noplot) then plot = 1 ; make plots


; check to see if we're using different arcs for the two sides
     diffarcs = n_elements(bluearc) NE n_elements(redarc)
     if diffarcs eq 0 then diffarcs=total(bluearc eq redarc) $
       ne n_elements(bluearc)

	chipnums=fix(strsplit(chiptext,',',/extract))
   	if NOT keyword_set(chiplist) then chiplist=chipnums

;  if n_elements(sciencenames) LT 3 then message, $
;    'You must have 3 science exposures!'

  vprint, 1,string('  DEIMOS pipeline version:  ', spec2d_version())

; limit to 3 flats.  When there are more, it is usually a mistake. 
  nflats = n_elements(flatnames)
  flatlimit=10-2*(quick gt 0)

  if nflats GT flatlimit then begin 
     vprint,3,''
     vprint, 3, $
       string('You have', nflats, ' flats; I will use the last ',flatlimit, $
       format='(A,I2,A,I2)')
     flatnames = flatnames[nflats-flatlimit:nflats-1]
     nflats = n_elements(flatnames)
  endif 

  maskdir = deimos_data+ rawdatadir
  vprint,3,''
  vprint, 3,string('Directory of raw data: ', maskdir)
  if (outdatadir ne '' and outdatadir ne '.') $
                  then outdatadir = outdatadir + '/'
  vprint, 3,string('directory of output data: ',  outdatadir)

  arc_header = headfits(arcnames[0]) ;header of 1st arc frame
  model_lambda = intarr(1000)
  vprint, 3,'==========================================================='
  deimos_grating, arc_header, grating, grangle, lambda_c
  anamorph = get_anamorph(floor(grating), lambda_c) 
         ;get correct anamorphic factor
  vprint, 3, $
   string('  Grating [line/mm]:', grating, '    Grating Angle:', grangle, $
    format='(A,I6,A,F12.5)')
  vprint, 3, $
    string('  Central Wave: ', lambda_c, '    Anamorphic factor:', anamorph, $
    format='(A,F10.4,A,F8.2)')
  vprint,3,''

; Read lamp file into a structure
  lampfilename = deep_dir+'spec2d/etc/'+linelist ;lamp_NIST.dat'

  bluelamps=process_arclists_sky(bluearc,lampfilename,grating,anamorph)

  if diffarcs then $
    redlamps=process_arclists_sky(redarc,lampfilename,grating,anamorph) $
    else redlamps = bluelamps


;keep only those lines from elements turned on in lamp spectrum


; perhaps things like this should be versioned???
; Generate badmask file if it is not there already
  badmaskname = concat_dir(calib_data,'deimos_badmask.fits.Z')
  pixflatname = concat_dir(calib_data,'processed_pix_mult_flat.fits.gz')

  if file_test(badmaskname) eq 0 then begin 
     ; generate file
     vprint, 1,' Cannot find badmask file - Generating now...'
     spawn, 'mkdir -p '+calib_data
     deimos_badchip, calib_data
     vprint, 1,' Done.'
  endif

; objectcat: grab table from FITS file on flatnames
  if nfiles gt 0 then $
    objectcat = deimos_tables(bfile[0], bluslits=slitcoords) $
     else objectcat = deimos_tables(flatnames[0], bluslits=slitcoords) 

; check for bogus slits
  good = (slitcoords.slitlength ge 0.1) and (slitcoords.slitwidth ge 0.1)
  w = where(good, ngood)
  if n_elements(good) NE ngood then begin 
     slitcoords = slitcoords[w]
     vprint, 2,'Removing bad entry in slitcoords!'
     vprint, 2,string(where(good eq 0))
  endif 
  slitcoords = slitcoords[sort(slitcoords.xmm)] ;order spatially
  nslitlets = n_elements(slitcoords) 


; LOOP OVER CHIPS

  for ichip=0, n_elements(chiplist) -1  do begin 
     chipno = chiplist[ichip]
     vprint,3, '-----------------------------------------------------------------------------'
     vprint,3,''
     vprint, 1,string('Chipno: ', chipno)
     vprint,3,''
     
     chipisred = chipno GT 4

; -------- read subimage
    deimos_badchip = mrdfits(badmaskname, chipno,  /silent) 

    ; will want throw out slits that are too vignetted at the sides
    vigcolumn= total(deimos_badchip AND 2B, 2) gt 3000
    vprint,3,string( total(vigcolumn), ' columns vignetted')

; identify nearly-vignetted pixels to help flatfield code
    nearvig=(deimos_badchip AND 6b)

;    if file_test(pixflatname) eq 0 OR quick gt 0 then begin 
    if file_test(pixflatname) eq 0  then begin 
                  vprint, 1, 'Setting pixflat to 1...'
                  deimos_pixflat_mul = 1.
    endif else begin
	vprint,4,'Using pixflat!'
          deimos_pixflat_mul = mrdfits(pixflatname, chipno, /silent)+1
          whzero = where(deimos_pixflat_mul eq 0, zeroct)
          if zeroct gt 0 then deimos_pixflat_mul[whzero] = 1.
    endelse

; process optical model for this chip
     model_lambda =  deimos_omodel(chipno, slitcoords, arc_header)

     nslits=total(model_lambda.xb gt 0 AND model_lambda.xt gt 0)



; process flats
    vprint, 3,'Processing flats...'
    flatimage = deimos_read_chip(flatnames[0], chipno, satmask=satmask, $
                                 quick=quick) 
    s = size(flatimage, /dim)
    sizexf=s[0]
    sizeyf=s[1]
    sumcts = 0.
    flatsat = satmask           ; flat saturated mask
    
    nsample=1E4-1
    sampleidx=lindgen(nsample)*(sizexf*1L*sizeyf)/(nsample+7)

; making flats: ~40 sec

; combine flats with avsigclip if there are more than 1
    if nflats gt 1 and quick le 0 then begin
       flatarr  = fltarr(s[0], s[1], nflats, /nozero) ; don't zero array
       flativar = fltarr(s[0], s[1], nflats, /nozero)
       for j=0,nflats-1 do begin
          if j NE 0 then flatimage=deimos_read_chip(flatnames[j], $
                                chipno, sat=satmask, quick=quick)
          medimg = float(median(flatimage[sampleidx]))
          sumcts = sumcts+medimg
          flatarr[*, *, j]  = flatimage/medimg
          flativar[*, *, j] = medimg^2/(flatimage > noisesq)
          flatsat = flatsat OR satmask
       endfor

       delvarx,satmask


       clipflat=avsigclip(flatarr, flativar)	
       delvarx,flatarr,flativar
       flatimage = clipflat.flux*(sumcts/nflats)
       nused=fltarr(sizexf,sizeyf)

       for j=0,nflats-1 do begin
          nused=nused+( (clipflat.mask AND 2^j) eq 0b)
       endfor	

; count frames that went into pixel for inverse variance array

    endif else if nflats gt 1 then begin
         for j=0,nflats-1 do begin
          if j NE 0 then flatimage=flatimage+deimos_read_chip(flatnames[j], $
                                chipno, sat=satmask, quick=quick)
;          medimg = float(median(flatimage[sampleidx]))
;          sumcts = sumcts+medimg
;          flatarr[*, *, j]  = flatimage/medimg
;          flativar[*, *, j] = medimg^2/(flatimage > noisesq)
          flatsat = flatsat OR satmask
         endfor
         flativar = 1/(flatimage > nflats*noisesq)
         nused = 1.
    endif else nused=1.



    delvarx, clipflat           ; save memory

    flatsig = stdev(flatimage[sampleidx])

    if flatsig GT 1 AND nslits gt 1 then begin 


     if n_elements(deimos_pixflat_mul) gt 1 then begin

        pixflatgood = (deimos_pixflat_mul GT 0.8) AND $
          (deimos_pixflat_mul LT 1.25)
           
        flativar = ((deimos_badchip AND 1b) eq 0b )/(flatimage > noisesq) 


        if quick le 0 then begin
            ; correct for pix-to-pix variation
            flatimage = flatimage*deimos_pixflat_mul 
            flativar = flativar/deimos_pixflat_mul^2*pixflatgood*nused
             flatinterp = djs_maskinterp(flatimage, flativar eq 0., iaxis=0)
        endif else begin
; for quicklook, just throw away deviant pixels in pixflat
            flativar=flativar*pixflatgood*nused
             flatinterp = djs_maskinterp(flatimage, flativar eq 0., iaxis=0)
 ;           flatinterp=flatimage
        endelse
      
     endif else begin
        flativar=((deimos_badchip AND 1b) eq 0b )/(flatimage > noisesq)
        flativar = flativar*nused
        flatinterp=flatimage
     endelse
        

      delvarx,flatimage,nused

; tracing: 7 sec

      deimos_traceall, chipno, flatinterp, flativar, $
          slitcoords, nslitlets, xpos1, xpos2, ypos1, ypos2, synth1, synth2, $
          indblu, ybin=ybin, plot=plot, pspath=pspath, badmatch=badmatch, $
        model=model_lambda



      if badmatch and quick gt 0 then begin
          message,'Bad slit traces - buckled mask????',/info
          openw,2,'buckledmask.txt'
          printf,2,'This mask appears to be buckled!'
          close,2
         return
      endif    

; fix the 'glow' at 3 edges of each chip

     flatinterp=deimos_fixglow(flatinterp, $
         mask=((deimos_badchip AND 1b) EQ 1b),fixbottom=chipisred)

     if quick le 0 then begin
   
; correct for the scattered light.  This is a 6-8% effect in our flats; the below code gets it right to 5-10%
;   (so the net result is a sub-percent error in our flats/slitfunctions, as opposed to 6-8%
; the correction is done by a 2-plane fit to the light in interslit minima on each chip.
	deimos_findmins,flatinterp,(deimos_badchip AND 1b) NE 0b,planepars
        flatinterp=deimos_submins(flatinterp,planepars)       
    endif     


    if chipisred then usenames=redarc else usenames=bluearc
 
    vprint,3,'Using arcs: ',usenames

    arcfile = usenames[0]

; processing arcs: 12 sec

    arcimage = deimos_read_chip(arcfile, chipno, sat=arcsat, quick=quick)
    dirtyimage=arcimage
    if n_elements(usenames) gt 1 then begin 
        fullarc=fltarr(s[0], s[1], n_elements(usenames), /nozero)
        fullarcivar=fltarr(s[0], s[1], n_elements(usenames), /nozero)
        fullarcsat=bytarr(s[0], s[1], n_elements(usenames), /nozero)
        fullarc[*,*,0]=arcimage
        fullarcsat[*,*,0]=arcsat
        pessimsat=arcsat

       for j=1,n_elements(usenames)-1 do begin
           singlearc=deimos_read_chip(usenames[j], chipno, $
                             sat = tmpsat, quick=quick)
           fullarc[*,*,j]=singlearc
           fullarcsat[*,*,j]=tmpsat
           fullarcivar[*, *, j] = ((deimos_badchip AND 1b) EQ 0b)/(singlearc > noisesq)
                 
 ;add arcs
         arcsat = arcsat*tmpsat
         pessimsat=pessimsat OR tmpsat
       endfor

; COMBINE arcs, throwing out saturated lines only in 1 arc
;       arcimage=combine_arcs(fullarc,fullarcsat,dirtyarc=dirtyimage)
       avearc=avsigclip(fullarc,fullarcivar)
       arcimage=avearc.flux
       arcivar=avearc.ivar
       delvarx,fullarc,fullarcsat,avearc

    endif else begin
      arcimage=deimos_read_chip(usenames[0], chipno, sat = tmpsat, quick=quick)
      arcsat=tmpsat
      arcivar= ((deimos_badchip AND 1b) EQ 0b)/(arcimage > noisesq)
      arcsat = arcsat*tmpsat
      pessimsat=arcsat
      pessimsat=pessimsat OR tmpsat
    endelse

;    arcivar = ((deimos_badchip AND 1b) EQ 0b)/(arcimage > noisesq)
    if n_elements(pixflatgood) gt 0 then arcivar = arcivar*pixflatgood

; interpolate for aesthetic reasons
    if quick le 0 then begin
        arcflatten = arcimage*deimos_pixflat_mul
        arcimage = djs_maskinterp(arcflatten, arcivar eq 0., iaxis=0, /const)
    endif    

; Find offset between arc and flat in spatial direction (should be small!)
; NB: we apply the offset only in the arc extraction.  It is the flat
;   positions that are written to the spCalib files. 

    dxarc1 = deimos_edge_offset(xpos1, xpos2, ypos1, ypos2, $
                                arcimage, arcivar,/SUM)


    vprint, 2,'Pass 1:  dxarc: ', dxarc1, ' [pix]'
; iterate once, just to be careful
    if quick le 0 then begin
       dxarc2 = deimos_edge_offset(xpos1+dxarc1, xpos2+dxarc1, $
                             ypos1, ypos2, arcimage, arcivar,/SUM)
       dxarc = dxarc1+dxarc2
       vprint,2, 'Pass 2:  dxarc: ', dxarc, ' [pix]'
    endif else dxarc = dxarc1  

;save,xpos1,xpos2,ypos1,ypos2,arcimage,arcivar,f='arcshift.sav'

; arc shifts: 21 sec  
 
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

;    arc_header = headfits(arcnames[0]) ;get header
    keylist = ['FRAMENO', 'OUTFILE', 'EXPTIME', 'DARKTIME', 'OBSERVER', $
         'OBJECT', 'OBSTYPE', 'ROTATVAL', 'DATE-OBS', 'UT', 'AIRMASS', $
         'TARGNAME', 'EPOCH', 'EQUINOX', 'DEC',  'RA', 'AZ', 'EL', $    
         'HA', 'ST', 'MJD-OBS', 'PARANG',  'SYNOPSIS', 'DWFILNAM', $
          'SLMSKNAM', 'GRATEPOS', 'HPLOGTIM']
     hdr = copy_keywords(arc_header, keylist, 'calibSlit')
     vers = spec2d_version() ;get version of code
     sxaddpar, hdr, 'COMMENT', 'DCS header cards taken from arc image'
     sxaddpar, hdr, 'SP2DVERS', vers, 'Version of spec2d' 
     sxaddpar, hdr, 'AUTHOR', 'Finkbeiner & Davis'
     sxaddpar, hdr, 'CHIPNO', chipno, 'chip number for this slitlet'


     gratepos = sxpar(arc_header, 'GRATEPOS')
; grab central wavelength from either slider 3 or 4 - already have - DPF
;     lambda_c = (gratepos eq 3) ? sxpar(arc_header, 'G3TLTWAV') : $
;                                  sxpar(arc_header, 'G4TLTWAV')
     sxaddpar, hdr, 'GRATING', grating,  'rule of grating'
     sxaddpar, hdr, 'GRTLTWAV', lambda_c,  'Nominal central wavelength'

;------------------------------------------------------------------------
     vprint, 2, 'Beginning loop over slitlets' 

     if chipisred then lamps = redlamps else lamps = bluelamps

     dline = fltarr(n_elements(lamps), nslitlets)
     dlinesig = fltarr(n_elements(lamps), nslitlets)

;         whereok = where(vigcolumn eq 0)

;         minok = min(whereok)
;         if minok ne 0 then minok = minok-15
;         maxok = max(whereok)
;         if maxok ne (sizexf-1) then maxok = maxok-15

         minok=0
         maxok=sizexf-1

         ;vprint,3,string( 'min/max columns used:', minok, maxok)
         vprint,3,''
         vprint,3,'-------------------------------------------------'
         vprint,3,''

   t1=systime(1)
  
    for slitno=0, nslitlets-1 do begin ;loop over all slitlets
       tstart = systime(1)

;retrieve true slitn, for use in output name and for possible selection
       bluslitno = slitcoords[indblu[slitno]].slitno
; get the slitpa and maskpa for this slit/mask.
       get_pa, bluslitno, slitpa=slitpa, maskpa=maskpa, $
         desitab=desitab, masktab=masktab
       nslitmatch = 0
       if keyword_set(slitlist) then $ ;check if on list of designated slitlets
       slityes = where(bluslitno eq slitlist, nslitmatch)
       slitwidth = slitcoords[indblu[slitno]].slitwidth ; mm size of slitlet
;do the analysis only if all objects to be reduced, or if on
;designated list, or if a wide slit.
    if (NOT keyword_set(onlyalign) and NOT keyword_set(slitlist)) OR $
          (keyword_set(onlyalign) and slitwidth gt 1.) OR $
          (keyword_set(slitlist) and $
              (nslitmatch gt 0 or (quick le 0 AND slitwidth gt 1.5))) $
          then begin

         x0 = xpos1[*, slitno] ; >0.
         x1 = xpos2[*, slitno] ;< (sizexf-1)

        ymid = n_elements(x0)/2

; testing 2/28/02 - keep more slits
        lowedge= (x0[ymid] < x1[ymid]) > 0
        highedge = (x0[ymid] > x1[ymid]) < (sizexf-1)

        slitisok=(highedge ge minok+25 AND lowedge LE maxok-25) $
          AND (highedge-lowedge) GE 15 AND (highedge-lowedge) LT 1024

        if (highedge-lowedge) LT 15 then $
          vprint,2,string( 'Bad slit - too narrow!! ',bluslitno)
         if (highedge-lowedge) GT 1024 then $
           vprint,2,string( 'Bad slit - too wide!! ',bluslitno)


        IF slitisok then begin 

           interp = (quick le 0)

;rectify arc 
         rect_arc = deimos_rectify_slit(arcimage,arcivar,x0+dxarc,x1+dxarc,$
            /interp,  xshift=xshift, npad=0, mask=arcmask) ;don't /recen

         rect_dirtyarc = deimos_rectify_slit(dirtyimage,arcivar,$
                                             x0+dxarc,x1+dxarc,$
            interp=interp,  xshift=xshift, npad=0, mask=arcmask) ;don't /recen

                                ; test to see if there's significant
                                ; light down this slit; if not, we'll
                                ; move on.  Do a cheesy CR avoidance first
      sizey = (size(rect_arc, /dimens))[1]
     
         testdata = rect_arc[100:sizeyf-100, sizey/2-2:sizey/2+2]
         testdata = testdata[sort(1E5-testdata)]
         
         nearmaxinarc = testdata[20]
         badslit = ((nearmaxinarc) lt 800)

; only proceed if badslit is not set
         if (badslit eq 0) then begin

            rect_arcsat=deimos_rectify_slit(float(arcsat),arcivar,x0+dxarc, $
                        x1+dxarc, interp=interp, xshift=xshift, npad=0) NE 0

            rect_arcivar = deimos_rectify_slit(arcivar, arcivar, x0+dxarc, $
                        x1+dxarc, interp=interp, xshift=xshift, npad=0)

; do not want to call lines saturated just because their middle
; crosses bad pixels
            rect_arcsat=rect_arcsat*(rect_arcivar NE 0.)

 
; -------- extract flat
            rect_flat = deimos_rectify_slit(flatinterp, flativar, x0, x1, $
                       interp=interp, npad=0,mask=rect_flatmask)
       
            rect_flatsat = deimos_rectify_slit(float(flatsat), flativar, $
                       x0, x1, interp=interp, npad=0) NE 0

            rect_nearvig = deimos_rectify_slit(float((nearvig AND 4b) eq 4b), nearvig, $
                          x0, x1, interp=interp, npad=0) NE 0

            rect_fullvig = deimos_rectify_slit(float((nearvig AND 2b) eq 2b), nearvig, $
                           x0, x1, interp=interp, npad=0) NE 0

            rect_nearvig = rect_fullvig*2b+rect_nearvig*4b
            delvarx,rect_fullvig
       
; before: used ivarmask=float(flativar eq 0)
; conserve memory instead of operations here
            rect_ivarmask = deimos_rectify_slit(float(flativar eq 0), float(flativar eq 0), x0, x1, $
                    interp=interp, npad=0)

            rect_flativar = deimos_rectify_slit(flativar, flativar, x0, x1, $
                      interp=interp, npad=0) * (rect_ivarmask eq 0)
	
	    rect_flativar=rect_flativar*(rect_flatmask)


            ncol = (size(rect_flat, /dimens))[1]
          
	
; -------- make the flat
        if slitwidth gt 1.95 then begin
	     vprint,2,'Setting aligment star flat to 1'
               flat2d = rect_flat*0.
               flat1d = fltarr(ncol)+1.
               flatmask = rect_ivarmask*0b+1b
               varslitfn=flat1d
               vigcorr = float(flatmask)
        endif else begin
           deimos_makeflat, rect_flat, flat2d, flat1d, vigcorr=vigcorr, $
              varslitfn=varslitfn, ivar=rect_flativar, $
              mask=flatmask, bitmask=rect_nearvig, quick=quick
       endelse
; if this is an alignment star, set slitfn & flat to 1...

     
; Note that we now use the 2d slit function...
            rect_arc = deimos_applyflat(rect_arc, flat2d, $
                        varslitfn, invvar=rect_arcivar, quick=quick)

            rect_dirtyarc = deimos_applyflat(rect_dirtyarc, flat2d, $
                        varslitfn, invvar=rect_arcivar, quick=quick)

; also, apply a vignetting correction to the arc

            whbadvig=where(vigcorr eq 0 OR finite(vigcorr) eq 0,vigbadct)
            if vigbadct gt 0 then vigcorr[whbadvig]=1.
            rect_arc=rect_arc/vigcorr
            rect_arcivar=rect_arcivar*vigcorr^2*flatmask*arcmask

         ENDIF  ; END SECTION CHECKED FOR BADSLIT

            if (badslit eq 0) then if total(finite(rect_arc) EQ 0b) NE 0 then badslit = 1
;--------- Fit for wavelengths across slitlet--------------

         polyx = 0.  ; define these
         s_coeff = 0.
	 nowave=0
         sigma=0.

         if (badslit eq 0) then begin 

; Note that arcsat now carries vignetting info:
            
;           if lambda_c gt 7000. AND grating ge 850 then begin
            if lambda_c gt 7000. AND grating ge 850 then begin   

                deimos_fitwavelength_withsky, rect_arc, rect_arcsat+rect_nearvig, $
                 rect_arcivar*(rect_nearvig EQ 0B), lamps, chipno, $
                 grating, lambda_c, slitcoords[indblu[slitno]], $
                 model_lambda[indblu[slitno]], $
                 wave, flat=flat, anamorph=anamorph, plot=plot,  $
                 dirtyarc=rect_dirtyarc, $
                 polyflag = polyflag,polyx=polyx, s_coeff=s_coeff, dlam=dlam, $
                 dline=dline1, dlsig=dlsig1, wset=wset,sigma=sigma
              endif else begin
                 deimos_fitwavelength_withsky, rect_arc, rect_arcsat, $
                 rect_arcivar, lamps, chipno, $
                 grating, lambda_c, slitcoords[indblu[slitno]], $
                 model_lambda[indblu[slitno]], $
                 wave, flat=flat, anamorph=anamorph, plot=plot,  $
                   dirtyarc=rect_dirtyarc, $
                 polyflag = polyflag,polyx=polyx, s_coeff=s_coeff, dlam=dlam, $
                 dline=dline1, dlsig=dlsig1, wset=wset,sigma=sigma
              endelse


            if n_elements(wave) EQ 1 OR n_elements(polyx) lt 2 $
		OR sigma gt 5.0 then begin 
               vprint,4,''
;               message, 'No wavelength solution, skipping', /info
               vprint,4, 'Poor wavelength solution, outputting diagnostics'
               save,rect_arcivar,rect_nearvig,lamps,chipno,grating,$
                 lambda_c,slitcoords,indblu,slitno,model_lambda,anamorph,$
                 rect_arc,rect_arcsat, rect_flat,rect_flatsat,rect_ivarmask,$
		 flat2d,flat1d,varslitfn,vigcorr, polyflag,flat,$
		 x0,x1,ymid,lowedge,highedge,flatmask,rect_flativar,$
                 file=strcompress('failure.'+string(bluslitno)+ $
                   	'.'+string(chipno)+'.sav',/REMOVE)
		if n_elements(wave) EQ 1 OR n_elements(polyx) lt 2 $
	           then nowave=1
               vprint,4,''
            endif 
;else begin
               
; -------- POLYFLAG structure format
               if polyflag then begin 

                  calibSlit = {flat: floatcompress(flat2d, ndig=12), $
                               ivar: floatcompress(rect_flativar, ndig=8), $
                               mask: byte(flatmask), slitfn: flat1d,  $
                               x0: x0, x1: x1, dlam: dlam, $
                               lambdax: polyx,  tiltx: s_coeff, $
                               rawflat:floatcompress(rect_flat),$
                               rawarc: floatcompress(rect_arc)}
                  
               endif else begin 
; -------- TRACESET structure format
                  ncol = (size(rect_arc, /dimens))[1]
                  dlam = fltarr(ncol)

                  calibSlit = {flat: floatcompress(flat2d, ndig=12), $
                               ivar: floatcompress(rect_flativar, ndig=8), $
                               mask: flatmask, slitfn: flat1d,  $
                               func: wset.func, xmin: wset.xmin, xmax: wset.xmax, $
                               coeff: wset.coeff, $
                               x0: x0, x1: x1, dlam: dlam, $
                               lambdax: polyx,  tiltx: s_coeff, $
                               rawflat:floatcompress(rect_flat),$
                               rawarc: floatcompress(rect_arc)}
               endelse
               dline[*, slitno]    = dline1
               dlinesig[*, slitno] = dlsig1

;            endelse


         endif 
         sxaddpar, hdr, 'SLITNO', bluslitno, 'slit number'
         sxaddpar, hdr, 'SLITX0', xpos1[sizeyf/2, slitno], 'center of trace'
         sxaddpar, hdr, 'SLITX1', xpos2[sizeyf/2, slitno], 'center of trace'
         sxaddpar, hdr, 'SLITWID', slitwidth, 'slit width in mm'
; write the maskpa and slitpa to the header if this is a DEEP2 mask.
         sxaddpar, hdr, 'MASKPA', maskpa[0], 'Mask PA on sky'
         sxaddpar, hdr, 'SLITPA', slitpa[0], 'Slit PA on sky'
         sxaddpar, hdr, 'X0SYNTH', synth1[slitno], $
           'Was X0 synthesized from other traces?'
         sxaddpar, hdr, 'X1SYNTH', synth2[slitno], $
           'Was X1 synthesized from other traces?'
         sxaddpar, hdr, 'WAVESIG', sigma, $
           'RMS error in lambda solution, in AA (1&10 are flags)'
         
         if polyflag then $
           sxaddpar, hdr,'WAVETYPE','POLYFLAG','wavelength solution method' $
         else sxaddpar,hdr,'WAVETYPE','TRACESET','wavelength solution method'
         
         
         
; use slit number from blueprint structure (slitcoords) for filename
         slitstr = string(bluslitno, format='(I3.3)')
         
         maskstr = STRCOMPRESS(maskname, /REMOVE_ALL)
; check if the mask is a DEEP2 mask. if so, trim to just the 4-digit
; mask number
         deimos_isdeep, isdeep, maskstr
        
         if chipno gt 4 then color = 'R' else $
           color = 'B'          ;red or blue side??

         if n_elements(badslit) eq 0 then badslit = 1
         
; if the slit is OK, write out the calibslit file
            IF (badslit EQ 0 AND nowave EQ 0) then begin
               fname = outdatadir+'calibSlit.'+maskstr+'.'+slitstr+ $
                 color+ '.fits'
               vprint, 3, 'Writing: ', fname
               
               mwrfits, calibSlit, fname, hdr, /create 
               
; Write the 2d slit function as HDU 2
               if quick le 0 then mwrfits, varslitfn, fname
               
               vprint, 3, string(systime(1)-tstart,  ' s -- slitlet ', bluslitno)
               vprint,3,'--------------------------------------------------'
               vprint,3,''
            ENDIF else begin
               vprint,3,''
               mstr = string('Slit', slitno, ' (blueprint #', bluslitno, $
                             ') had insufficient arc flux - Skipping!', $
                             format='(A,I4,A,I4,A)')
               message, mstr, /INFO
               vprint,3,mstr
               vprint,3,''
            ENDELSE

         endif else begin 
            print
            if n_elements(badslit) eq 0 then badslit = 1
            if (badslit ne 0) then vprint, 3, $
              string('Slit ', slitno, ' (blueprint # ', bluslitno, ') off edge of chip - Skipping!') else $
		vprint,3, $
              string('Slit ', slitno, ' (blueprint # ', bluslitno, ') had unsuccessful wavelength solution - Skipping!')
            vprint,3,''
         endelse 
  
        endif ;end of loop on selected slitlets
      endfor  ;end loop over slitlets

      dfname = 'arcwave_qa'+string(chipno, format='(I1)')+'.fits'

      if total(dline ne 0) ne 0 then begin
        writefits, dfname, dline
        mwrfits, dlinesig, dfname
        mwrfits, lamps, dfname
      endif


    endif else begin 
       vprint,2, 'Skipping chip!'
    endelse 
  endfor      ;end loop over chips

 
; verbosity levels

; 0 nothing
; 1  per mask
; 2  per chip
; 3  per slit
; 4  everything!
  return
end
