

;-----------------------------------------
;--------------- SETCOLORS ---------------
;-----------------------------------------
;-----------------------------------------
pro setcolors, NAMES=names, VALUES=ctindx, $
               NPLOTCOLORS=nplotcol, TABLESIZE=size, $
               DECOMPOSED=decomposed, $
               SYSTEM_VARIABLES=system_variables, $
               GETINFO=getinfo, $
               SILENT=silent, $
               TEST=test, $
               START=start, $
               NOCOLORS=nocolors,$
               PLOTCOLORS=plotcolors, $
               R_PLOT=r_plot, $
               G_PLOT=g_plot, $
               B_PLOT=b_plot, $
               PSEUDO256=pseudo256, $
               _REF_EXTRA=extra

  on_error, 0

; FIND IDL VERSION NUMBER... 5.2 OR HIGHER...
  if (float(!version.release) lt 5.2) then begin
      message,'This routine is only supported on IDL ' + $
        'version 5.2 or higher', /INFO
      return
  endif
  
; IF NOT IN X WINDOWS, DECOMPOSED DOESN'T MEAN ANYTHING!
  if (!d.name eq 'X') then begin
      
                                ; SET COLOR DECOMPOSITION EXPLICITLY?
      if (N_elements(DECOMPOSED) gt 0) then device, DECOMPOSED=decomposed
      
                                ; A WINDOW NEEDS TO HAVE BEEN CREATED TO ESTABLISH THE VISUAL TYPE...
      if (!d.window eq -1) then begin
          window, /FREE, /PIXMAP, COLORS=256*keyword_set(PSEUDO256)
          wdelete, !d.window
      endif
      
                                ; WHAT VISUAL STATES HAVE BEEN SET? 
      device, get_visual_name=visual, $
        get_visual_depth=depth, $
        get_decomposed=decomposed
      
                                ; CAN'T DISPLAY COLORS IF WE HAVE A GRAYSCALE VISUAL...
      if (strupcase(visual) eq 'GRAYSCALE') OR $
        (strupcase(visual) eq 'STATICGRAY') then begin
          message, 'Cannot set colors with '+strtrim(visual,2)+' visual!', /INFO
          if keyword_set(SILENT) $
            then return $
          else goto, output
      endif
      
  endif else decomposed=0
  
  size = !d.table_size          ; WHAT IS THE COLOR TABLE SIZE?
  ncol = !d.n_colors            ; HOW MANY COLORS ARE AVAILABLE?
  
; DO YOU JUST WANT THE CURRENT VISUAL INFO...
  if keyword_set(GETINFO) then goto, output
  
; THE IDL-DEFINED COLOR COMMON BLOCK...
  common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr
  
; LOAD THE ORIGINAL COLOR TABLE...
  if (N_elements(r_orig) ne 0) $
    then tvlct, r_orig, g_orig, b_orig $
  else loadct, 0, /silent
  
; DO WE USE PREDEFINED PLOTCOLORS OR USER-SPECIFIED...
  if (N_elements(PLOTCOLORS) eq 0) then begin
      
                                ; DO WE NOT WANT TO ADD LINE PLOT COLORS...
      if decomposed or (not decomposed AND not keyword_set(NOCOLORS)) then begin
          
                                ; NAMES OF THE PLOTCOLORS...
          names = [ 'black',  $
                    'red',    $
                    'orange', $
                    'green',  $
                    'forest', $
                    'yellow', $
                    'cyan',   $
                    'blue',   $
                    'magenta',$
                    'purple', $
                    'gray',   $
                    'white']
          
; GENERATE RGB COLORS FOR COLORED LINE PLOTS...
;                        b  r   o  g   f  y  c  b  m   p   g  w
          r_plot = 255b*[0, 1,  1, 0,  0, 1, 0, 0, 1,  0,  0, 1] + $
            byte([0, 0,  0, 0, 35, 0, 0, 0, 0,153,127, 0])
          g_plot = 255b*[0, 0,  0, 1,  0, 1, 1, 0, 0,  0,  0, 1] + $
            byte([0, 0,127, 0,142, 0, 0, 0, 0, 50,127, 0])
          b_plot = 255b*[0, 0,  0, 0,  0, 0, 1, 1, 1,  0,  0, 1] + $
            byte([0, 0,  0, 0, 35, 0, 0, 0, 0,205,127, 0])
          
      endif

  endif else begin

    ; ONLY NEEDS TO BE DONE IF UNDECOMPOSED!!!!!
    ; COLORS PROTECTED OTHERWISE!!!

    ; MAKE SURE RIGHT # OF RGB & NAMES...
      if (N_elements(R_PLOT)+N_elements(G_PLOT)+N_elements(B_PLOT))/3. $
        ne N_elements(PLOTCOLORS) then $
        message, 'PLOTCOLORS, R_PLOT, G_PLOT, & B_PLOT must have same size!'

    ; PLOTCOLORS HAD BETTER BE A STRING ARRAY...
      if (size(PLOTCOLORS,/TNAME) ne 'STRING') then $
        message, 'PLOTCOLORS must be of type STRING!'
      
    ; IF YOU'RE SILLY AND PUT MORE THAN ONE VALUE OF A COLOR IN
    ; THE TABLE... OH WELL, YOUR PROBLEM...
      if not keyword_set(START) then begin

        ; FORCE COLOR 'WHITE' TO BE AT TOP OF COLOR TABLE...
          good = where(strupcase(plotcolors) ne 'WHITE', ngood)
          if (ngood ne 0) then begin
              plotcolors = plotcolors[good]
              names  = [plotcolors,'white']
              r_plot = [r_plot[good],size-1]
              g_plot = [g_plot[good],size-1]
              b_plot = [b_plot[good],size-1]
          endif else names = plotcolors
          
      endif else begin

        ; FORCE COLOR 'BLACK' TO BE AT BOTTOM OF COLOR TABLE...
          good = where(strupcase(plotcolors) ne 'BLACK', ngood)
          if (ngood ne 0) then begin
              plotcolors = plotcolors[good]
              names  = ['black',plotcolors]
              r_plot = [234,r_plot[good]]
              g_plot = [123,g_plot[good]]
              b_plot = [211,b_plot[good]]
          endif else names = plotcolors
          
      endelse
  endelse

; HOW MANY PLOTCOLORS HAVE WE CREATED?
  nplotcol = N_elements(names)
  
; DO YOU EVEN HAVE THIS MANY COLORS AVAILABLE???
; AS CRAZY AS IT SOUNDS, THIS HAS HAPPENED!
  if (ncol le nplotcol) then begin
      message, 'Close some windows! Not enough available colors!', /info
      return
  endif
  
  if not decomposed then begin 
    ; WE HAVE AN 8-BIT DISPLAY, OR 24-BIT WITH COLOR DECOMPOSITION OFF
      
    ; STRETCH THE CURRENT COLOR TABLE OVER THE RANGE THAT 
    ; WON'T BE USED BY OUR PLOT COLORS...
    ; THE CURRENT COLOR TABLE IS STORED IN "ORIG" COLOR VECTORS.
    ; AFTER THE STRETCH, THE "CURR" COLOR VECTORS WILL CONTAIN
    ; THE NEW COLOR TABLE...
      if not keyword_set(START) then begin
          start = size-nplotcol
          stretch, 0, start-1
      endif else begin 
          start = 0
          stretch, nplotcol, size-1 
      endelse

    ; IF WE WANT PLOT COLORS, ADD THEM TO COLOR TABLE...
      if not keyword_set(NOCOLORS) then begin
        ; THE COLOR TABLE INDICES OF THE LINE PLOT COLORS...
          ctindx = bindgen(nplotcol)+byte(start)

        ; STORE THEM IN THE SPECIFIED ROWS OF THE NEW COLOR TABLE...
          r_curr[ctindx] = r_plot
          g_curr[ctindx] = g_plot 
          b_curr[ctindx] = b_plot
        
        ; LOAD THE NEW COLOR TABLE...
          tvlct, r_curr, g_curr, b_curr
      endif

    ; WE WANT TO LEAVE THE "ORIG" VALUES AS THEY WERE WHEN SENT IN!
    ; IF THEY WERE SET TO "CURR" VALUES NOW, IF THIS ROUTINE WERE CALLED
    ; AGAIN, THE PLOTCOLORS WOULD BE STRETCHED INTO THE COLOR TABLE!

  endif $
; OTHERWISE WE HAVE A 24-BIT DISPLAY...
  else ctindx = long(r_plot) + ishft(long(g_plot),8) + ishft(long(b_plot),16)

; SET THESE COLOR NAMES AS SYSTEM VARIABLES...
  if keyword_set(SYSTEM_VARIABLES) then $
    for i = 0, N_elements(names)-1 do begin
      success = execute('defsysv, "!'+names[i]+'", long(ctindx[i]), 0')
      if not(success) $
        then message, '!'+strupcase(names[i])+' not changed.', /cont
  endfor
  
; DISPLAY THE COLOR TABLE AND SEE IF THE COLORS ARE CORRECT...
  if keyword_set(TEST) then begin
      
      if (!d.name eq 'X') then $
        window, 28, ysize=2*size, xsize=size, TITLE='IDL 28: Plot Colors'
      if (!d.name eq 'PS') then $
        device, xsize=3, ysize=6, xoff=(8.5-3)/2., yoff=(11.-6)/2., /inches
      device, set_font='Helvetica Bold', /tt_font
      tv, [indgen(size),indgen(size)] ## (intarr(size)+1)
      if keyword_set(SYSTEM_VARIABLES) $
        then prefix = '!!' $
      else prefix = ''
      for i = 0.0, nplotcol-1.0 do begin
          printcolor = 'xyouts, 0.45, '+strtrim((i+1)/(nplotcol+1),2)$
            +', '''+prefix+names[i]+''', font=1, /norm, charsize=2.0, co='$
            +strtrim(long(ctindx[i]),2)
          useless = execute(printcolor)
      endfor
  endif
  
; OUTPUT THE COLOR TABLE INFO...
output:
  if not keyword_set(SILENT) or keyword_set(GETINFO) then begin
      print
      print, 'Display Device  : ', !d.name
      if (!d.name eq 'X') then begin
          print, 'Visual Class    : ', visual
          print, 'Visual Depth    : ', strtrim(depth,2)+'-Bit'
      endif
      print, 'Color Table Size: ', strtrim(!d.table_size,2)
      print, 'Number of Colors: ', strtrim(!d.n_colors,2)
      print, 'Decomposed Color: ', strtrim(decomposed,2)
      print
  endif
  
end ; setcolors

;-----------------------------------------
;--------------- GETCOLOR ----------------
;-----------------------------------------
;-----------------------------------------
function getcolor, colorname, colornames, tableindices, $
                   DECOMPOSED=decomposed

; A LOT OF ERROR CATCHING INVOLVED!
  on_error, 2

  case 1 of 
; FIND IDL VERSION NUMBER... 5.2 OR HIGHER...
      float(!version.release) lt 5.2 : $
        message,'This routine is only supported on ' + $
        'IDL version 5.2 or higher', /INFO
      
; ONLY ACCEPTS ONE COLOR NAME...
      size(colorname, /n_dimensions) ne 0 : $
        message, 'COLORNAME cannot be an array!', /INFO

; EITHER SEND IN BOTH COLORNAMES AND TABLE INDICES WITH THE COLOR NAME
; OR SEND IN JUST THE COLOR NAME...
      n_params() eq 2 : begin
          message, 'Syntax:', /INFO
          message, ' color = getcolor(colorname, names, tableindx)', /INFO
          message, '   OR', /INFO
          message, ' color = getcolor(colorname [, DECOMPOSED=decomposed])', /INFO 
      end

; THE REQUESTED COLOR MUST BE A STRING...
      size(colorname, /tname) ne 'STRING' : $
        message, 'COLORNAME must be a string!', /INFO
      
; THE COLORNAMES MUST BE DEFINED...
      N_params() eq 3 AND N_elements(colornames) eq 0 : $
        message, 'COLORNAMES must be defined!', /INFO
      
; TABLEINDICES MUST BE DEFINED...
      N_params() eq 3 AND N_elements(tableindices) eq 0 : $
      message, 'TABLEINDICES must be defined!', /INFO
      
; THE COLORNAMES MUST BE STRINGS...
      N_params() eq 3 AND size(colornames, /tname) ne 'STRING' : $
        message, 'COLORNAMES must be strings!', /INFO
      
; IF ONLY A COLORNAME IS SENT...
      N_params() eq 1 AND N_elements(colorname) ne 0 : $
        begin
          
; CAN'T USE DECOMPOSED IF POSTSCRIPT!!!
; do we just set dec=0 if not x-win?
          if (!d.name ne 'X') then decomposed=0 $
          else $
; IF DECOMPOSED NOT EXPLICITLY SET, USE CURRENT VALUE...
          if (N_elements(DECOMPOSED) eq 0) $
            then device, get_decomposed=decomposed
          
; IS SETCOLORS ON THE IDL PATH...
;          if float(!version.release) ge 5.4 then begin
;              found = file_which('setcolors.pro')
;              if (found eq '') then begin
;                  message, 'SETCOLORS.PRO not found on IDL path!', /INFO
;                  goto, error
;              endif
;          endif
          
; RUN SETCOLORS TO ESTABLISH THE LINE PLOT COLORS...
          setcolors, NAMES=colornames, VALUES=tableindices, $
            DECOMPOSED=decomposed, /SILENT
          
          goto, noerror
      end
      
; NO PROBLEMS HERE...
      else : goto, noerror
  endcase
  
error:
; IF THERE WAS AN ERROR, SPLIT...
  message, 'Using system plot color: !p.color = '+strtrim(!p.color,2), /INFO
  return, !p.color
  
noerror:
; IS THE REQUESTED COLOR DEFINED?
  colorname = strupcase(strtrim(colorname,2))
  colorindx = where(strupcase(colornames) eq colorname, match)
  if (match eq 0) then begin
      message, 'Color '+colorname+' not found in list of colors:', /INFO
      message, strjoin(colornames,' '), /INFO
      goto, error
  endif
  
; RETURN THE COLOR TABLE INDEX OF THIS COLOR...
  return, (tableindices[colorindx])[0]
end

;-----------------------------------------
;------------ MAIN ROUTINE ---------------
;------------   qa_check   ---------------
;-----------------------------------------

pro qa_check, result=result, doplot=doplot, silent=silent

; find the object info file.
  objfile = findfile('obj_info*.fits*', count=nfiles)
  if nfiles gt 0 then objfile = objfile[0] $
  else message, 'ERROR: No obj_info.xxxx.fits file found!'
; read in the object info file.
  objinfo = mrdfits(objfile, 1, objhdr, /silent)

; find the bin table file.
  binfile = findfile('*.bintabs.*fits*', count=nfiles)
  if nfiles eq 0 then $
    message, 'ERROR: No xxxx.bintabs.fits file found!'
  binfile = binfile[0]
  fits_open, binfile, fcb
  extnames = strcompress(fcb.extname, /rem)
  fits_close, fcb

; read in the ObjectCat table from the bin table file.
  catdex = where(extnames eq 'ObjectCat', catcnt)
  if catcnt eq 0 then $
    message, 'ERROR: No ObjectCat table found!'
  objcat = mrdfits(binfile, catdex[0], /silent)
; determine the mask number and date observed.
  maskname = sxpar(objhdr, 'SLMSKNAM')
  date = strcompress(sxpar(objhdr, 'DATE-OBS'), /rem)
  airmass = sxpar(objhdr, 'AIRMASS')
; if this is a DEEP2 mask, then trim the maskname.
  deimos_isdeep, isdeep, maskname
; extract the positions of all of the slits for this mask.
  bludex = where(extnames eq 'BluSlits', blucnt)
  if blucnt eq 0 then $
    message, 'ERROR: No BluSlits table found!'
  desidex = where(extnames eq 'DesiSlits', desicnt)
  if desicnt eq 0 then $
    message, 'ERROR: No BluSlits table found!'
  bluslits = mrdfits(binfile, bludex[0], /silent)
  desislits = mrdfits(binfile, desidex[0], /silent)
; calculate the x and y positions in units of mm.
  xpos = (bluslits.slitx1 + bluslits.slitx2 + $
          bluslits.slitx3 + bluslits.slitx4)/4.
  ypos = (bluslits.slity1 + bluslits.slity2 + $
          bluslits.slity3 + bluslits.slity4)/4
; match the bluslits structure to the desislits via the dslitid tag.
  xpos = xpos[sort(bluslits.dslitid)]
  ypos = ypos[sort(bluslits.dslitid)]
  slitname = desislits[sort(desislits.dslitid)].slitname

; read in the zcat, if it is there.
  if isdeep then begin
      d2dir = getenv('D2_RESULTS')
      deep2products = getenv('DEEP2PRODUCTS')
      zfile = concat_dir(deep2products, 'zcat.*.fits*')
      zfile = findfile(zfile, count=zcnt)
      if zcnt gt 0 then begin
          zfile = zfile[reverse(sort(zfile))]
          print, 'Reading from zcat: ' + zfile[0]
          zz = mrdfits(zfile[0], 1, /silent)
          dex = where(strcompress(zz.maskname, /rem) eq maskname and $
                      strcompress(zz.date, /rem) eq date, cnt)
          if cnt gt 0 then begin
              zz = zz[dex]
          endif else print, '(qa_check.pro) No zresult info found for this mask!'
      endif else print, '(qa_check.pro) No zcat file found!'
  endif

; the obj_info file contains the s/n numbers for each object as
; derived from their spatial profiles (recall that the obj_info file
; is created by the routine do_extract.pro). let's also find the s/n
; for each object in a more empirical manner by determining the mean
; and standard deviation of the 1-d spectra in a given lambda
; range. this second means for estimating s/n is slow because the
; spec1d files must be read.

; construct a list of all of the spec1d file names using the
; information stored in the objinfo file.
  files = 'spec1d.' + maskname + '.' + $
    string(objinfo.slitno, format='(I3.3)') + '.' + $
    strcompress(objinfo.objno, /rem) + '.fits'
  nfiles = n_elements(files)
; set a threshold value determining the number of pixels needed in the
; wavelength region (given by blam) in order to do statistical
; test. if there isn't enough pixels, then the s/n for that slitlet
; will simply be taken to be 0.
  thresh = 500. ;pixels
; define the wavelength regions in which we will measure s/n. these
; regions are chosen to avoid sky lines but still provide a large
; enough number of pixels (6600 < lambda < 6830 and 6960 < lambda <
; 7250).
  blam = [6600., 6830., 6960., 7250.] ;angstroms 

; define the output array of structures and fill with the flag value.
  flagval = -999.
; s/n from mean/sigma of 1-d spectra: sn_spec1d.
; s/n from spatial profile (over full fwhm): sn_sprof.
; s/n from spatial profile (over 5 pixels): sn_detect.
  tmp = {mask:maskname, objno:'', sn_spec1d:flagval, $
         sn_sprofB:flagval, sn_sprofR:flagval, $
         sn_detectB:flagval, sn_detectR:flagval, $
         mag:flagval, xpos:flagval, ypos:flagval, $
         z:flagval, zquality:-1, edge_dist:flagval, $
         slitlength:flagval, astar_seeing:flagval}
  ss = replicate(tmp, nfiles)

; loop through the spec1d files and determine the various s/n values.
  for i=0,nfiles-1 do begin

; --------------- determine sn_spec1d --------------
      ss1d = fill_gap(files[i], /silent)
; check that fill_gap worked.
      if size(ss1d, /tname) eq 'STRUCT' then begin
; extract the portion of the spectrum in the selected lambda range.
          region = where((ss1d.lambda gt blam[0] and $
                          ss1d.lambda lt blam[1]) or $
                         (ss1d.lambda gt blam[2] and $
                          ss1d.lambda lt blam[3]), cnt)
          if cnt gt thresh then begin
; calculate the signal-to-noise.
              djs_iterstat, ss1d.spec[region], mean=signal, $
                sigma=noise, sigrej=3.0
              ss[i].sn_spec1d = signal / noise
          endif
      endif else $ 
        if not(keyword_set(silent)) then $
        print, 'ERROR: fill_gap.pro error for file ' + files[i]

; ---------- determine sn_sprof & sn_detect ----------
; extract the object number from the spec1d header.
      hdr = headfits(files[i], ext=1, /silent)
      objno = sxpar(hdr, 'OBJNO')
      objno = strcompress(objno, /rem)
      ss[i].objno = objno
; also extract the slit number.
      slitno = sxpar(hdr, 'SLITNO')
      slitno = string(slitno, format='(I3.3)')
; find the entries in the objinfo file and the ObjectCat table for
; this object.
      dex = where(strcompress(objinfo.objno, /rem) eq objno, cnt)
      if cnt gt 0 then begin
          for j=0,cnt-1 do begin
              dd = dex[j]
              hue = objinfo[dd].color
              if hue eq 'B' then begin
                  ss[i].sn_sprofB = objinfo[dd].s2n_fwhm
                  ss[i].sn_detectB = objinfo[dd].s2n_window
                  ss[i].slitlength = objinfo[dd].nrows
              endif
              if hue eq 'R' then begin
                  ss[i].sn_sprofR = objinfo[dd].s2n_fwhm
                  ss[i].sn_detectR = objinfo[dd].s2n_window
                  ss[i].slitlength = objinfo[dd].nrows
              endif
              ss[i].edge_dist = $
                (objinfo[dd].nrows - objinfo[dd].cat_objpos) < $
                (objinfo[dd].cat_objpos)
          endfor
      endif

      dex = where(strcompress(objcat.object, /rem) eq objno, cnt)
      if cnt gt 0 then $
        ss[i].mag = objcat[dex[0]].mag

      dex = where(string(slitname, format='(I3.3)') eq slitno, cnt)
      if cnt gt 0 then begin
          ss[i].xpos = xpos[dex[0]]
          ss[i].ypos = ypos[dex[0]]
      endif

      if keyword_set(zz) then begin
          dex = where(string(zz.slitname, format='(I3.3)') eq slitno, cnt)
          if cnt gt 0 then begin
              dd = dex[0]
              ss[i].z = zz[dd].z
              ss[i].zquality = zz[dd].zquality
          endif
      endif
  endfor

; determine the z-completeness of the mask.
  zdex = where(ss.zquality ge 3, zgood)
  zgood = zgood / n_elements(ss)

; get the seeing value as derived from the alignment stars
  astar = where(objinfo.objtype eq 'A', acnt)
  if acnt eq 0 then begin
      print, '(qa_check.pro) No alignment stars found!!'
  endif else begin
      pixscl = 0.117371         ;arcsec/deimos pixel
      astar_seeing = median(objinfo[astar].fwhm) * pixscl
  endelse

; also get the seeing from the distribution of fwhm's.
  cfht_pixscl = 0.207 ;arcsec/cfht pixel
  psee = sxpar(objhdr, 'PCAT_SEE') * cfht_pixscl ;in arcseconds
  seediff = sxpar(objhdr, 'SEEDIFFB') ;in arcseconds
  if seediff lt 0 then pcat_seeing = sqrt(psee^2 - seediff^2) $
  else pcat_seeing = sqrt(psee^2 + seediff^2)


; define a set of good objects which we will use to characterize the
; overall properties of the mask.
  good = where(ss.sn_spec1d ne flagval and finite(ss.sn_spec1d) eq 1 and $
               ss.sn_sprofB ne flagval and finite(ss.sn_sprofB) eq 1 and $
               ss.sn_sprofR ne flagval and finite(ss.sn_sprofR) eq 1 and $
               ss.sn_detectB ne flagval and finite(ss.sn_detectB) eq 1 and $
               ss.sn_detectR ne flagval and finite(ss.sn_detectR) eq 1 and $
               ss.mag ne flagval, gnum)

; only proceed if there are enough "good" objects.
  if gnum lt 50 then begin
      print, '(qa_check.pro) Not enough points to continue!'
      result = flagval
      return
  endif
  gg = ss[good]

; make sure that qa directory exists.
  d2dir = getenv('D2_RESULTS')
  qadir = concat_dir(d2dir, 'qa')
  spawn, 'mkdir -p ' + qadir, spres, sperr

;----------------------------------
; plot the results in a single post-script file.
  if keyword_set(doplot) and isdeep then begin
      set_plot, 'ps'
      qapsdir = concat_dir(d2dir, 'qa/ps')
      spawn, 'mkdir -p ' + qapsdir, spres, sperr
      plotfile = 'qa/ps/qa.' + maskname + '.' + date + '.ps'
      plotfile = concat_dir(d2dir, plotfile)
      device, file=plotfile, /landscape, /color
      tname = 'mask ' + maskname

; PLOT #1
; plot histogram of sn_detect distribution. recall that sn_detect is
; the s/n over a 5 pixel window in the spatial profile.
      xname = 'detection signal-to-noise'
      subname = 'signal-to-noise from spatial profile'
      binsz = 5.
      plothist, gg.sn_detectB, xhistB, yhistB, bin=binsz, /noplot
      plothist, gg.sn_detectR, xhistR, yhistR, bin=binsz, /noplot
      ymax = max([yhistB, yhistR])
      xmax = max([xhistB, xhistR])
      plothist, gg.sn_detectB, bin=binsz, /xsty, xr=[0, 200], $
        xtitle=xname, thick=2, xthick=1.5, ythick=1.5, /ysty, $
        yr=[0, 1.1*ymax], title=tname, subtitle=subname
      plothist, gg.sn_detectR, bin=binsz, line=2, thick=2, /overplot
      legend, ['Solid Line = Blue side of slit', $
               'Dashed Line = Red side of slit'], /right_legend
; PLOT #2
; plot histogram of s/n (fwhm) distribution...for both the blue and
; red portions of slit.
      xname = 'signal-to-noise (per pixel)'
      subname = 'signal-to-noise from spatial profile'
      binsz = 0.25
; catch for places where the signal-to-noise might be a non-finite
; value...just in case.
      plothist, gg.sn_sprofB, xhistB, yhistB, bin=binsz, /noplot
      plothist, gg.sn_sprofR, xhistR, yhistR, bin=binsz, /noplot
      ymax = max([yhistB, yhistR])
      xmax = minmax([xhistB, xhistR])
      plothist, gg.sn_sprofB, bin=binsz, /xsty, $
        xr=[xmax[0], 20], xtitle=xname, $
        thick=2, xthick=1.5, ythick=1.5, /ysty, yr=[0, 1.1*ymax], $
        title=tname, subtitle=subname
      plothist, gg.sn_sprofR, bin=binsz, line=2, thick=2, /overplot
      legend, ['Solid Line = Blue side of slit', $
               'Dashed Line = Red side of slit'], /right_legend

; PLOT #3
; plot histogram of s/n (empirical) distribution...calculated on the
; blue side of slit (6600A - 7300A).
      xname = 'signal-to-noise (per pixel)'
      subname = 'signal-to-noise from 1-d spectra'
      binsz = 0.5
      plothist, gg.sn_spec1d, bin=binsz, xtitle=xname, thick=2, $
        xthick=1.5, ythick=1.5, title=tname, subtitle=subname, $
        /xsty, xr=[-0.5, 10]

; PLOT #4
; plot s/n (from spec1d) versus s/n (from spatial profile).
      plot, gg.sn_spec1d, gg.sn_sprofB, psym=2, xthick=1.5, $
        ythick=1.5, xtitle='s/n from 1-d spectra (per 1d pixel)', $
        ytitle='s/n from spatial profile (per 1d pixel)', $
        /xsty, /ysty, xr=[-0.5, 2.5], yr=[0.0, 3.0], title=tname
      oplot, findgen(10)-5, findgen(10)-5, line=2, thick=2


; PLOT #5
; plot s/n (sn_detect) versus magnitude.
      xname = strcompress(objcat[0].pband, /rem) + ' magnitude'
      yname = 'signal-to-noise (of detection)'
      plot, gg.mag, gg.sn_detectB, psym=6, xtitle=xname, $
        ytitle=yname, /xsty, /ysty, xr=[20, 24.6], $
        yr=[0, 120], xthick=1.5, ythick=1.5, title=tname
      oplot, gg.mag, gg.sn_detectR, psym=4
      legend, ['Square = Blue side of slit', $
               'Diamond = Red side of slit'], /left_legend

; PLOT #6
; plot s/n (from spatial profile) versus magnitude.
      yname = 'signal-to-noise (per pixel)'
      plot, gg.mag, gg.sn_sprofB, psym=6, xtitle=xname, $
        ytitle=yname, /xsty, /ysty, xr=[20, 24.6], $
        yr=[-0.25, 6], xthick=1.5, ythick=1.5, title=tname
      oplot, gg.mag, gg.sn_sprofR, PSYM=4
      legend, ['Square = Blue side of slit', $
               'Diamond = Red side of slit'], /left_legend

; PLOT #7
; plot s/n (from spec1d files) verus magnitude.
      yname = 'signal-to-noise per pixel (empirical)'
      plot, gg.mag, gg.sn_spec1d, psym=5, xtitle=xname, $
        ytitle=yname, /xsty, /ysty, xr=[20, 24.6], $
        yr=[-0.25, 5], xthick=1.5, ythick=1.5, title=tname

  endif

;---------------------------
; now fit the data s/n versus magnitude data with a median-median
; fit. use the s/n data from the blue side of the slit since it is
; closer to the R band magnitudes in the objcat.
; bin the magnitude data into 2 bins.
  dex1 = where(gg.mag gt 23.1 and gg.mag le 23.6, cnt1)
  dex2 = where(gg.mag gt 23.6 and gg.mag le 24.1, cnt2)
  if cnt1 gt 0 and cnt2 gt 0 then begin
; take median values in both X and Y directions.
      x1 = median(gg[dex1].mag)
      x2 = median(gg[dex2].mag)
      y1 = median(gg[dex1].sn_sprofB)
      y2 = median(gg[dex2].sn_sprofB)
      z1 = median(gg[dex1].sn_spec1d)
      z2 = median(gg[dex2].sn_spec1d)
; fit a line to the median-median points (for sn_sprof).
      mcc_polyfit, [x1, x2], [y1, y2], [0, 1], a=a
      xfit = gg[sort(gg.mag)].mag
      yfit = a[0] + xfit*a[1]
; do the same for the sn_spec1d data.
      mcc_polyfit, [x1, x2], [z1, z2], [0, 1], a=b
      zfit = b[0] + xfit*b[1]
; evaluate the linear fit at R=23.5.
      result = [a[0] + (23.5)*a[1], $
                b[0] + (23.5)*b[1]]

; now determine distance of point from the fit line (at that
; magnitude). that is, find distance (vertically) in s/n.
      if keyword_set(doplot) and isdeep then begin
          subset = where(gg.mag gt 23., subcnt)
          if subcnt lt 5 then $
            message, 'ERROR not enough points to continue!'
          fitvals = a[0] + gg[subset].mag*a[1]
          disty = gg[subset].sn_sprofB - fitvals
; bin the distance data.
          binsz = 0.25
          plothist, disty, xhist, yhist, bin=binsz, /noplot
          nbins = n_elements(xhist)

; PLOT #8 
; plot s/n vs. mag for all objects showing the median-median fit.
          yname = 'signal-to-noise (per pixel)'
          plot, gg.mag, gg.sn_sprofB, psym=6, xtitle=xname, $
            ytitle=yname, /xsty, /ysty, xr=[20., 24.6], $
            yr=[-0.25, 6], xthick=1.5, ythick=1.5, title=tname
          oplot, xfit, yfit, thick=2
          legend, ['Square = Blue side of slit', $
                   's/n (R=23.5) = ' + $
                   string(result[0], format='(F6.2)')], $
            /left_legend

; PLOT #9
; plot s/n vs. mag for all object with the size of the symbol scaled
; according to distance from the fit line. 
          plot, xfit, yfit, thick=2, xtitle=xname, ytitle=yname, $
            /xsty, /ysty, xr=[22.8, 24.2], yr=[-0.25, 4.5], xthick=1.5, $
            ythick=1.5, title=tname
          maxsz = 5.
; make a user-defined plotting symbol: an open circle.
          res = 32
          foo = findgen(res) * (!pi * 2 / (res-1))
          usersym, cos(foo), sin(foo)
          for i=0,nbins-1 do begin
              dex = where(disty lt xhist[i] and $
                          disty ge xhist[i]-binsz, cnt)
              if cnt gt 0 then $
                oplot, (gg[subset].mag)[dex], $
                (gg[subset].sn_sprofB)[dex], psym=8, $
                symsize=(i*maxsz/nbins+0.5)
          endfor
          
; PLOT #9
; plot of deimos x position vs. deimos y position for all the
; objects....use the same symbol size coding.
          xname = 'Deimos x position'
          yname = 'Deimos y position'
          maxsz = 5.
          pflag = 0
          for i=0,nbins-1 do begin
              dex = where(disty lt xhist[i] and $
                          disty ge xhist[i]-binsz, cnt)
              if cnt gt 0 then begin
                  if i ge nbins/4. and i lt nbins/3.*2. then hue = getcolor('red')
                  if i ge nbins/3.*2. then hue = getcolor('green')
                  if i lt nbins/4. or pflag eq 0 then hue = getcolor('black')
                  if pflag eq 0 then begin
                      plot, (gg[subset].xpos)[dex], (gg[subset].ypos)[dex], $
                        xtitle=xname, ytitle=yname, $
                        /xsty, /ysty, xr=[-400, 400], yr=[0, 200], $
                        xthick=1.5, ythick=1.5, title=tname, $
                        symsize=(i*maxsz/nbins+1), psym=8, color=hue
                      pflag = 1
                  endif else $
                    oplot, (gg[subset].xpos)[dex], (gg[subset].ypos)[dex], $
                    psym=8, symsize=(i*maxsz/nbins+0.5), color=hue 
              endif
          endfor
          
; PLOT #10
; plot z quality versus s/n (sn_sprof).
;          xm = minmax(gg.sn_sprofB)
          xm = [0, 10]
          xname = 'signal-to-noise (per pixel)'
          yname = 'z quality'
          plot, gg.sn_sprofB, gg.zquality, psym=4, xtitle=xname, $
            ytitle=yname, /ysty, /xsty, yr=[-0.5, 4.5], $
            xr=[xm[0], xm[1]], title=tname, xthick=1.5, $
            ythick=1.5
          
; PLOT #11
; plot z quality versus DEIMOS x,y.
          xname = 'DEIMOS x position'
          yname = 'DEIMOS y position'
          plot, [-1], [-1], psym=3, /xsty, /ysty, $
            xr=[-420, 400], yr=[0, 200], xtitle=xname, $
            ytitle=yname, title=tname, xthick=1.5, ythick=1.5
          for i=0,n_elements(gg)-1 do begin
              if gg[i].zquality eq -1 then $
                oplot, [gg[i].xpos], [gg[i].ypos], $
                psym=2, color=getcolor('black')
              if gg[i].zquality eq 0 then $
                oplot, [gg[i].xpos], [gg[i].ypos], $
                psym=6, color=getcolor('purple')
              if gg[i].zquality eq 1 then $
                oplot, [gg[i].xpos], [gg[i].ypos], $
                psym=4, color=getcolor('red')
              if gg[i].zquality eq 2 then $
                oplot, [gg[i].xpos], [gg[i].ypos], $
                psym=5, color=getcolor('blue')
              if gg[i].zquality eq 3 then $
                oplot, [gg[i].xpos], [gg[i].ypos], $
                psym=6, color=getcolor('green')
              if gg[i].zquality eq 4 then $
                oplot, [gg[i].xpos], [gg[i].ypos], $
                psym=1, color=getcolor('green')
          endfor
              oplot, [-410], [190], psym=2, color=getcolor('black')
              oplot, [-410], [180], psym=6, color=getcolor('purple')
              oplot, [-410], [170], psym=4, color=getcolor('red')
              oplot, [-410], [160], psym=5, color=getcolor('blue')
              oplot, [-410], [150], psym=6, color=getcolor('green')
              oplot, [-410], [140], psym=1, color=getcolor('green')
              xyouts, fltarr(6)-405, [188, 178, 168, 158, 148, 138], $
                ['Stars', 'Q=0', 'Junk', 'Q=2', 'Q=3', 'Q=4']
          
; PLOT #12
; plot zquality versus distance from edge of slit.
          xm = minmax(gg.edge_dist)
          xname = 'distance from edge of slit (pixels)'
          yname = 'z quality'
          plot, gg.edge_dist, gg.zquality, psym=4, xtitle=xname, $
            ytitle=yname, /ysty, /xsty, yr=[-0.5, 4.5], $
            xr=[xm[0], xm[1]], title=tname, xthick=1.5, $
            ythick=1.5
          

; PLOT #13
; plot zquality versus slit length.
          xm = minmax(gg.slitlength)
          xname = 'slit length (pixels)'
          yname = 'z quality'
          plot, gg.slitlength, gg.zquality, psym=4, xtitle=xname, $
            ytitle=yname, /ysty, /xsty, yr=[-0.5, 4.5], $
            xr=[xm[0]-5, xm[1]], title=tname, xthick=1.5, $
            ythick=1.5

; close the post-script file.
          device, /close
          set_plot, 'x'
      endif
  endif




; write the results of this routine to file.
  nrows = 1
  fxbhmake, outhdr, nrows, /date, /initialize
  sxaddpar, outhdr,  'maskname', maskname
  sxaddpar, outhdr,  'date-obs', date
  sxaddpar, outhdr, 'airmass', airmass
  sxaddpar, outhdr, 'asee', astar_seeing, 'seeing from alignment stars'
  sxaddpar, outhdr, 'psee', pcat_seeing, 'seeing from sprof'
  sxaddpar, outhdr, 'sn_val1', result[0], 's/n (sprof) at R=23.5'
  sxaddpar, outhdr, 'sn_val2', result[1], 's/n (spec1d) at R=23.5'
  sxaddpar, outhdr, 'zge3', zgood, 'z >= 3'
  outfile = 'qa/qa.' + maskname + '.' + date + '.fits'
  outfile = concat_dir(d2dir, outfile)
  print, 'Writing qa output file: ' + outfile
  mwrfits, ss, outfile, outhdr, /silent, /create



end


