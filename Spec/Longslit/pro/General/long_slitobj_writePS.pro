;+
; NAME:
;   long_slitobj_writePS
;
; PURPOSE:
;  Script for writing to Postscript the slits and extracted objects for
;  all reduced exposures within the local directory, similar to the image plotted
;  on ATV in LONG_LOOK, but now in easily printable form. Writes
;  out to files named 'slitobj*.ps'. 
;
;  For brevity we only plot the central 1/sqrt(2) in the y-direction
;
;  Adapted from LONG_LOOK.  
;
; CALLING SEQUENCE:
;   IDL> cd, 'Science' 
;   IDL> .run long_slitobj_writePS
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   
; PROCEDURES CALLED:
;   
; REVISION HISTORY:
;   21-Jul-2014  - Adapted from LONG_LOOK by K.G. Lee (MPIA)
;-  
;------------------------------------------------------------------------------

set_plot, 'ps'
!p.font=0

fils = findfile('sc*.fits*', count = nfil)
stdfil = findfile('std*.fits*', count = nstd)
IF nstd GT 0 THEN BEGIN 
    fils = [fils, stdfil]
    nfil = nfil + nstd
ENDIF
IF nfil GT 0 THEN path = '../'
IF nfil EQ 0 THEN BEGIN
    fils = findfile('Science/s*.fits*', count = nfil) 
    path = './'
ENDIF
IF nfil EQ 0 THEN message, 'Cannot find files'

; Loop over all science exposures
for ifil=0, nfil-1 do begin

   file_fits = fils[ifil]
   print, 'Plotting slits & objects for ', file_fits

   ;; Open PS device
   fstr_split = strsplit(file_fits, '-.',/extract)
   fname = 'slitobj-'+fstr_split[1]+'.ps'

   device, file=fname, /color, xsize=10, ysize= 0.9*10./sqrt(2.), $
           /inches,/encap 

   sciimg = xmrdfits(file_fits, 0, scihdr,/silent)
   sciivar= xmrdfits(file_fits,1,/silent)
   sky_model = xmrdfits(file_fits, 2,/silent)
   obj       = xmrdfits(file_fits, 3,/silent)
   outmask   = xmrdfits(file_fits, 4,/silent)
   objstruct = xmrdfits(file_fits, 5,/silent)
   dims = size(sciimg, /dim)
   nx = dims[0]
   ny = dims[1]

If TAG_EXIST(objstruct, 'FLX_SHFT_SPA') THEN $
  xshift = objstruct[0].FLX_SHFT_SPA $
ELSE xshift = 0.0

If xshift NE 0.0 AND NOT KEYWORD_SET(NOSHIFT) THEN BEGIN
    slitfile =    findfile(path + 'slits*.fits*', count = nslitfil)
    IF nslitfil GT 1 THEN BEGIN
        print, 'Multiple slit files found. Using ' ,slitfile[0] 
        slitfile = slitfile[0]
    ENDIF
    tset_slits = xmrdfits(slitfile[0], 1,/silent)
    tset_slits = long_shiftslits(tset_slits, xshift)
    slitmask = long_slits2mask(tset_slits)
    wavefile = findfile(path + 'wave*.fits*', count = nwavefil)
    IF nwavefil GT 1 THEN BEGIN
       print, 'Multiple wavelength files found. Using ', wavefile[0]
       wavefile = wavefile[0]
    ENDIF
    pixset  = xmrdfits(wavefile[0], 1,/silent)
    wavesvfile = repstr(wavefile[0], '.gz', '')
    wavesvfile =  repstr(wavesvfile, '.fits', '.sav')
    restore, wavesvfile
    piximg = long_wpix2image(pixset, tset_slits, XFIT = xfit $
                             , waveimg = waveimg)
ENDIF ELSE BEGIN
    slitfile =    findfile(path + 'slits*.fits*', count = nslitfil)
    IF nslitfil GT 1 THEN BEGIN
        print, 'Multiple slit files found. Using ' ,slitfile[0] 
        slitfile = slitfile[0]
    ENDIF
    IF nslitfil GT 0 THEN BEGIN
        slitmask = xmrdfits(slitfile[0], 0,/silent) 
        tset_slits = xmrdfits(slitfile[0], 1,/silent)
    ENDIF ELSE slitmask = fltarr(nx, ny) + 1.0D
    wavefile = findfile(path + 'wave*.fits*', count = nwavefil)
    IF nwavefil GT 1 THEN BEGIN
        print, 'Multiple wavelength files found. Using ', wavefile[0]
        wavefile = wavefile[0]
    ENDIF
    IF nwavefil GT 0 THEN waveimg = xmrdfits(wavefile[0], 0,/silent) $
    ELSE waveimg = fltarr(nx, ny)
ENDELSE

plot, findgen(nx), findgen(round(ny/sqrt(2.))), /nodata, $
       charsize=0.5, charthick=2, $
      xran=[150, 1950], xsty=1, ysty=1, pos=[0.04,0.04,0.98,0.93] 

xyouts, 0.45, 0.97, file_fits, charsize=0.8, charthick=2,/norm

IF KEYWORD_SET(tset_slits) THEN BEGIN
; ------
;; Expand slit set to get left and right edge
    traceset2xy, tset_slits[0], rows, left_edge
    traceset2xy, tset_slits[1], rows, right_edge
    ;oplot, left_edge, rows, psym = 3, color = djs_icolor('black')
    ;oplot, right_edge, rows, psym = 3, color = djs_icolor('black')

    ; Plot traced slits as shaded regions
    nobj = (size(rows))[2]
    loadct, 0,/silent
    for ii=0, nobj-1 do begin
       left_tmp = reform(left_edge[*,ii])
       right_tmp = reform(right_edge[*,ii])
       row_tmp = reform(rows[*,ii])
       ;; Trim to only the y-values on the plot axes
       ycut = where(row_tmp GE !y.crange[0] AND row_tmp LE !y.crange[1],npix)
       left_tmp = left_tmp[ycut]
       right_tmp = right_tmp[ycut]
       row_tmp = row_tmp[ycut]

       ;Greyscale is a random value within some range
       greycolor= round(96.*randomu(seed)) + 84L

       ; Downsample number of vertices to stop IDL complaining
       downsamp = 4L*lindgen(ceil(float(npix)/4.))
       downsamp = [downsamp, npix-1]

       polyfill, [left_tmp[downsamp], right_tmp[downsamp]], $
                 [row_tmp[downsamp], reverse(row_tmp[downsamp])],/data, $
                 color=greycolor

       ; Plot slit numbers along top axis
       ymax = max(row_tmp, i_ymax)
       xmid = (left_tmp[i_ymax]+right_tmp[i_ymax])/2.

       xyouts, xmid, !y.crange[1]+10, strtrim(ii+1,2), color=djs_icolor('blue'), $
               charsize=0.6, charthick=2

    endfor

 ENDIF

oplot, objstruct.xpos, objstruct.ypos, psym = 3, color=djs_icolor('red')


objid = strtrim(objstruct.slitid, 2) + '-' + strtrim(objstruct.objid, 2)
nobj = n_elements(objstruct)
FOR j = 0L, nobj-1L DO BEGIN
    xpix = djs_median(objstruct[j].XPOS) - 25.0
    ypix = (!y.crange[1]-!y.crange[0])*float(j)/float(nobj)+!y.crange[0]+20.
    ;print, xpix, ypix, objid[j]
    xyouts, xpix, ypix, strcompress(objid[j], /rem) $
               , color = djs_icolor('red'), charsize = 0.6, charthick=2
ENDFOR

device, /close

endfor 

set_plot, 'x'

END
