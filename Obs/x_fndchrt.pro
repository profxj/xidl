;+ 
; NAME:
; x_fndchrt   
;    Version 1.1
;
; PURPOSE:
;    Given an array of Obj name, RA, and DEC create a set of postscript
;  finding charts. 
;
; CALLING SEQUENCE:
;  x_fndchrt, targlist, OUTDIR=, IMSIZE=, SURVEY=, /ESO
;
; INPUTS:
;  targlist  -- ASCII file containing  (QSO,  RA,  DEC)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  imsize - Arcmin of image [default is 5']
;  circ  -- Radius of green to draw about the target [default: 5"]
;  /ESO -- Use the DSS at ESO
;  /radec -- Input is ['Name', 'RA:RA', 'DEC:DEC']
;  EPOCH=  -- Epoch of RA/DEC [default: 2000]
;  SKIP=   -- Skip lines at the start of the file
;  TWOCIRC = (x,y) offset in arcmin from field center
;  ADDCIRC = [x,y] offset in arcmin from field center for multiple circles
;  CONNECTX, CONNECTY -- Draw a line through these [x,y] coordinates
;  NOTES -- can add an extra string of text 
;
; OPTIONAL OUTPUTS:
;  OUTDIR=  -- Name of output directory
;
; COMMENTS:
;
; EXAMPLES:
;   x_fndchrt, 'targets.list'
;
; PROCEDURES/FUNCTIONS CALLED:
;  showfits
;  querydss
;  sdss_queryimage
;
; REVISION HISTORY:
;   21-Nov-2003 Written by JXP
;-
;------------------------------------------------------------------------------
pro x_fndchrt, targlist, OUTDIR=outdir, imsize=imsize, survey=survey, $
               FITS=fits, SHOWFITS=showfits, RADEC=radec, CIRC=circ, $
               EPOCH=epoch, DECI=deci, ESO=eso, SDSS=sdss, $
               CONNECTX=connectx, CONNECTY=connecty, COLOR=color, $
               TWOCIRC=twocirc, ADDCIRC=ADDCIRC, SKIP=skip, _EXTRA=extra, NOTES=NOTES

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'x_fndchrt, targlist, /ESO, EPOCH=, OUTDIR=, IMSIZE=, /SHOWFITS, /RADEC, /DECI, CIRC=, SURVEY=, SKIP= [v1.1]'
      return
   endif 

;  resolve_routine, 'sdss_queryimage'

  common xcommon_color, r_vector, g_vector, b_vector
  xicmm_colors

  if not keyword_set( OUTDIR ) then outdir='./'
  if not keyword_set( IMSIZE ) then imsize = 5.
  if not keyword_set( CIRC ) then circ = 5.
  if not keyword_set( SURVEY ) then survey = '2r'

  if keyword_set(SDSS) then begin
      ;; Check for wget
      spawn, 'which wget', blah
      if strmid(blah, 0, 1) NE '/' then begin
          print, 'x_fndchrt:  No wget on your machine. SDSS will not work'
          print, 'x_fndchrt:  Continue only if you are sure it is in your path.'
      endif
  endif

  ;; Readlist
  if not keyword_set( RADEC ) then readcol, targlist, nam, ra, dec, SKIPLINE=skip, FORMAT='A,A,A', COMME='#' $
    else begin
      if n_elements(targlist[*,0]) NE 3 then begin
          print, 'x_fndchrt: Dont forget the name!!'
          return
      endif
      nam = targlist[0,*]
      ra = targlist[1,*]
      dec = targlist[2,*]
  endelse

  ;; Loop
  nobj = n_elements(nam)

  for q=0L,nobj-1 do begin
     ;; Grab ra, dec
     if not keyword_set(DECI) then x_radec, ra[q], dec[q], rad, decd $
     else begin
        rad = double(ra[q])
        decd = double(dec[q])
        x_radec, ras, decs, rad, decd, /flip
        ra[q] = ras
        dec[q] = decs
     endelse
     

     ;; Precess if necessary
     if keyword_set(EPOCH) then begin
        precess, rad, decd, epoch, 2000.
        ;; For labeling of J2000
        x_radec, ras, decs, rad, decd, /flip
        ra[q] = ras
        dec[q] = decs
     endif
    
     
     ;; Grab dss image
      if not keyword_set(FITS) then begin
          if keyword_set(SDSS) then begin
              npix = round(imsize*60./0.39612)
;              scale = imsize*60. / isz[0]  ;; arcsec/pix
              if not keyword_set(COLOR) then gray = 1
              img = sdss_queryimage(RAD, DECD, xs=npix, ys=npix, GRAY=GRAY)
          endif else begin
              ;; DSS
              querydss, [rad,decd], img, hdr, imsize=imsize, survey=survey, ESO=eso 
          endelse
      endif else img = x_readimg(fits, head=hdr,/fscale)

;      querydss, [rad,decd], img, hdr, imsize=imsize, survey=survey 
      ;; Write to fits
      flg = 1
      if n_elements(img) GT 1 then begin
          if keyword_set(SVFITS) then imnm=nam[q]+'.fits' else imnm='tmp.fits'
          mwrfits, img, imnm, hdr, /create 
      endif else begin
          print, 'Image not found!  Try survey=''1'''
          flg=0
      endelse
      ;; Showfits
      if flg NE 1 OR keyword_set( NOPS ) then continue
      if keyword_set( SHOWFITS ) then begin
          spwncmd = 'showfits -objnm='+strtrim(nam[q],2)+' -ps -fi=' $
            +outdir+strtrim(nam[q],2)+'.ps '+imnm 
          spawn, spwncmd
      endif else begin
          ;; PSFILE
          psfil = outdir+strtrim(nam[q],2)+'.ps'
          x_psopen, psfil, /portrait
          state = { $
                    ncolors: 0L, $
                    brightness: 0.7, $
                    contrast: 0.5 $
                  }
          loadct, 0, /silent
          ncolors = !d.table_size - 9
          state.ncolors=ncolors
          r_vector = bytarr(ncolors)
          g_vector = bytarr(ncolors)
          b_vector = bytarr(ncolors)
          ximgd_getct, state, 0, /CLR
          ;; Invert
          r_vector = reverse(r_vector)
          g_vector = reverse(g_vector)
          b_vector = reverse(b_vector)
          ximgd_stretchct, state
          
          ;; Display image
          mx = max(img, min=mn)
          med = median(img)
          sig = stddev(img)
          pltmax = (med + (3.0*sig)) < mx
          pltmin = (med - (3.0*sig)) > mn
          print, pltmin, pltmax
          display_image = bytscl(img, min=pltmin, max=pltmax, /nan, $
                                 top=state.ncolors-1) + 8B

          ;; TV
          ntv = 200L
          ;ntv = 700L
          tv_image = congrid(display_image, ntv, ntv)
          tv, tv_image, 0.8, 1., xsize=6.0, ysize=6.0, /inches

          ;; Label
          clr = getcolor(/load)
          plot, [0], [0], /nodata, xtitle='Arc-minutes', ytitle='EAST',$ 
            /noerase, pos=[0.108,0.11,0.905,0.778], xstyle=1,ystyle=1, $
            xmargin=[0,0], ymargin=[0,0], charsize=1.5, $
            xr=[-imsize/2., imsize/2.], yr=[-imsize/2., imsize/2.]
          xyouts, 0.5, 0.95, nam[q], alignment=0.5, charsize=3., $
            color=clr.black, /normal
          xyouts, 0.5, 0.9, ra[q]+'   '+dec[q],$
            alignment=0.5, charsize=2.5, color=clr.black, /normal
          xyouts, 0.5, 0.87, '(J2000)', $
            alignment=0.5, charsize=2., color=clr.black, /normal
          xyouts, 0.5, 0.8, 'NORTH', alignment=0.5, charsize=2.3, /normal, $
            color=clr.black
          if keyword_set(NOTES) then begin  ;; KHRR added -- for extra text
             xyouts, 0.5, 0.83, NOTES, alignment=0.5, charsize=1.5, color=clr.blue, /normal
          endif 

          ;; Circle
          x_oplotcirc, circ/(60.), color=clr.green

          ;; Second circle?
          if keyword_set(TWOCIRC) then begin
              x_oplotcirc, circ/(60.), x0=twocirc[0], y0=twocirc[1], $
                           color=clr.blue
          endif

          ;; additional circles   (MF)
          if keyword_set(ADDCIRC) then begin
             nadd_size=size(ADDCIRC)
             ;;loop over
             for ac=0, nadd_size[1]-1 do begin
                x_oplotcirc, circ/(60.), x0=addcirc[ac,0], y0=addcirc[ac,1], $
                             color=clr.blue
             endfor
          end
          
          ;; Draw a line on the plot
          ;; added by KHRR
          if keyword_set(CONNECTX) then begin
             plots, connectx[0], connecty[0]
             plots, connectx[1], connecty[1], color=clr.red, linestyle=1, /continue
          endif 

          
          x_psclose
      endelse

  endfor

  print, 'x_fndchrt: All done!!'
  
  return
end


