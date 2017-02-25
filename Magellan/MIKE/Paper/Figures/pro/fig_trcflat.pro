;+ 
; NAME:
; fig_trcflat
;    Version 1.1
;
; PURPOSE:
;    Given an array of Obj name, RA, and DEC create a set of postscript
;  finding charts using Barlows showfits routine. 
;
; CALLING SEQUENCE:
;  x_fndchrt, targlist, OUTDIR=, IMSIZE=, SURVEY=
;
; INPUTS:
;  targlist  -- ASCII file containing  (QSO,  RA,  DEC)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   imsize - Arcmin of image (default is 5')
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
;
; REVISION HISTORY:
;   08-Dec-2004 Written by JXP
;-
;------------------------------------------------------------------------------
;fig_galimg, [650, 1350, 750, 1450]
pro fig_trcflat, infil, sub

  common xcommon_color, r_vector, g_vector, b_vector
  xicmm_colors

  ;; info
  if not keyword_set(PSFIL) then psfil='fig_trcflat.ps'
  if not keyword_set(TRC_FIL) then trc_fil = $
    getenv('MIKE_PAP')+'Flats/Flat_B_01_T.fits'
  if not keyword_set(OSTRFIL) then ostrfil = $
    getenv('MIKE_PAP')+'Flats/OStr_B_01.fits'

  ordr_str = xmrdfits(ostrfil, 1, /silent)
  nordr = n_elements(ordr_str)

  ;; FITS Img
  img = xmrdfits(trc_fil,/silent)
;  img = img*20.
  sz = size(img,/dimensions)

  ;; Normalize the orders
  omsk = x_ordermask(sz[0], sz[1], ordr_str, trim=3.)
  ex_flat =  x_extractarc(img, ordr_str)

;  img[where(omsk LT 0)] = 0.
;  img[where(omsk GT 110)] = 0.

;  for qq=0L, nordr-1 do begin
;      ;; Fit
;      bset = bspline_iterfit(findgen(sz[1]),ex_flat[*,qq], yfit=yfit,$
;                             everyn=50) 
;      for jj=0L,sz[1]-1 do begin
;          gd = where(omsk[*,jj] EQ ordr_str[qq].order, ngd)
;          if ngd NE 0 then img[gd,jj] = img[gd,jj] / yfit[jj]
;      endfor
;  endfor
  sv_img = img
      

  sub = [1, 1000, 1248L, 2047]
  if not keyword_set( SUB ) then sub = [0L, sz[0]-1, 0, sz[1]-1]
  img = img[sub[0]:sub[1],sub[2]:sub[3]]
  sz = size(img,/dimensions)


  set_plot, 'x'

  ;; PSFILE
;  psfil = 'Figures/fig_pks0405.ps'
;  x_psopen, psfil, /portrait
  state = { $
            ncolors: 0L, $
            brightness: 0.01, $
            contrast: 0.05 $
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

  ;; PS
  aspect = 1.
  tvlct, rr, gg, bb, 8, /get
  
  forminfo = cmps_form(cancel = canceled, create = create, $
                       aspect = aspect, $;parent = state.base_id, $
                       /preserve_aspect, $
                       xsize = 10.0, ysize = 8.0 * aspect, $
                       xoffset=1.0, $
                       /color, /encapsulated, /landscape, $
                       /nocommon, papersize='Letter', $
                       bits_per_pixel=8, $
                       filename = psfil, $
                       button_names = ['Create PS File'], /initialize)

  tvlct, rr, gg, bb, 8
  screen_device = !d.name


  set_plot, 'ps'
;  
  device, _extra = forminfo
  device, /times,isolatin=1

  tvlct, rr, gg, bb, 8, /get
  
  rn = congrid(rr, 248)
  gn = congrid(gg, 248)
  bn = congrid(bb, 248)
;
  tvlct, temporary(rn), temporary(gn), temporary(bn), 8

;  newdisplay = bytarr(xsize, ysize)

          
  ;; Display image
;  mx = max(img, min=mn)
;  med = median(img)
;  sig = stddev(img)
;  pltmax = (med + (3.0*sig)) < mx
;  pltmax = 1.83
;  pltmin = (med - (3.0*sig)) > mn
  pltmin = 0.
;  pltmax = 1.2
  pltmax = 10000.
;  print, pltmin, pltmax
  display_image = bytscl(img, min=pltmin, max=pltmax, /nan, $
                         top=state.ncolors-1) + 8B
;  dimage = display_image
  dimage = bytscl(display_image, top = 247, min=8, max=(!d.table_size-1)) + 8

  ;; TV
;  ntv = 200L
;  tv_image = congrid(display_image, ntv, ntv)
  tv, dimage, 0., 0., xsize=1., /norma
  plot, [0], [0], /nodata, position=[0,0,1,1], $
    xrange=[sub[0],sub[1]], yrange=[sub[2],sub[3]], xstyle=5, ystyle=5, /noerase


  ;; Traces
  dumy = findgen(n_elements(ordr_str[0].lhedg))
  ;; Plot
  xoff = 0.
  for q=0L, n_elements(ordr_str)-1 do begin
      oplot, ordr_str[q].lhedg + xoff, dumy, color=1, thick=3
      oplot, ordr_str[q].rhedg + xoff, dumy, color=2, thick=3
  endfor

  ;; Zoom in 
  sub = [411, 470, 1001, 1160L]
  dx = sub[1]-sub[0] + 1
  dy = sub[3]-sub[2] + 1
  if not keyword_set( SUB ) then sub = [0L, sz[0]-1, 0, sz[1]-1]
  img = sv_img[sub[0]:sub[1],sub[2]:sub[3]]
;  img = congrid(img, dx*3, dy*3)
  sz = size(img,/dimensions)
  pltmin = 0.
  pltmax = 5000.
;  print, pltmin, pltmax
  display_image = bytscl(img, min=pltmin, max=pltmax, /nan, $
                         top=state.ncolors-1) + 8B
;  dimage = display_image
  dimage = bytscl(display_image, top = 247, min=8, max=(!d.table_size-1)) + 8
  xpos = 0.05
  ypos = 0.05
  xsz = 0.3
  ysz = xsz*(dy/dx)
  tv, dimage, xpos, ypos, xsize=xsz, ysize=ysz, /normal

  plot, [0], [0], $
        position=[xpos, ypos, xpos+xsz, ypos+ysz], $
;        position=[xpos,ypos, xpos+xsz,ypos+(dy/dx)*xsz], /device, $
        xmarg=[0,0], ymarg=[0,0], $
        color=1, $
        xrange=[sub[0],sub[1]], yrange=[sub[2],sub[3]], $
        xstyle=5, ystyle=5, /noerase

  ;; Boundary
  oplot, [sub[0],sub[0]], [sub[2],sub[3]], color=0, thick=2, linestyl=2
  oplot, [sub[1],sub[1]], [sub[2],sub[3]], color=0, thick=2, linestyl=2 
  oplot, [sub[0],sub[1]], [sub[2],sub[2]], color=0, thick=2, linestyl=2
  oplot, [sub[0],sub[1]], [sub[3],sub[3]], color=0, thick=2, linestyl=2

  ;; Traces
  dumy = findgen(n_elements(ordr_str[0].lhedg))
  ;; Plot
  xoff = 0.
  for q=0L, n_elements(ordr_str)-1 do begin
      oplot, ordr_str[q].lhedg + xoff, dumy, color=1, thick=4
      oplot, ordr_str[q].rhedg + xoff, dumy, color=2, thick=4
  endfor

  close, /all
  device, /close
  set_plot, 'x'
  device, decompose=1
  !p.font = -1

  print, 'fig_trcflat: All done!!'

  return
end


