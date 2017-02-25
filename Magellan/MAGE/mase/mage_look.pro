
;; identify the mage structure
magefils = findfile('mage*.fits*', count = nmagefil)
mage_fits = x_guilist(magefils)
;; read in the mage strucutre
mage = xmrdfits(mage_fits,1)

;; find all the science object files
indx = where((strcompress(mage.exptype,/rem) EQ 'SCIENCE' OR  strcompress(mage.exptype,/rem) EQ 'BRIGHT') AND $
             mage.obj_id GT 0,nindx)
mage = mage[indx]
finfil = 'Final/f_' + mage.fitsfile
file_exist = lonarr(nindx)
FOR kk = 0L, nindx-1L DO BEGIN
   temp = file_search(finfil[kk] + '*', count = this_count)
   file_exist[kk] = this_count
ENDFOR
keep = where(file_exist, nkeep)
IF nkeep EQ 0 THEN message, 'No reduced files found' $
ELSE finfil = finfil[keep]
mage = mage[keep]

gui_string = finfil + '  --  ' + mage.OBJECT
obj_fits = x_guilist(gui_string, indx = ithis)
finfil=finfil[ithis]
tmp=strcompress(strsplit(mage[ithis].ARCFILE,'mage',/extract),/rem)
wavefile='Arcs/ArcImg' + tmp[0]
waveimg   = 10.d^xmrdfits(wavefile)
tmp2=strcompress(strsplit(mage[ithis].FITSFILE,'mage',/extract),/rem)
objfile = 'Object/ObjStr' + tmp2[0]
objstruct = xmrdfits(objfile, 1)
;; read in ordermask and slit structure
ordermask=mrdfits('Orders.fits',0)
tset_slits = mrdfits('Orders.fits',1)
traceset2xy, tset_slits[0], rows, left_edge
traceset2xy, tset_slits[1], rows, right_edge

nx = tset_slits[0].DIMS[0]
ny = tset_slits[0].DIMS[1]

sciimg   = xmrdfits(finfil, 0, hdr)
skyimg   = xmrdfits(finfil, 2)
modelivar = xmrdfits(finfil, 4)
objimg   = xmrdfits(finfil, 5)
outmask  = xmrdfits(finfil, 6)


IF KEYWORD_SET(RESIDS) THEN $
  xatv, sqrt(modelivar)*(sciimg-skyimg - objimg)*(ordermask GT 0.0)*outmask $
  , wvimg = waveimg, min = -6.0, max = 6.0 $
ELSE BEGIN     
    inds = WHERE(ordermask GT 0.0 AND outmask GT 0.0, ngood)
    sky, (sciimg-skyimg)[inds], skymode, skysig $
         , /silent
    xatv, (sciimg-skyimg)*(ordermask GT 0.0)*outmask, wvimg = waveimg $
          , min = -2.5*skysig, max = 25.0*skysig
   ; xatv, (sciimg-skyimg), wvimg = waveimg $
   ;       , min = -2.5*skysig, max = 25.0*skysig
ENDELSE

IF KEYWORD_SET(PLOTMASK) THEN BEGIN
    ximg = findgen(nx) # replicate(1.0, ny)
    yimg = replicate(1.0, nx) # findgen(ny)
    maskpix = WHERE(ordermask GT 0 AND outmask EQ 0, nmask)
    PLOTSYM, 0, 1.0, THICK = 1.0
    IF nmask GT 0 THEN xatvplot, ximg[maskpix], yimg[maskpix], psym = 8 $
      , color = 4
ENDIF

nord = n_elements(objstruct)
rows = findgen(ny) # replicate(1.0, nord)
xatvplot, left_edge, rows, psym = 3, color = 2
xatvplot, right_edge, rows, psym = 3, color = 2
xatvplot, objstruct.trace[0:ny-1L], rows, psym = 3

END
