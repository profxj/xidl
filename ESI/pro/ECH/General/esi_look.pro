
;; identify the esi structure
esifils = findfile('esi*.fits*', count = nesifil)
esi_fits = x_guilist(esifils)
;; read in the esi strucutre
esi = esi_ar(esi_fits)

;; find all the science object files
indx = where(esi.flg_anly NE 0 AND esi.mode EQ 2 AND $
             strtrim(esi.type, 2) EQ 'OBJ', nindx)
esi = esi[indx]
IF KEYWORD_SET(FRINGE) THEN BEGIN 
   imgfil = esi_getfil('fringe_fil', subfil = esi.img_root, /name)
   objfiles = 'Extract/FringeObj_'+esi.img_root
ENDIF ELSE BEGIN 
   imgfil = esi_getfil('fin_fil', subfil = esi.img_root, /name)
   objfiles = esi.OBJ_FIL
ENDELSE
file_exist = lonarr(nindx)
FOR kk = 0L, nindx-1L DO BEGIN
    temp = findfile(imgfil[kk] + '*', count = this_count)
    file_exist[kk] = this_count
ENDFOR
keep = where(file_exist, nkeep)
IF nkeep EQ 0 THEN message, 'No reduced files found' $
ELSE BEGIN
   imgfil = imgfil[keep]
   objfiles = objfiles[keep]
   esi = esi[keep]
ENDELSE

gui_string = imgfil + '  --  ' + esi.OBJ
obj_fits = x_guilist(gui_string, indx = ithis)

waveimg   = xmrdfits(esi[ithis].ARC_FIL)
objstruct = xmrdfits(objfiles[ithis], 1)
;; read in the extensions we need
sedg_fil = esi_getfil('sedg_fil', SLIT = esi[ithis].slit $
                      , cbin = esi[ithis].cbin, rbin = esi[ithis].rbin, /name)
slit_edg = xmrdfits(sedg_fil, /silent)
tset_slits = xmrdfits(sedg_fil, 1)
nx = tset_slits[0].DIMS[0]
ny = tset_slits[0].DIMS[1]
ordermask = esi_echordermask(tset_slits)
;;fin_fil  = esi[ithis].IMG_FINAL
fin_fil = imgfil[ithis]
sciimg   = xmrdfits(fin_fil, 0, hdr)
skyimg   = xmrdfits(fin_fil, 2)
modelvar = xmrdfits(fin_fil, 4)
objimg   = xmrdfits(fin_fil, 5)
outmask  = xmrdfits(fin_fil, 6)

sciivar = (modelvar GT 0.0)/(modelvar + 3*(modelvar LE 0.0))

IF KEYWORD_SET(RESIDS) THEN $
  xatv, sqrt(sciivar)*(sciimg-skyimg - objimg)*(ordermask GT 0.0)*outmask $
  , wvimg = waveimg, min = -6.0, max = 6.0 $
ELSE IF KEYWORD_SET(ZEROMASK) THEN BEGIN
    inds = WHERE(ordermask GT 0.0 AND outmask GT 0.0, ngood)
    sky, (sciimg-skyimg)[inds], skymode, skysig $
         , /silent
    xatv, (sciimg-skyimg)*(ordermask GT 0.0)*(outmask GT 0.0), wvimg = waveimg $
          , min = -2.5*skysig, max = 25.0*skysig
ENDIF ELSE BEGIN
    inds = WHERE(ordermask GT 0.0 AND outmask GT 0.0, ngood)
    sky, (sciimg-skyimg)[inds], skymode, skysig $
         , /silent
    xatv, (sciimg-skyimg)*(ordermask GT 0.0), wvimg = waveimg $
          , min = -2.5*skysig, max = 25.0*skysig
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
xatvplot, slit_edg[*, *, 0], rows, psym = 3, color = 2
xatvplot, slit_edg[*, *, 1], rows, psym = 3, color = 2
xatvplot, objstruct.trace[0:ny-1L], rows $
          , psym = 3



END
