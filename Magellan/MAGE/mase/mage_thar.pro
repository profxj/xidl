pro mage_thar, img_num, std_num, sensfunc=sensfunc

number = strtrim(img_num,2)
len = strlen(number)
while len lt 4 do begin
   number = '0'+number
   len = strlen(number)
endwhile
number_std = strtrim(std_num,2)
len = strlen(number_std)
while len lt 4 do begin
   number_std = '0'+number_std
   len = strlen(number_std)
endwhile

if not keyword_set(sensfunc) then sensfunc = 'std.sav'

 tset_slits = mrdfits("Orders.fits",1)
;; ------
;; Expand slit set to get left and right edge
traceset2xy, tset_slits[0], rows, left_edge
traceset2xy, tset_slits[1], rows, right_edge
trace = (left_edge + right_edge)/2.0D

arcimg=xmrdfits('Arcs/Arc'+number+'.fits',hdr)
sz = size(arcimg, /dimen)
nx = sz[0]
ny = sz[1]

waveimg = xmrdfits('Arcs/ArcImg'+number+'.fits')
ordr_str=mrdfits("OStr_mage.fits", 1)
;; read in std trace
obj_strct=xmrdfits('Object/ObjStr'+number_std+'.fits',1) 
scihdr=xheadfits('Raw/mage'+number+'.fits')
dum0=0*arcimg
dum1=dum0 + 1.0d
mage_echextobj,arcimg,dum1,scihdr,dum0,dum1,waveimg,tset_slits,obj_strct $
               ,/BOXCAR,BOX_RAD=2.0,outfil='Final/debug_arc'+number+'.fits' $
               ,/NOVACHELIO
obj_strct.exp=float(sxpar(scihdr,'EXPTIME'))

mage_flux, sensfunc, obj_strct, rej=0.05

res = 299792.458d/4100.d*float(sxpar(scihdr,'SLITNAME'))
mage_combspec, obj_strct, fspec
outflux = 'FSpec/ThAr_F.fits'
outerr  = 'FSpec/ThAr_E.fits'
combname = 'FSpec/ThAr_comb.fits'
mage_1dspec, fspec, outflux, outerr, combname $
             , hdr=scihdr, resvel=res
END
