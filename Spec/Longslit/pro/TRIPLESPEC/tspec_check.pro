;PRO TSPEC_CHECK

fils = findfile('sci*.fits*', count = nfil)
tels = findfile('tel*.fits*', count = ntel)
IF nfil GT 0 AND ntel GT 0 THEN  fils = [fils, tels] $
ELSE IF nfil GT 0 AND ntel EQ 0 THEN fils = fils $
ELSE IF nfil EQ 0 And ntel GT 0 THEN fils = tels

if nfil EQ 0 then message, 'Cannot find files'
file_fits = x_guilist(fils)

orderfile = findfile('../tspec-orders*.fits', count = nord)
IF nord GT 1 THEN thisorder = x_guilist(orderfile) $
ELSE thisorder = orderfile[0]
tset_slits = xmrdfits(thisorder, 1)
traceset2xy, tset_slits[0], yy1, xx1
traceset2xy, tset_slits[1], yy2, xx2

abba = mrdfits(file_fits, 0, scihdr)
sky_resids = mrdfits(file_fits, 1)
ivar_abba = mrdfits(file_fits, 2)
waveimg = mrdfits(file_fits, 3)
objstruct = mrdfits(file_fits, 4)

mask = double(ivar_abba GT 0.0)
xatv, mask*(abba-sky_resids), min = -200.0, max = 200.0, wvimg = waveimg
xatvplot, objstruct.xpos, objstruct.ypos, psym = 3
xatvplot, xx1, yy1, psym = 3, color = 2
xatvplot, xx2, yy2, psym = 3, color = 2


END
