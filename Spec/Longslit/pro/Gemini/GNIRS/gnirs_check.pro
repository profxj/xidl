;PRO NIRI_CHECK

dirs = findfile('Science/', count = nfil)
scidir = x_guilist(dirs)
path =  'Science/' + strcompress(scidir, /rem) + '/'
fils1 = findfile(path + '*.fits*', count = nfil)
fils = fileandpath(fils1, path = path1)
if nfil EQ 0 then message, 'Cannot find files'
file_fits = path + x_guilist(fils)
psfile = repstr(file_fits, '.fits', '.ps')
spawn, 'gv --orientation=Seascape  ' + psfile + '&'

flatfile = findfile('superflat*.fits*', count = nflat)
IF nflat GT 1 THEN thisflat = x_guilist(flatfile) $
ELSE thisflat = flatfile[0]
tset_slits = xmrdfits(thisflat, 1)
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
