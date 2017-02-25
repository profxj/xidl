path = '/Users/joe/DATA/SOFI_DATA/2011-09-22/'
TELLURIC = 0                  ;path + 'SOFI_0091.fits'
;;TELLURIC = path + 'SOFI_0099.fits'
;init = 107
;nseq = 1
init = 95
;init = 145
nseq = 4
IF KEYWORD_SET(TELLURIC) AND nseq NE 1 THEN message, 'Do not stack Tellurics'
even = 2*lindgen(nseq)
odd = 2*lindgen(nseq) + 1
anum = init + even
bnum = init + odd

afiles = path + 'SOFI_' + string(anum, FORMAT = '(I4.4)') + '.fits'
bfiles = path + 'SOFI_' + string(bnum, FORMAT = '(I4.4)') + '.fits'

IF KEYWORD_SET(TELLURIC) THEN scifile = 'tel-' + fileandpath(afiles[0]) $
ELSE scifile = 'sci-' + fileandpath(afiles[0]) 
;;wavefile = 'wave-' + fileandpath(afiles[0]) 
reduxpath = '//Users/joe/SOFI_pipeline/'
hdr = headfits(afiles[0])
object = strcompress(sxpar(hdr, 'OBJECT'), /rem)
targdir = reduxpath + '/' + object
IF FILE_TEST(targdir, /DIR) EQ 0 THEN spawn, 'mkdir ' + targdir

sofi_proc, afiles[0], imag
dims = size(imag, /dim)
nx = dims[0]
ny = dims[1]
slitfile = reduxpath + 'slitmask.fits'
tset_slits = niri_slitset(nx, ny)
slitmask = long_slits2mask(tset_slits)
mwrfits, slitmask, slitfile, /create
mwrfits, tset_slits, slitfile
waveqafile = reduxpath + 'waveQA.ps'

chk = 1
luci_reduce_work, afiles, bfiles, slitfile, scifile, waveqafile $
                  ,  TELLURIC = TELLURIC, CHK = CHK, WVCHK = WVCHK $
                  , FILESTD = FILESTD1 $
                  , pixflatfile = pixflatfile $
                  , illumflatfile = illumflatfile $
                  , darkfile = darkfile, /SOFI


END
