
reduxpath = '/Users/joe/LUCI_redux_0415/'
slitfile = reduxpath + 'slits-luci.20150413.0060.fits'
path = '/Users/joe/LUCI_data_april15/20150414/'
TELLURIC = 0                  ;path + 'SOFI_0091.fits'
;  TELLURIC = path + 'SOFI_0099.fits'
;init = 107
;nseq = 1
init = 136
nseq = 4
IF KEYWORD_SET(TELLURIC) AND nseq NE 1 THEN message, 'Do not stack Tellurics'
even = 2*lindgen(nseq)
odd = 2*lindgen(nseq) + 1
anum = init + even
bnum = init + odd

prefix = 'luci.20150414.'
afiles = path + prefix + string(anum, FORMAT = '(I4.4)') + '.fits'
bfiles = path + prefix + string(bnum, FORMAT = '(I4.4)') + '.fits'

hdr = headfits(afiles[0])
object = strcompress(sxpar(hdr, 'OBJECT'), /rem)

targdir = reduxpath + '/' + object + '/'
IF FILE_TEST(targdir, /DIR) EQ 0 THEN spawn, 'mkdir ' + targdir

IF KEYWORD_SET(TELLURIC) THEN scifile = tardir + $
                                        'tel-' + fileandpath(afiles[0]) $
ELSE scifile = targdir + 'sci-' + fileandpath(afiles[0]) 
;;wavefile = 'wave-' + fileandpath(afiles[0]) 

chk = 1
tset_slits = mrdfits(slitfile, 1)
slitmask = long_slits2mask(tset_slits)

diff = luci_skysub(afiles, bfiles, tset_slits $
                   , ivar = ivar, waveimg = waveimg, sky = sky $
                   , obj_pos = obj_pos, obj_neg = obj_neg, chk = chk)
;; Loop over objects and extract

plate_scale = 0.250D

IF KEYWORD_SET(obj_pos) THEN npos = n_elements(obj_pos) ELSE npos = 0L
IF KEYWORD_SET(obj_neg) THEN nneg = n_elements(obj_neg) ELSE nneg = 0L
final_struct = 0

FOR iobj = 0L, npos -1L DO BEGIN
   extract = gnirs_extract(diff, ivar, waveimg, slitmask $
                           , sky, obj_pos[iobj], plate_scale)
   final_struct = struct_append(final_struct, extract)
ENDFOR
FOR iobj = 0L, nneg -1L DO BEGIN
   extract = gnirs_extract(-diff, ivar, waveimg, slitmask $
                           , sky, obj_neg[iobj], plate_scale)
   final_struct = struct_append(final_struct, extract)
ENDFOR
isort = sort(final_struct.xfracpos)
final_struct = final_struct[isort]
final_struct.OBJID = lindgen(n_elements(final_struct)) + 1L

;----------
; Write output file
splog, 'Writing FITS file ', scifile
mwrfits, float(diff), scifile, scihdr, /create
mwrfits, float(ivar), scifile
mwrfits, float(waveimg), scifile
mwrfits, final_struct, scifile

;IF npos NE 0 OR nneg NE 0 THEN $
;   niri_plotsci, scifile, hard_ps = repstr(scifile, '.fits', '.ps') $
;  , box = keyword_set(TELLURIC)

;splog, 'Elapsed time = ', systime(1)-t0, ' sec'


END
