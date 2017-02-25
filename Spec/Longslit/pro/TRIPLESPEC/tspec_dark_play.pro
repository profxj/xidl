
path = '/Users/joe/DATA/P200_DATA/2011418/'
prefix = 'tspec110418_'
darkind = 45L + lindgen(15)
dark_str = string(darkind, FORMAT = '(I4.4)')

filenames = path + prefix + dark_str + '.fits'
superdark = niri_superdark(filenames, /TSPEC)


image = superdark
dims = size(image, /dim)
nx = dims[0]
ny = dims[1]
xarr = findgen(nx) # replicate(1.0, ny)
scatt_shift = 2L
tset_slits = xmrdfits('tspec-orders.fits', 1)
scatt_model = tspec_scattlight(image, tset_slits, scatt_shift = scatt_shift)

slitmask = long_slits2mask(tset_slits)
slitmask_left  = long_slits2mask(tset_slits, xshift = -abs(scatt_shift))
slitmask_right = long_slits2mask(tset_slits, xshift = abs(scatt_shift))
scattmask = (slitmask_left EQ 0 AND slitmask_right EQ 0 AND slitmask EQ 0)

scattpix = WHERE(scattmask AND image GT 1.0 AND image LT 5d4 AND $
                 (xarr GT 0 AND xarr LT (nx-1)))
scatt_amp = total(image[scattpix]*scatt_model[scattpix])/ $
            total(scatt_model[scattpix]^2)
image2 = image - scatt_amp*scatt_model


END
