pro slitview, slitno

  maskno = 1145

; -------- Get environment variables, set paths
  deimos_data = getenv('DEIMOS_DATA')+'/'
  if deimos_data eq '/' then message, 'You need to set $DEIMOS_DATA!'

  maskstr = string(maskno, format='(I4.4)')
  maskdir = deimos_data+maskstr+'/'

  slitstr = string(slitno, format='(I3.3)')

  fname = maskdir+'spSlit.'+maskstr+'.'+slitstr+'R.fits'

  a = mrdfits(fname, 1)
  med = median(a.flux)

  im = [[a.ivar], [a.flux-med], [a.skymodel-med], [a.flux-a.skymodel], [a.lambda], [100*a.mask]]

  atv, im, max=800, min=-200
return
end

