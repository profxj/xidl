function gnirs_fit_patternnoise, array, noutput

  sz = size(array)


  indx = findgen(sz[1])
  inflect = fltarr(sz[1])

  sset = bspline_iterfit(indx, array, everyn=25, /SILENT)
  fit = bspline_valu(indx, sset)

  array2 = array - fit

  for i=1, sz[1]-2 do begin
      if (array2[i] LT array2[i-1] AND array2[i] LT array2[i+1] AND $
          array2[i] LT 0) then begin
          inflect[i] = 1
      endif
  endfor
      
  negative = where(inflect EQ 1 AND array2 LT mean(array2) - 0.5 * stddev(array2))

  phase0 = indx[negative[0]]

  pi  = 3.14159
  period = findgen(4) + 7
  amp = fltarr(4)

  for i=0,3 do begin
      sinmodel = sin((2*pi*(indx-phase0))/period[i]+3*pi/2)
      amp[i] = total(sinmodel * array2)
  endfor

  phase = period[where(amp EQ max(amp))]

  phase = [8]

;  plot, indx[0:50], array2[0:50]
;  oplot, indx, 0.5*sin((2*pi*(indx-phase0))/phase[0]+3*pi/2), color=getcolor('green')

;  print, "phase period = ", phase[0], " pixels"

  sz2 = size(negative)

  phase = phase[0]

  phasearr = (indx - phase0) - (floor(indx / phase)*phase)
  minphase = min(phasearr)
  maxphase = max(phasearr)

  phasefit = fltarr(maxphase-minphase+1)

  for i=0, maxphase-minphase do begin
      phasefit[i] = median(array[where(phasearr EQ minphase+i)])
  endfor

  model = fltarr(noutput)
  x = findgen(noutput)
  modelphase = (x - phase0) - (floor(x / phase)*phase)


  for i=0, maxphase-minphase do begin
      model[where(modelphase EQ minphase+i)] = phasefit[i]
  endfor
  
  model = model - median(model)

  return, model

end


