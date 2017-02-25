pro kwater

;  openr, 10, "HITRAN.txt"
  openr, 10, "atrans.txt"

  wavelength  = fltarr(28000)
  radiance = fltarr(28000)

  i = 0L
  v1=0.0d
  v2=0.0d
  while (NOT EOF(10)) do begin
     readf, 10, v1, v2
;     wavelength[i]= 1./ v1 * 1e8
     wavelength[i]= v1 * 1e4
     radiance[i] = v2
     i++
  endwhile
  
  close, 10

  set_plot, "ps"
  device, file="FIRE_K_HITRAN.ps", /color, /landscape
  !p.multi=[0,1,2]
;  plot, wavelength, 1-radiance, xrange=[23500,24000], /xsty

;;;;;;;;;;;;;;;;;;;

  openr, 10, "HITRAN.txt"
;  openr, 10, "atrans.txt"

  wavelength  = fltarr(200000)
  radiance = fltarr(200000)

  i = 0L
  v1=0.0d
  v2=0.0d
  while (NOT EOF(10)) do begin
     readf, 10, v1, v2
     wavelength[i]= 1./ v1 * 1e8
;     wavelength[i]= v1 * 1e4
     radiance[i] = v2
     i++
  endwhile
  
  close, 10

  dlam = 500.
  totlam = 2000.
  nlam = totlam/dlam

  colors=getcolor(/load)

  ymax=[2e-11,5e-11,6e-10,1e-9]

  for ilam=0, nlam-1 do begin

     if (ilam EQ 0) then begin
        plot, wavelength, radiance, xrange=[23000.+ilam*dlam,23000.+(ilam+1)*dlam], /xsty, yrange=[0,ymax[ilam]], title="HITRAN Emission Spectrum"
     endif else begin
        plot, wavelength, radiance, xrange=[23000.+ilam*dlam,23000.+(ilam+1)*dlam], /xsty, yrange=[0,ymax[ilam]]
     endelse

     spec = radiance[where(wavelength GT 23000.+ilam*dlam AND $
                     wavelength LT 23000.+(ilam+1)*dlam)]

     wv = wavelength[where(wavelength GT 23000.+ilam*dlam AND $
                     wavelength LT 23000.+(ilam+1)*dlam)]

     x_fndpeaks, spec, pks_i, /all, nsig=1.5
     pkwv = wv[pks_i]
     pkfx = spec[pks_i]

     oplot, pkwv, pkfx, psym=1, color=colors.red

     xyouts, pkwv, pkfx*1.05, strtrim(pkwv,2), orientation=90

  endfor

  device, /close
  set_plot, "x"

  !p.multi=0
  stop

end
