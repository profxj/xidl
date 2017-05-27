; code to analyze DEIMOS data taken with Greg's arcstep shell script.  We
; step the grating through a range of tilts and watch arclines go
; across the chip to look for defects. 
;
; The conclusion was that there are no defects at the 1-3 mA level. 
; The chips look very good!

; 2002-Oct-04 Written by Doug Finkbeiner at Keck.
;
pro arcstep_planfile

  for frame=51, 63 do begin 
     
     str = string(frame, format='(I4.4)')
     spawn, 'mkdir -p '+str

     openw, wlun, str+'/arcstep.plan', /get_lun
     printf, wlun, 'MASK: 3205 - 1200 line'
     printf, wlun, 'RAWDATA: 2002oct04'
     printf, wlun, 'FLATNAME:  d1004_0018.fits'
;  printf, wlun, 'FLATNAME:  d1004_0019.fits'
;  printf, wlun, 'FLATNAME:  d1004_0020.fits'
     printf, wlun, 'ARCNAME:   d1004_'+str+'.fits'
  endfor 
  return
end

pro arcstep

  for frame=50, 63 do begin 
     str = string(frame, format='(I4.4)')
     cd, str
     deimos_mask_calibrate, 'arcstep.plan', chip=6, /noplo
     cd, '../'
  endfor 

  return
end

pro arcstep_dline, d

  deep_dir = getenv('DEEP_DIR')+'/'
  if deep_dir eq '/' then message, 'You need to set $DEEP_DIR!'
  lampfilename = deep_dir+'spec2d/etc/lamp_NIST.dat'
  lamps = read_lampfile(lampfilename)

;  head_arc = headfits(arcnames[0])  ;get FITS header for arc-lamp file
;  lamps_on = sxpar(head_arc, 'LAMPS') ; e.g. LAMPS='Ne Ar Kr Xe'

  lamps_on='Ne Ar Kr Xe'
  element = strsplit(lamps_on, /extract) ;list of lamps turned on
  print, 'Arclamps: ', element

  lamp_on = lonarr(n_elements(lamps))
  for i=0, n_elements(lamps)-1 do begin 
    j = where(lamps[i].element eq element,yes)
    lamp_on[i] =  yes
  endfor
  lamps = lamps[where(lamp_on gt 0)]  
  help, lamps

  frames = lindgen(14)+50
  d = dblarr(537, 34, n_elements(frames))
  for i=0, n_elements(frames)-1 do begin 
     frame = frames[i]
     str = string(frame, format='(I4.4)')
     dline = readfits(str+'/dline6.fits')
     pixm = dline*0
     mask = dline NE 0
     flist = findfile(str+'/calib*.fits', count=ncal)
     for j=0, ncal-1 do begin 
        a = mrdfits(flist[j], 1)
        wavmeas = lamps.lambda-dline[*, j]
        pix = traceset2pix(a, wavmeas, /silent)
        ymid = (size(pix))[2]/2
        pix = pix[*, ymid]
        pixm[*, j] = pix
     endfor 
     d[*, *, i] = pixm
  endfor 
; dline is linelist minus measured
  return
end


pro arcstep_compare, d, h, line

  nslit = 34
  nexp = 14
  dd = d[[426, 428, 429, 437, 439, 456, 459, 460, 477], *, *]

  e = reform(d[line,*,*], nslit, nexp)
  f = e
  for i=0, nslit-1 do begin 
     foo = linfit(dindgen(nexp), e[i, *], yfit=yfit)
     f[i, *] = e[i, *]-yfit
  endfor 

;  m = total(f, 1)/nslit
;  g = f-(fltarr(nslit)+1)#m
;window, 0
;plot, -reverse(m)

;OVERRIDE
  g = e-total(dd,1)/9
  h = g
  for i=0, nslit-1 do begin 
     foo = linfit(dindgen(nexp), g[i, *], yfit=yfit)
     h[i, *] = g[i, *]-yfit
  endfor 
  per = 68.5
  splot,e,h,ps=1, yr=[-.1, .1], /nodata
  for i=0, nslit-1 do soplot, e[i, *], h[i, *], color=(i mod 6)+1
  for i=0, 60 do soplot, [1, 1]*i*per, [-1, 1]

return
stop
  mm = total(e, 2)/nexp
  help, mm#(fltarr(nexp)+1)
  g = f-mm#(fltarr(nexp)+1)


  return
end
