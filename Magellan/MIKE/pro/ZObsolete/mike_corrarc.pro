pro mike_corrarc, xoff, yoff, xopt, yopt, Magf=Magf, i1nm=i1nm, i2nm=i2nm, $
                  plateau=plateau, SLOW=slow, REGION=region, CHK=chk

  if not keyword_set( i1nm ) then i1nm = 'Raw/mb0428.fits'
  if not keyword_set( i2nm ) then i2nm = 'Raw/mb0435.fits'
  if not keyword_set( REGION ) then region = [ 575, 830, 1015L, 1270]

  print, 'Correlating: ', i1nm, ' ', i2nm
  print, 'Using region: ', region

  i1 = xmrdfits(i1nm, /fscale, /silent)
  i2 = xmrdfits(i2nm, /fscale, /silent)
  i1 = i1[region[0]:region[1],region[2]:region[3]]
  i2 = i2[region[0]:region[1],region[2]:region[3]]

  sz = size(i1, /dimensions)

  ;; Cross correlate

  if not keyword_set( MAGF ) then Magf = 8L

  ;; REBIN
  r1 = rebin(i1, sz[0]*Magf, sz[1]*Magf)
  r2 = rebin(i2, sz[0]*Magf, sz[1]*Magf)

  if keyword_set( SLOW ) then begin

      ;; LOOP
      corr_val = fltarr(4*Magf, 2*Magf)
      
      x1_1 = 40 
      x2_1 = sz[0]*Magf - 1 - 40
      y1_1 = 20
      y2_1 = sz[1]*Magf - 1 - 20
      
      sub_r1 = r1[x1_1:x2_1,y1_1:y2_1]
      tot11 = total( sub_r1 * sub_r1 )
      
      for i=0L, 4*Magf-1 do begin  ;; Loop on 2 pix offset in x
          for j=0L, 2*Magf-1 do begin  ;; Loop on 1 pix offset in y
              
              x1_2 = 40 + xoff*Magf + (i - 2*Magf)
              x2_2 = sz[0]*Magf - 1 - 40 + xoff*Magf + (i - 2*Magf)
              y1_2 = 20 + yoff*Magf + (j - Magf)
              y2_2 = sz[1]*Magf - 1 - 20 + yoff*Magf + (j - Magf)
              
              ;; Auto Corr
              sub_r2 = r2[x1_2:x2_2,y1_2:y2_2]
              tot22 = total( sub_r2 * sub_r2 )
              
              ;; Cross-Corr
              corr_val[i,j] = total( sub_r1 * sub_r2 ) / sqrt( tot11 * tot22)
          endfor
      endfor
          
      if keyword_set( CHK ) then xatv, corr_val, /block
      mx = max( corr_val, imx)
      xopt = xoff + ((imx mod (4*Magf)) - 2*Magf) * (1./float(Magf))
      yopt = yoff + (imx / (4*Magf) - Magf) * (1./float(Magf))
      print, 'Offset: ', xopt, yopt
  endif else begin
      ;; FFT
      fft_1 = fft(r1)
      fft_2 = fft(r2)

      corr = fft( fft_1 * conj(fft_2), /inverse)

      ans = double(corr)
      ans2 = shift(ans, 64, 64)

      if keyword_set( CHK ) then $
        xatv, ans2[64-2*Magf:64+2*Magf,64-Magf:64+Magf], /block

      ;; Find max
      mx = max( ans2[64-2*Magf:64+2*Magf,64-Magf:64+Magf], imx)

      ;; Find offset
      xopt = xoff + ((imx mod (4*Magf + 1)) - 2*Magf) * (1./float(Magf))
      yopt = yoff + (imx / (4*Magf+1) - Magf) * (1./float(Magf))

      ;; Sign flip (not exactly sure why)
      xopt = -xopt
      yopt = -yopt

      print, 'Offset: ', xopt, yopt

  endelse

  return
end
