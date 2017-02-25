pro gemini_proc, fil 


  ;; Chip1
  d1 = xmrdfits(fil, 1, head1, /fscale, /silent)
  d2 = xmrdfits(fil, 2, head2, /fscale, /silent)
  d3 = xmrdfits(fil, 3, head3, /fscale, /silent)

  ;; OV subtract
  sz3 = size(d3, /dimensions)
  ov3 = djs_median(d3[1028:*,*],1)
  ovfit = x1dfit(ov3)
  ovi_3 = d3 - replicate(1.,sz3[0]) # ovfit

  sz1 = size(d1, /dimensions)
  ov1 = djs_median(d1[0:20,*],1)
  ovfit = x1dfit(ov1)
  ovi_1 = d1 - replicate(1.,sz1[0]) # ovfit

  sz2 = size(d2, /dimensions)
  ov2 = djs_median(d2[0:20,*],1)
  ovfit = x1dfit(ov2)
  ovi_2 = d2 - replicate(1.,sz2[0]) # ovfit

  ;; output
  sl = strlen(fil)
  outfil = 'OV/ov_'+'1_'+strmid(fil,sl-9)
  mwrfits, ovi_1[30:*,*], outfil, /create
  outfil = 'OV/ov_'+'2_'+strmid(fil,sl-9)
  mwrfits, ovi_2[30:*,*], outfil, /create
  outfil = 'OV/ov_'+'3_'+strmid(fil,sl-9)
  mwrfits, ovi_3[0:1023,*], outfil, /create

  return
end
