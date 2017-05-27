function process_arclists,arcnames,lampfilename,grating,anamorph,depth=depth


  if n_elements(depth) eq 0 then depth=3

  scale_pix_per_asec = 8.52

  wid_asec = 1. ; try 0.75?
  wid_pix = wid_asec*scale_pix_per_asec/anamorph

  wid_ang=wid_pix*0.33*(1200./grating)

; playing around...
  if grating lt 1200 then mindist=wid_ang*1.6 else $
      mindist=wid_ang*2.35


; 2.35 yields results that match Doug's original - <4 AA spacing for blends 
; in the 1200-line grating	

  lamps = read_lampfile(lampfilename)
  longxe=0

	oldfile=strpos(lampfilename,'lamp_NIST.dat') ge 0

  narc = n_elements(arcnames)
  for i = 0, narc-1 do begin

    head_arc = headfits(arcnames[i]) ;get FITS header for arc-lamp file
    lampstmp = sxpar(head_arc, 'LAMPS') ; e.g. LAMPS='Ne Ar Kr Xe'
    exptime = sxpar(head_arc, 'EXPTIME')   
; check for long xenon arc
    if exptime gt 20 and strpos(lampstmp,'Xe') ge 0 then longxe=1

	 if i eq 0 then lamps_on = lampstmp else lamps_on = lamps_on+' '+lampstmp
  endfor

; if long xenon arc is in the set, change 'Xe' to 'Xel'
	if longxe then begin
		split=strsplit(lamps_on,/extract)
		wh=where(split eq 'Xe',badct)
		if badct gt 0 then split[wh]='Xel'
		lamps_on=strjoin(split,' ')
	endif

  element = strsplit(lamps_on, /extract) ;list of lamps turned on
  print, 'Arclamps: ', element

  lamp_on = lonarr(n_elements(lamps))
  for i=0, n_elements(lamps)-1 do begin
    j = where(lamps[i].element eq element,yes)
    lamp_on[i] =  yes
  endfor

  lamps = lamps[where(lamp_on gt 0)]

  whlong=where(lamps.element eq 'Xel',longct)
  if longct gt 0 then lamps[whlong].element = 'Xe'



; find blends/contaminated lines.  First, classify lines by strength:
;   'GOOD/BLEND' = 1, 'FAIR/TILT/FBLEND'- 2, etc.
	strength=lamps.intensity*0

	if oldfile eq 0 then begin

	   whc1=where(lamps.quality eq 'BLEND' or $
		  lamps.quality eq 'GOOD' or lamps.quality eq 'BAD' or $
		  lamps.quality eq 'TILT' or lamps.quality eq 'TBLEND',c1ct)
	   if c1ct gt 0 then strength[whc1]=1
	   whc2=where(lamps.quality eq 'FAIR' OR $
		lamps.quality eq 'FBLEND',c2ct)
	   if c2ct gt 0 then strength[whc2]=2
	   whc3=where(lamps.quality eq 'POOR' OR $
		lamps.quality eq 'PBLEND',c3ct)
	   if c3ct gt 0 then strength[whc3]=3
	   whc4=where(lamps.quality eq 'WEAK' OR $
		lamps.quality eq 'WBLEND',c4ct)
	   if c4ct gt 0 then strength[whc4]=4
	   prefixes=['','F','P','W']	

	for i=1,(depth < 4) do begin
		whclass=where(strength eq i,nlines)
		for j=0,nlines-1 do begin
			cenlambda=lamps[whclass[j]].lambda
			whinterfere=where(lamps.lambda gt $
    (cenlambda-mindist) AND lamps.lambda lt (cenlambda+mindist) $
    AND strength le i+1,blendct)	
			if blendct gt 1 then begin
;				print,'old quality: ',lamps[whclass[j]].quality
				lamps[whclass[j]].quality = $
					prefixes[i-1]+'BLEND'
;				print,'new quality: ',prefixes[i-1]+'BLEND'
;				print,lamps[whclass[j]].lambda
			endif
		endfor
	endfor

;	lamps.good = (lamps.quality eq 'GOOD' OR lamps.quality eq 'FAIR')
	lamps.good = (lamps.quality eq 'GOOD'); OR lamps.quality eq 'FAIR')

	return,lamps[where(strength le depth)]

      endif else return,lamps

end













