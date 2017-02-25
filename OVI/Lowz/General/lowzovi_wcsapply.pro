;; Written by K. Cooksey
;; Modified by JXP  (12/12/2010)
;;  Only modifies the Galaxy structure
pro lowzovi_wcsapply,imgfil, galfil, clobber=clobber,chip3=chip3


;     prs = strsplit(wcsfits[ii],'/',/extract,count=nprs)
     
     ;; Backup Image file
;     origfil = strmid(prs[nprs-1],0,strpos(prs[nprs-1],'wcs.fits'))+'.fits'
;     imgdir = rootdir+prs[0]+'/'+prs[1]+'/'
;     cd,imgdir
;     if not keyword_set(galstronly) then begin
;         test = file_search(origfil+'.bkp',count=ntest)
;         if ntest eq 0 or keyword_set(clobber) then begin
;             spawn,'cp '+origfil+' '+origfil+'.bkp' 
;             
;             ;; Copy new WCS-ed image up
;             nwfil = prs[nprs-2]+'/'+prs[nprs-1]
;             spawn,'cp '+nwfil+' '+origfil
;             img = x_readimg(origfil,head=hdimg,/fscale)
;         endif else stop,'wcsapply: already have a *.bkp file'
;     endif else img = x_readimg(origfil,head=hdimg,/fscale) ;/galstronly

  ;; Read in image
  img = x_readimg(imgfil,head=hdimg,/fscale) ;/galstronly

  ;; Since image flipped and rotated , change pixel location
  gal = mrdfits(galfil,1,hd,/silent)
  ngal = n_elements(gal)
  if keyword_set(chip3) then begin
     ;; flip to make east left and rotate 180d to make north up
     xnw = gal.xypix[0]
     ynw = n_elements(img[0,*]) - gal.xypix[1]
  endif else begin
     ;; Flip to make east left and rotate 90d to make north up
     xnw = n_elements(img[*,0]) - gal.xypix[1]
     ynw = n_elements(img[0,*]) - gal.xypix[0]
  endelse 

  gal.xypix[0] = xnw
  gal.xypix[1] = ynw

  ;; Apply astrometry to pixel locations
  extast,hdimg,astr ; extract astrometry info
  xy2ad,gal.xypix[0],gal.xypix[1],astr,amod,dmod ; convert pix to cdt
  gal.ra = amod
  gal.dec = dmod
  
  aqso = fxpar(hdimg,'CRVAL1') ; qso ra
  dqso = fxpar(hdimg,'CRVAL2') ; qso dec

  ;; Now for the complicated taking
  qso = create_struct({ galsurveystrct },'dra',-1.d,'ddec',-1.d)
  nwgal = replicate(qso,ngal)
  for jj=0L,ngal-1 do begin
     gcirc,1,amod[jj]/15,dqso,aqso/15,dqso,dra
     gcirc,1,aqso/15,dmod[jj],aqso/15,dqso,ddec
     ;; account for sign
     dra = (amod[jj]-aqso)/abs(amod[jj]-aqso)*dra ; arcsec
     ddec = (dmod[jj]-dqso)/abs(dmod[jj]-dqso)*ddec ; arcsec
     tmp = create_struct(gal[jj],'dra',dra,'ddec',ddec)
     nwgal[jj] = tmp
  endfor 

  ;; Add info about QSO
  qso.gal_type = 'QSO'
  qso.ra = aqso
  qso.dec = dqso
  qso.dra = 0.
  qso.ddec = 0.
  nwgal = [nwgal,qso]
  ngal = ngal+1
  
  ;; Notes
  if keyword_set(chip3) then $
     hist=['KLC modify XYPIX (flip, rot.)',$
           'x_curr=x_old ; y_curr=ymx-y_old',$
           'Convert (x,y) to (a,d) with xy2ad',$
           'Current d(ra,dec) relative to QSO CRPIX# w/GCIRC'] $
  else  hist=['KLC modify XYPIX (flip, rot.)',$
              'x_curr=ymx-y_old; y_curr=xmx-x_old',$
              'Convert (x,y) to (a,d) with xy2ad',$
              'Current d(ra,dec) relative to QSO CRPIX# w/GCIRC']
  sxaddhist,hist,hd,/comment
  
  ;; Back up original file
  test = file_search(galfil+'.bkp',count=ntest)
  if ntest eq 0 or keyword_set(clobber) then begin
     spawn,'cp '+galfil+' '+galfil+'.bkp' 

     ;; Write new galaxy file
     mwrfits,nwgal,galfil,hd,/create
  endif else stop,'wcsapply: galaxy file already modified'
  
  return
end
