;name of the images
;array with the images to stack
;thier mask
;the out image
;the out mask
;the headet to add notes


;;;;;;;;;;;;
;NB This procedure uses shiftf
;A test shows that it preserves flux at 0.001 mag level
;-----------------------
;Aperture list:
;[3,6,8,10,20,25,30]
;
;------------------------------------------------
;test.fits    102.41   100.82  -0.09882   -9.343   -9.538   -9.563 -9.576   -9.602   -9.608   -9.616  ok
;test.fits    102.37   100.76  -0.09901   -9.342   -9.538   -9.563 -9.576   -9.602   -9.608   -9.616  ok
;test.fits    102.31   100.80  -0.09893   -9.343   -9.538   -9.564 -9.576   -9.602   -9.608   -9.615  ok
;-----------------------------------------------
;DLA1_R_cut.fits     91.52    90.40  -0.09725   -9.345   -9.538   -9.563   -9.576   -9.601   -9.608   -9.614  ok
;DLA1_R_cut.fits     91.52    90.40  -0.09725   -9.345   -9.538   -9.563   -9.576   -9.601   -9.608   -9.614  err
;DLA1_R_cut.fits     91.52    90.40  -0.09725   -9.345   -9.538   -9.563   -9.576   -9.601   -9.608   -9.614  err
;--------------------------------------------------------------------------;
;;;;;;;;;;;;;;;;;


pro lbtc_getcombined, stackname, outimg, outmask, header, $
                      path=path, noskysub=noskysub, medianout=medianout, $
                      creject=creject, currside=currside, currchip=currchip, $
                      manshift=manshift



  ;*****************************************
  ;Do some introductury stuff
  ;******************************************

  nstack=n_elements(stackname)
  
   
  ;define useful arrays you will use 
  weight=fltarr(nstack)
  skymean=fltarr(nstack)
  skysigma=fltarr(nstack)
  scale=fltarr(nstack)
  
  



  ;*****************************************
  ;Get the sky properties (no need for shifts here)
  ;******************************************

  ;----------------------------------------
  ;find weigths (background)
  ;----------------------------------------

    
  splog, 'Find weights...'
  for st=0, nstack-1 do begin
     ;open
     fits=mrdfits(stackname[st],/silent)
     thisfits=reform(fits[currchip,*,*])
     undefine, fits
     ;use sky median for sky level (mean is biased sometimes)
     djs_iterstat, thisfits, median=i_skymean, sigma=i_skysigma
     skymean[st]=i_skymean
     skysigma[st]=i_skysigma
     weight[st]=1/(i_skysigma)^2
     splog, stackname[st], ":"
     splog, "Sky value ", i_skymean, '+/-', i_skysigma
     
     undefine, thisfits
     
  endfor
  
  splog, 'Weights ', weight
  
  sxaddpar, header, "COMMENT", $
            strjoin(strcompress(string('SKYMEAN', skymean))),after='COMMENTS'
  sxaddpar, header, "COMMENT", $
            strjoin(strcompress(string('SKYSIGMA', skysigma))),after='COMMENTS'
  


  ;*****************************************
  ;Get the shifts
  ;******************************************

  splog, 'Compute the shifts...'


  if keyword_set(manshift) then begin
        ;go for a manual shift, marking stars by hand
          
     lbtc_manual_align, stackname, xref_star=xref_star,yref_star=yref_star,$
                        xshift=shiftx, yshift=shifty, chipnum=currchip
  
  
  endif else begin
                                ;this is to try an automated
                                ;aligment. Sometimes results are poor.
                                ;Maybe for focus or distorsion at the
                                ;edge of the field
     spawn, 'cp '+getenv('LBT_IMG')+'/SExtract/default.* .'
     fwhm=replicate(1.0,nstack)
        ;force a non zero shift with a median filter
     lbtc_offset, stackname, fwhm, shiftx, shifty, chipnum=currchip, /real, /nozero
  endelse
  
  

  ;*****************************************
  ;Now deal with the scales
  ;******************************************
  
  scale[0]=1.
  splog, 'Work on reference image ', stackname[0]
  
   ;open the reference image and skysub
  fits=mrdfits(stackname[0],/silent)
  refimage=reform(fits[currchip,*,*])-skymean[0]
  undefine, fits

  ;find scales from sources 
  splog, 'Extract reference catalogue...'
  spawn, 'cp '+getenv('LBT_IMG')+'/SExtract/default.* .'

  mwrfits, refimage, 'ref_frame_tmpBJHL.fits', /create, /silent
  spawn, 'sex ref_frame_tmpBJHL.fits'
  spawn, 'rm -fr ref_frame_tmpBJHL.fits'
  ;sort by magnitude
  spawn, 'sort -nk4 tmp.dat -o tmp.dat'

  readcol, 'tmp.dat', x_ref, y_ref, a1, b1, c1, d1, e1, f1, g1,$
           h1, i1, j1, k1, l1, class_star, xg1, yg1,$
           format='F,F,I,F,F,F,F,F,F,F,F,I,I,F,F,F,F', /silent

  

  ;take the 30 brightest stars 
  nsextra=n_elements(x_ref)
  if(nsextra gt 30) then nsextra=30

  x_ref=x_ref[0:nsextra-1]
  y_ref=y_ref[0:nsextra-1]
  
  ;compute reference magnitude
  X_APER, refimage, x_ref, y_ref, flux_ref, errap, sky, $
        skyerr, 1., 10., [15.,20.],[0.,600000.], /flux, SETSKYVAL=0., /silent

  ;free stuff
  undefine, refimage

  ;now go for all the other and find scales
  
  for st=1, nstack-1 do begin
     
        ;open the image, skysub and shift
     fits=mrdfits(stackname[st],/silent)
     refimage=reform(fits[currchip,*,*])-skymean[st]
     refimage=shiftf(refimage,shiftx[st],shifty[st])
     undefine, fits


     X_APER, refimage, x_ref, y_ref, flux_comp, errap, sky, $
        skyerr, 1., 10., [15.,20.],[0.,600000.],/flux, SETSKYVAL=0., /silent
     compar=flux_ref/flux_comp
     ;get good values
     good_nan=where(finite(compar,/nan)-1,ngd)
     if(ngd gt 0) then scale[st]=median(compar[good_nan]) else begin
        splog, 'No good stars found. Set scale to 1 for ', stackname[st]
        scale[st]=1.
     endelse
  endfor

  splog, 'Computed scales ', scale
  sxaddpar, header, "COMMENT",strjoin(strcompress(string('SCALES', scale))),$
            before='DATASUM'

          
 ;****************************************
 ;get the final median (use it late for cosmic)
 ;***************************************  
  
  splog, 'Compute median...'
  fits=mrdfits(stackname[0],/silent)
  nxf=n_elements(fits[0,*,0])
  nyf=n_elements(fits[0,0,*])
  store=fltarr(nstack,nxf,nyf)
  undefine, fits

  for nim=0, nstack-1 do begin
     fits=mrdfits(stackname[nim],/silent)
     ;shift and sky sub
     shiftimg=shiftf(reform(fits[currchip,*,*]),shiftx[nim],shifty[nim])
     undefine, fits
     if ~keyword_set(noskysub) then shiftimg=temporary(shiftimg)-skymean[nim]
     ;apply scale here for median only
     store[nim,*,*]=shiftimg[*,*]*scale[nim]
     undefine, shiftimg
  endfor
  medianout=djs_median(store,1)
  undefine, store
  
  ;xatv, medianout, /block

  ;***************************************************
  ;Build the mean and mask
  ;The mean is built accumulating terms
  ;****************************************************
  
  mean_num=fltarr(nxf,nyf)
  mean_den=fltarr(nxf,nyf)
  outmask=fltarr(nxf,nyf)
  
  ;make a check image for clening CR rejection
  splog, 'Get object mask...'
  mwrfits, medianout, 'tmpMKDEONKOE.fits', /create
  spawn, 'sex tmpMKDEONKOE.fits'
  tmpfits=mrdfits('check.fits')
  good_obj=where(tmpfits GT 0.,ngd)
  file_delete,  'tmpMKDEONKOE.fits'
  undefine, tmpfits


  ;loop over images
  for nim=0, nstack-1 do begin
     splog, 'Build mean/mask for ', stackname[nim]
     ;open, shift and sky sub
     fits=mrdfits(stackname[nim],/silent)
     ;shift and sky sub
     shiftimg=shiftf(reform(fits[currchip,*,*]),shiftx[nim],shifty[nim])
     undefine, fits
     if ~keyword_set(noskysub) then shiftimg=temporary(shiftimg)-skymean[nim]
     ;xatv,  shiftimg, /block
 
     ;now get the mask
     mask=shiftimg-shiftimg+1.
     lbtc_cosmicworld, shiftimg, mask, creject=creject, currside=currside, $
                       currchip=currchip, medianimg=medianout,$
                       scale=scale[nim], shift=[shiftx[nim],shifty[nim]],$
                       objectpix=good_obj
 
     ;xatv,  shiftimg, /block
     ;xatv,  mask, /block
     

     ;add to the mean, mask
     
     mean_num=temporary(mean_num)+shiftimg*mask*scale[nim]*weight[nim]
     mean_den=temporary(mean_den)+mask*weight[nim]
     outmask=temporary(outmask)+abs(mask-1)
     
     ;xatv, mean_num, /block
     ;xatv, mean_den, /block
     ;xatv, outmask, /block


     ;at the end kill what allocated in the for loop
     undefine, shiftimg, mask


  endfor

  
 ;create final mean 
 ;avoid nan
 outimg=mean_num-mean_num
 nozero=where(mean_den gt 0, numnz)
 if(numnz gt 0) then outimg[nozero]=mean_num[nozero]/mean_den[nozero]
  
 ;free mem
 undefine, mean_num, mean_num, tota_mas
  
   
 ;;Feb 2010. Test. This code preserves the photometric properties of
 ;;the images.
 ;;
 ;;


end


  
