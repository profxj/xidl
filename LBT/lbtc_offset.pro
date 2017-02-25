
;chipnum --> number of the chip to consider
;smooth --> if set, take the smooth version of the histo image
;nozero --> if set, the code try to reject a null shift using median filter

pro lbtc_offset, files, fwhm, outx, outy, real=real, chipnum=chipnum, nozero=nozero

  
  !EXCEPT=0                     ;Turn off annoying illegal operand errors

   ;close, /all                   ; close all previous stuff!

;---------------------------------------------------------------------------
;  get ready to start
;---------------------------------------------------------------------------


  FWHM_str = strtrim(string(FWHM, '(F4.2)'),2)
   

; set pixel scale
  pixscale = 0.224              ;"/pix
  
  pix = strtrim(string(pixscale),2) 

; prepare output array

  outx = fltarr(n_elements(FWHM))
  outy = fltarr(n_elements(FWHM))
  outx[0] = 0.0
  outy[0] = 0.0

;---------------------------------------------------------------------------
;  Setup file sizes, precision, and midpoints
;---------------------------------------------------------------------------

; This should be bigger than anything you need

  ax = -1024.                   ;relative reference
  ay = -1024.                   ;relative reference
  dx = 1.0
  dy = 1.0
  nx = 4000                     ; dimension of the histo image
  ny = 2000                     ; (maximum offset)
  nstep = 1000
  xedge=2040                    ;this should be approx the size of your image
  yedge=4600

;---------------------------------------------------------------------------
;  Run Source Extractor
;---------------------------------------------------------------------------

;prepare first file
  splog, 'Extract ', files[0]
  fits=mrdfits(files[0],/silent)
  fits=reform(fits[chipnum,*,*])
  mwrfits, fits, 'this_frame_tmpBJHL_1.fits', /create, /silent
  undefine, fits
  

;Run on first file
  input = 'sex this_frame_tmpBJHL_1.fits -SEEING_FWHM '+ FWHM_str[0] + ' -PIXEL_SCALE ' + pix
  spawn, input, errResult, EXIT_STATUS=exitStatus
  if exitStatus ne 0 then begin
     splog, 'Error! SExtractor had a problem'
     return
  endif
  
  FILE_DELETE, 'this_frame_tmpBJHL_1.fits'

;clean from junk. Read only first 200 brightest objects
;sort by mag
  spawn, 'sort -nk4 tmp.dat -o tmp.dat'


; Read in the first file data
  readcol, 'tmp.dat', x1, y1, a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, xg1, yg1, format='F,F,I,F,F,F,F,F,F,F,F,I,I,F,F,F,F', /silent, numline=200
  
  noedg=where(xg1 lt xedge and yg1 lt yedge, ned)
  if(ned gt 0) then begin
     xg1=xg1[noedg]
     yg1=yg1[noedg]
     x1=x1[noedg]
     y1=y1[noedg]
  endif

; Run on the remaining files

  for idz=1, n_elements(FWHM)-1 do begin

   ;prepare file  
     splog, 'Extract ', files[idz]
     fits=mrdfits(files[idz],/silent)
     fits=reform(fits[chipnum,*,*])
     mwrfits, fits, 'this_frame_tmpBJHL.fits', /create, /silent
     undefine, fits
     
     input = 'sex this_frame_tmpBJHL.fits -SEEING_FWHM '+ FWHM_str[idz] + ' -PIXEL_SCALE ' + pix
   
     spawn, input, errResult, EXIT_STATUS=exitStatus
     if exitStatus ne 0 then begin
        splog,  'Error! SExtractor had a problem'
        return
     endif
     
     FILE_DELETE, 'this_frame_tmpBJHL.fits'
    
    ;sort by mag
     spawn, 'sort -nk4 tmp.dat -o tmp.dat'

    ; Read in next files data
     readcol, 'tmp.dat', x2, y2, a2, b2, c2, d2, e2, f2, g2, h2, i2, j2, k2, l2, m2, xg2,yg2, format='F,F,I,F,F,F,F,F,F,F,F,I,I,F,F, F, F', /silent, numline=200
     
     noedg=where(xg2 lt xedge and yg2 lt yedge and xg2 gt 15 and  yg2 gt 15, ned)
     if(ned gt 0) then begin
        xg2=xg2[noedg]
        yg2=yg2[noedg]
        x2=x2[noedg]
        y2=y2[noedg]
     endif
     
    ;plot, xg2, yg2, psym=4
    ;oplot, xg1, yg1, psym=1
    ;stop


;---------------------------------------------------------------------------
;  Make pairs and calculate integer offset
;---------------------------------------------------------------------------

; array holding hit cts for good offset
     z = fltarr(nx,ny)
     z[*,*] = 0
; number of points in two data files
     np1 = n_elements(xg1)
     np2 = n_elements(xg2)
     
     for idx=0, np1-1 do begin
        for idy=0, np2-1 do begin
           zx = (xg1[idx]-xg2[idy])
           zy = (yg1[idx]-yg2[idy])
           
           k = round((zx-ax)/dx)
           l = round((zy-ay)/dy)
           
           if (k ge 1) and (k lt nx) and (l ge 1) and (l lt ny) then begin
              z[k,l]=z[k,l]+1
           endif
        endfor
     endfor



;---------------------------------------------------------------------------
;  Find  best integer offset zmax
;---------------------------------------------------------------------------


     zmax = max(z, zloc)
     indx = array_indices(z, zloc)
     imax = indx[0]
     jmax = indx[1]
     
; integer offset

     offx = ax+imax*dx
     offy = ay+jmax*dy
     
     
     if keyword_set(nozero) then begin
        if(offx lt 1 and offy lt 1) then begin
           splog, 'Shift is zero. Check if others are allowed...'
     ;I the shift is at zero, try to see if there are others 
     ;median reject single pixel and fwhm smooth.
           zsmo=filter_image(z,/median,fwhm=20)
           
           ;get new offset
           zmax = max(zsmo, zloc)
           indx = array_indices(zsmo, zloc)
           imax = indx[0]
           jmax = indx[1]
           
           if(imax eq 0 or jmax eq 0) then begin
              splog, 'I cannot found other offset. Set to zero '
              offx=0.
              offy=0.
              real=0
           endif else begin
              offx = ax+imax*dx
              offy = ay+jmax*dy
              splog, 'I found this new offset ', offx, offy
           endelse
           
        endif 
     endif
     

     if not keyword_set(real) then begin
        splog, 'Interger offsets for :' + files[idz]+ ' '+strtrim(string(offx),2)+ ' ' +strtrim(string(offy),2)
        outx[idz] = offx
        outy[idz] = offy
     endif
     
;---------------------------------------------------------------------------
;  Floating point offset (REAL)
;---------------------------------------------------------------------------

     if keyword_set(real) then begin

;---------------------------------------------------------------------------
;  Find all pairs within 1 pixel
;---------------------------------------------------------------------------

; nsv is the number of matched pairs
        nsv = 0
        sv1 = fltarr(np1*np2)
        sv2 = sv1
        
        match1 = intarr(np1*np2)
        match2 = match1
        
        for idx=0, np1-1 do begin
           for idy=0, np2-1 do begin
              zx = xg1[idx]-xg2[idy]-offx
              zy = yg1[idx]-yg2[idy]-offy
              if abs(zx) lt 2.0 and abs(zy) lt 2.0 then begin
                 nsv = nsv+1
                 sv1[nsv-1] = zx
                 sv2[nsv-1] = zy
                 
; lets match extra pairs 
                 match1[nsv-1] = idx
                 match2[nsv-1] =idy
                 
              endif
           endfor
        endfor
        
           
        sv1 = sv1[0:nsv-1]
        sv2 = sv2[0:nsv-1]

; average sv1 and sv2 - average offsets

        dum = median(sv1)
        zmax = median(sv2)
        
; loop on displacement and calculate sum di

        zmax = 99999999.0
        
        for k=0, nstep-1 do begin
           deltx = -1 +2*float(k)/float(nstep)
           for l=0, nstep-1 do begin
              delty = -1 +2*float(l)/float(nstep)
              
              dum = total(sqrt( (sv1-deltx)^2. + (sv2-delty)^2.))
              
              if(dum lt zmax) then begin
                 zmax = dum
                 svx = deltx
                 svy = delty
              endif
           endfor
        endfor
        

; The sign needs to be this way!
        roffx = offx+svx
        roffy = offy+svy
        
        
;print, svx, svy
;print, roffx, roffy
;print, medoffx, medoffy

; matched pairs are below
;for idx=0, nsv-1 do begin
;    x1[idx] = xg1[match1[idx]]
;    x2[idx] = xg2[match2[idx]]
;    y1[idx] = yg1[match1[idx]]
;    y2[idx] = yg2[match2[idx]]
;endfor

; correct matched pairs with offsets, and look at residuals

;dist1 = total(sqrt((x1-x2-roffx)^2 + (y1-y2-roffy)^2))
;dist2 = total(sqrt((x1-x2-medoffx)^2 + (y1-y2-medoffy)^2))

;print, dist1, dist2


; ----------

;print, roffx, roffy
        splog, string('Real offsets for: '+files[idz]+' '+strtrim(string(roffx),2)+' '+strtrim(string(roffy),2))
        
        outx[idz] = roffx
        outy[idz] = roffy
        
     endif
        
;stop
        
  endfor                        ;loop pver all the images


;forprint, outx, outy, /nocomment, /silent, format='(F9.4,4x,F9.4)', TEXTOUT=listfile+'_shift'
  FILE_DELETE, 'tmp.dat'
  FILE_DELETE, 'check.fits'
  
  
;stop
  
end
