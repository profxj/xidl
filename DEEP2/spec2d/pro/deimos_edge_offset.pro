;+
; NAME:
;   deimos_edge_offset
;
; PURPOSE:
;   Compute global spatial offset caused by flexure
;
; CALLING SEQUENCE:
;   offs = deimos_edge_offset(xpos1, xpos2, ypos1, ypos2, specimage, $
;                  specivar[,/SUM])
;
; INPUTS:
;   xpos1,2    - xposition of slit edges (spatial direction)
;   ypos1,2    - yposition of slit edges (spectral direction)
;   specimage  - 2d spectrum
;   specivar   - 2d spectrum inverse variance
;
; OUTPUTS:
;   offs       - global spatial offset [pix]
;
; OPTIONAL KEYWORDS:
;   /SUM        - set to use the sum, rather than median, over the
;                wavelength direction for the cross-correlation shift
;   /LONGSLITS       - treat as long slit data (median & total)
; 
; REVISION HISTORY:
;   Written by Finkbeiner, 2002-Jun-05
;-
;------------------------------------------------------------------------------

function deimos_edge_offset, xpos1, xpos2, ypos1, ypos2, specimage, $
                             specivar, sum=sum, longslits=longslits

  if NOT keyword_set(xpos1) then message, 'xpos1, ypos1 must be set!'
  
  xposdim = size(xpos1, /dim)
  sizeimage=size(specimage,/dim)

  if n_elements(sum) eq 0 then sum=0
  if n_elements(longslits) eq 0 then longslits=0

  if n_elements(xposdim) gt 1 then nslits = xposdim[1] else nslits = 1
  dx = fltarr(nslits)
  mid=n_elements(xposdim)/nslits/2
  isgood=dx*0+1
  nx=sizeimage[0]
  ny=sizeimage[1]

  fakeprof=fltarr(nx)
  

;  vprint, 2,'Finding global flexure offset for this chip'


  for i=0, nslits-1 do begin 
     x0 = (xpos1[ny/2, i] < xpos2[ny/2, i] ) < nx-1
     x1 = (xpos2[ny/2, i] > xpos1[ny/2, i] ) > 0

     if (x0+2) gt (x1-2) then begin
         x0=x0-2
         x1=x1+2
     endif

     fakeprof[round(x0+2) > 0: round(x1-2) < (nx-1)] = 1
  endfor


  if longslits gt 0 then begin
      vprint,2,'Treating as longslit data!'
      tempimage=djs_median(specimage,width=5,boundary='reflect')
      realprofile=total(tempimage[*,ny/2-800:ny/2+800],2)
      delvarx,tempimage
  endif else begin
      if sum then realprofile=total(specimage[*,ny/2-500:ny/2+500],2) $
        else realprofile=djs_median(specimage[*,ny/2-100:ny/2+100],2)
  endelse

  maxlag=40
  lags=findgen(2*maxlag+1)-maxlag

; make sure alignment boxes don't pull us too much
  if sum eq 0 then realprofile=realprofile < 1.25*djs_median(realprofile, $
                                           width=501,boundary='reflect')

; do cross-correlation, subtract off BG to clean up 
  cc=c_correlate(fakeprof,realprofile,lags)
  cc=cc-djs_median(cc,width=15,boundary='reflect')
  maxcc=max(cc[5:2*maxlag-4],bestidx)
  shift=lags[bestidx+5]
;  vprint,2,'whole pixel shift: ',shift

  npoly = 3
  mintofit=(bestidx+5) - 2 > 0
  maxtofit=(bestidx+5) + 2 < (n_elements(lags) -1)

  shiftfit=poly_fit(lags[mintofit:maxtofit], $
                               cc[mintofit:maxtofit], npoly-1, $
                               sigma=err)


  dx0= - 0.5*shiftfit[1]/shiftfit[2]
  shifterr = abs(dx0)*sqrt((err[1]/shiftfit[1])^2 +$
                                            (err[2]/shiftfit[2])^2)
  
  vprint,2,'Global cross-correlation shift: ',dx0,' +/- ',shifterr
  
  dx=dx+dx0

  for i=0, nslits-1 do begin 
     x0 = xpos1[*, i] 
     x1 = xpos2[*, i] 
     length=abs(x1[mid]-x0[mid])
     if length lt 200 then begin 
        rect_spec = deimos_rectify_slit(specimage, specivar, x0+dx0, x1+dx0, $
                               /interp, xshift=xshift, npad=2, /recen) 

         dx[i] = xshift + dx0
     endif else isgood[i]=0

;     print, i, dx[i]

  endfor

  whgood=where(isgood AND xpos1[mid,*] gt 0 AND xpos2[mid,*] lt nx-1,goodct)

  if goodct lt 3 then whgood=where(isgood,goodct)
  if goodct lt 3 then whgood=findgen(nslits)

  offs = djs_median(dx[whgood])
  
;  dfpsplot,'ps/offset.ps',/square
;	plot,(xpos1[mid,*]+xpos2[mid,*])/2.,dx
;  dfpsclose

  vprint, 2,'Median centroid offset in spatial direction; ', offs

;  if abs(offs) lt stdev(dx[whgood])/sqrt(n_elements(whgood)) then begin
;      vprint,2,'Offset not significant (<1 sigma); using 0 instead!'
;      offs=0.
;  endif


  if finite(offs) eq 0 then message, 'failed!'
  return, offs
end










