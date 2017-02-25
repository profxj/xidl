;+
; NAME:
;   gmos_trnxy1to3
;
; PURPOSE:
;
;   Transforms a single vector of x y coordinates into 3-ccd image
;   coordinates. 
;
; CALLING SEQUENCE:
;  
;
; INPUTS:
;  transimg -- Transformed image. 
;
; OPTIONAL INPUTS:
; /NOFLEX  -- Do not apply flexure correction [necessary if your setup
;             has not been calibrated.  Contact JH or JXP for help if
;             this is the case.]
;  HAND_FWHM -- Set the FWHM of the object profile to this value (in
;               pixels)
; /NOHELIO -- Do not correct to heliocentric velocities
; /NOZAP   -- Do not flag CRs
;
; OUTPUTS:
;   outimg = mosaiced 3-ccd image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   
; PROCEDURES CALLED:
;   
; REVISION HISTORY:
;   11-Mar-2005  Written by JH + SB
;-  
;  -----------------------------------------------------------------------------
pro gmos_getradec, n, ra, dec, dra, ddec, dx, dy, kx, ky, xf, yf
  
  ra = ra-dra
  dec = dec-ddec
  
  xf = ra-ra
  yf = xf
  for i = 0, n do begin
     for j = 0, n do begin
        xf = xf+kx(i, j)*ra^j*dec^i
        yf = yf+ky(i, j)*ra^j*dec^i
     endfor
  endfor
  
  xf = xf+dx
  yf = yf+dy
  ra = ra+dra
  dec = dec+ddec
  
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro gmos_tranny, n, ra, dec, x, y, dra, ddec, dx, dy, kx, ky
  
  dx = total(x)/n_elements(x)
  dy = total(y)/n_elements(x)
  dra = total(ra)/n_elements(ra)
  ddec = total(dec)/n_elements(dec)
  
  x = x-dx
  y = y-dy
  ra = ra-dra
  dec = dec-ddec
  
  polywarp, x, y, ra, dec, n, kx, ky
  
  xf = x-x
  yf = xf
  for i = 0, n do begin
     for j = 0, n do begin
        xf = xf+kx(i, j)*ra^j*dec^i
        yf = yf+ky(i, j)*ra^j*dec^i
     endfor
  endfor
  
  mx = moment(x-xf)
  mxc = 2.5*sqrt(mx(1))+1e-8
  my = moment(y-yf)
  myc = 2.5*sqrt(my(1))+1e-8
  
  for Z = 0, 3 do begin
     qui = where(abs(x-xf) lt mxc and abs(y-yf) lt myc)
     if(n_elements(qui) ge (n+1)^2)then begin
        polywarp, x(qui), y(qui), ra(qui), dec(qui), n, kx, ky 
        
        xf = x-x
        yf = xf
        for i = 0, n do begin
           for j = 0, n do begin
              xf = xf+kx(i, j)*ra^j*dec^i
              yf = yf+ky(i, j)*ra^j*dec^i
           endfor
        endfor
        
     endif else begin
        polywarp, x(qui), y(qui), ra(qui), dec(qui), n-1, kx, ky
        xf = x-x
        yf = xf
        for i = 0, n-1 do begin
           for j = 0, n-1 do begin
              xf = xf+kx(i, j)*ra^j*dec^i
              yf = yf+ky(i, j)*ra^j*dec^i
           endfor
        endfor
     endelse
     
     mx = moment(x(qui)-xf(qui))
     mxc = 2.5*sqrt(mx(1))+1e-8
     my = moment(y(qui)-yf(qui))
     myc = 2.5*sqrt(my(1))+1e-8
  endfor
  
  dist = sqrt((x-xf)^2+(y-yf)^2)
  
  x = x+dx
  y = y+dy
  ra = ra+dra
  dec = dec+ddec
  RETURN
end


PRO GMOS_TRNXY1TO3, xin, yin, x3 = x3, y3 = y3, BTAG=BTAG 

  npt = n_elements(xin)
  x3 = fltarr(npt)
  y3 = fltarr(npt)
  
  IF NOT KEYWORD_SET(BTAG) THEN message, 'Must specify binning'

  bin = long(strsplit(btag,'*x*',/extract))
  specbin = bin[0]
  spatbin = bin[1]
  
  ygap1 = 36L/SPECBIN
  ygap2 = 36L/SPECBIN
  nx = 4608L/SPATBIN
  nmy = 3L*(2048L/SPECBIN) + ygap1 + ygap2

  ny = (nmy-ygap1-ygap2)/3
  fny = float(ny)
  
  ;; Convert my x,y convention to Gemini x,y (transpose and reverse)
  x = float(nmy) - yin
  y = xin
  
  IF SPECBIN EQ 2 THEN BEGIN
     c1ind = WHERE(x LE 1033, nc1)
     c2ind = WHERE(x GT 1033 AND X LE 2075, nc2)
     c3ind = WHERE(x GT 2075, nc3)
  ENDIF ELSE IF SPECBIN EQ 1 THEN BEGIN
     c1ind = WHERE(x LE 2066, nc1)
     c2ind = WHERE(x GT 2066 AND X LE 4151, nc2)
     c3ind = WHERE(x GE 4151, nc3)
  ENDIF 
  noscan = 64.0d/double(SPECBIN)
  ;; transform ccd1
  path = getenv('LONGSLIT_DIR') + '/calib/transform/GMOS/'
  ;; determine transformation
  rdfloat, path + 'chip1-'+btag+'.coo', tx, ty, ux, uy, /double, /silent
  gmos_tranny, 1, tx, ty, ux, uy, odra, oddec, odx, ody, kx, ky
  IF nc1 GT 0 THEN BEGIN
     gmos_getradec, 1, double(x[c1ind]), double(y[c1ind]) $
                    , odra, oddec, odx, ody, kx, ky, xf, yf
     ;;xf = (xf > 0.0) < float(ny)
     ;;yf = (yf > 0.0) < float(nx)
     x3[c1ind] = yf 
     y3[c1ind] = (fny - (xf - noscan)) + 2.0*fny + float(ygap1) + float(ygap2)
  ENDIF
  ;; transform ccd2
  ;; determine transformation
  rdfloat, path + 'chip2-'+btag+'.coo', tx, ty, ux, uy, /double, /silent
  gmos_tranny, 1, tx, ty, ux, uy, odra, oddec, odx, ody, kx, ky
  IF nc2 GT 0 THEN BEGIN
     gmos_getradec, 1, double(x[c2ind]), double(y[c2ind]) $
                    , odra, oddec, odx, ody, kx, ky, xf, yf
     ;;xf = (xf > 0.0) < float(ny)
     ;;yf = (yf > 0.0) < float(nx)
     x3[c2ind] = yf
     y3[c2ind] = (fny - (xf-noscan)) + fny + float(ygap1)
  ENDIF
  ;; no oscan subtraction for chip 3, as oscan is on other side
  ;; transform ccd3
  ;; determine transformation
  rdfloat, path + 'chip3-'+btag+'.coo', tx, ty, ux, uy, /double, /silent
  gmos_tranny, 1, tx, ty, ux, uy, odra, oddec, odx, ody, kx, ky
  IF nc3 GT 0 THEN BEGIN
     gmos_getradec, 1, double(x[c3ind]), double(y[c3ind]) $
                    , odra, oddec, odx, ody, kx, ky, xf, yf
     ;;xf = (xf > 0.0) < float(ny)
     ;;yf = (yf > 0.0) < float(nx)
     x3[c3ind] = yf
     y3[c3ind] = fny - xf
  ENDIF
  
  
  RETURN
END









