;+ 
; NAME:
; mike_fitgap
;     Version 0.1
;
; PURPOSE:
;    Returns a 2d model fit of the flux associated with the inter-order
;      gaps.  The goal is to remove scattered light with a poor-man's
;      2d bspline.
;
; CALLING SEQUENCE:
;   
;   modelfit  = mike_fitgap(image, invvar, maskimage, $
;           [nxbkpt=nxbkpt, nybkpt=nybkpt, maxdev=maxdev])
;
;
; INPUTS:
;   image    - image to fit inter-order flux
;   invvar   - associated weights
;   maskimage- integer mask returned by mike_ordermask
;
; RETURNS:
;   modelfit - result of smooth 2d fit to image inter-order gap 
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   nxbkpt   - Number of breakpoints along columns (default 15)
;   nybkpt   - Number of breakpoints along rows    (default 10)
;   maxdev   - Used in rejection, pixels with flux-model residuals 
;                 greater than this are rejected (default 1000)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
;   The fit towards the edges can be improved if partial orders are
;     identified.  Without identifying partial orders, fitgap will try
;     to fit the flux in the partial orders, not helpful at all.
;   So we should add partial orders to ordr_str, but give them a flag
;
;   The current number of default breakpoints is not set based on binning
;   But should be about right for most frames
;
; EXAMPLES:
;   ncol = (size(image))[1]
;   nrow = (size(image))[2]
;   maskimage      = mike_ordermask(ncol, nrow, ordr_str, trim=1.5)
;   scatteredlight = mike_fitgap(tflat, tflativar, maskimage)
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   15-Jul-2003 Checked in by SB
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
       
function mike_fitgap, tflat, tflativar, maskimage, $
           nxbkpt=nxbkpt, nybkpt=nybkpt, maxdev=maxdev

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'fit = mike_fitgap(tflat, tflativar, maskimage, NXBKPT=, ' + $
        'NYBKPT=, MAXDEV= [v1.1]'
      return,-1
  endif


;      Number of breakpoints for row fitting
  if NOT keyword_set(nxbkpt) then nxbkpt = 5L

;      Number of breakpoints for column fitting
  if NOT keyword_set(nybkpt) then nybkpt = 10L

;
;     Maxdev is the maximum allowed deviation from the fit along columns
;
  if NOT keyword_set(maxdev) then maxdev=  1000.0


;
;     This is a strange setting to give  double the bkpt spacing in
;      the x direction at the left and right edges of the CCD.
;
  sz = size(tflat, /dimen)
  xbkpt = (findgen(nxbkpt)+1)/(nxbkpt+1) * sz[0]

;
;     Some simple initializations
;
  xset = 0
  yset = 0
  yp = -1
  x_base = findgen(sz[0])
  yp = findgen(sz[1])

;
;     Some not so simple initializations
;
  nord = 4L
  colmask= total(maskimage LT 0,2)
  x_here = where(smooth(colmask,15) GT 50,nhere)
  everyn = (nhere)/nxbkpt-1
  x_here = [0,x_here[1:nhere-2],sz[0]-1]
  tempfullbkpt = bspline_bkpts(x_here, everyn=everyn, nord=nord, /silent)
;     tempfullbkpt = bspline_bkpts(x_base, bkpt=xbkpt, nord=nord, /silent)
  tempset = create_bsplineset(tempfullbkpt, nord)
  xset = replicate(tempset, sz[1])

;
;      Perfrom a simple bspline fit to each row, only including pixels
;        in the order gaps.
;      

     for i=0,sz[1]-1 do begin
         inthegap = where((maskimage[*,i] LT 0) AND (tflativar[*,i] GT 0), n)
         if n GT 2.0*nxbkpt then begin
             fullbkpt = tempfullbkpt
             temp = bspline_iterfit(float(inthegap), tflat[inthegap,i], $
                                    invvar=tflativar[inthegap,i], $
                                    fullbkpt=fullbkpt, $
                                    nord=nord, /groupbadpix, maxrej=3, $
                                    upper=5, lower=5)
             if (size(temp,/tname) NE "INT") then xset[i] = temp
         endif
     endfor

;
;     Find a smooth evaluation of the bspline coefficients along columns
;       by, whatelse, bspline fit of the coefficients
;
  xsmoothset = xset
  for i=0, n_elements(xset[0].coeff) - 1 do begin  
      coeffivar = xset.icoeff[i]^2 * (xset.icoeff[i] GT 0.0005)
      fitset = bspline_iterfit(yp,  xset.coeff[i], yfit=tempfit, $
                               invvar=coeffivar, upper=4, lower=4, $
                               maxdev=maxdev, $
                               everyn=sz[1]/nybkpt, /silent) 
      xsmoothset.coeff[i] = tempfit
      xsmoothset.bkmask = 1
  endfor

;
;   Close your eyes here
;   This is just repeating bspline_valu to work on a full-2d image,
;    it's fast!  Action is called only once and looping is over bkpts only.
;
  gap_image = tflat * 0.0
  action = bspline_action(findgen(sz[0]), xsmoothset[0], $
                          upper=upper, lower=lower)
  spot = lindgen(nord)
  
  for i=0, n_elements(xsmoothset[0].fullbkpt) - 2*nord do begin
      ict = upper[i] - lower[i] + 1
      if (ict GT 0) then begin
          gap_image[lower[i]:upper[i],*] = action[lower[i]:upper[i],*] # $
            (xsmoothset.coeff)[i+spot,*]
      endif
  endfor 
     
  return, gap_image
end
