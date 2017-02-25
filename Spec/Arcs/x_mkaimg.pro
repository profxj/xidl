;+ 
; NAME:
; x_mkaimg
;     Version 1.1
;
; PURPOSE:
;   Given the 2D solution for the slope of the lines as a function of
;   position this code creates a wavelength image (i.e. assigns a
;   unique wavelength to each pixel in each order).  The xoffset is 
;   input in order to properly determine the edges of each order.  A
;   simple spline interpolation is used to determine the values.
;
;  Note, you should have shifted the order structure as necessary
;  prior to calling this routine.
;
; CALLING SEQUENCE:
; x_mkaimg, arc_fil, ordr_str, arc2d_fil, fil_fittrc, $
;                  out_fil, /CHK, /CLOBBER, BAD_ANLY=
;
; INPUTS:
;  arc_fil  -- Name of arc file
;  ordr_str -- Order strucure describing the echelle footprint
;  arc2d_fil -- Name of 2D wavelength solution
;  fil_fittrc --  Name of FITS file containing the 2D fit to the
;                tilted arc lines [No longer used!!]
;
; RETURNS:
;
; OUTPUTS:
;  2D wavelength image with name like 'Arcs/Arc_mb0439I.fits'
;
; OPTIONAL KEYWORDS:
;   /CLOBBER  - Overwrite previous image
;   /CHK      - Display the final image
;   /ZERO_SLOPE -- Do not tilt the arc lines [kludge!!]
;
; OPTIONAL OUTPUTS:
;  BAD_ANLY=  - Set to 1 if the code finds double valued solutions.
;
; COMMENTS:
;
; EXAMPLES:
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   15-May-2003 Written by SB
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_mkaimg, arc_fil, ordr_str, arc2d_fil, fil_fittrc, $
                   out_fil, CHK=chk, CLOBBER=clobber, NOEXTRA=noextra, $
                   BAD_ANLY=bad_anly, ZERO_SLOPE=zero_slope
;

  if  N_params() LT 5  then begin 
      print,'Syntax - ' + $
        'rslt = x_mkaimg( arc_fil, ordr_str, arc2d_fil, fil_fittrc, ' + $
        'out_fil, /CHK, ' + $
        '/CLOBBER, SHFTPRM= ) [v1.1]'
      return, -1
  endif 

  ;; Optional keywords
  bad_anly = 0
 
  ;; Order structure
  nordr = n_elements(ordr_str)

  ;; Grab Arc 2D Fit
  arc_2dfit = xmrdfits(arc2d_fil,1,/silent)

  ;; Check arc_fil for size
  head = xheadfits(arc_fil)
  sz = lonarr(2)
  sz[0] = sxpar(head, 'NAXIS1')
  sz[1] = sxpar(head, 'NAXIS2')
  ximage = lindgen(sz[0]) # replicate(1,sz[1])
  yimage = lindgen(sz[1]) ## replicate(1,sz[0])

  ;; Grab Slope fit
  ;mkaimg_str = xmrdfits(fil_fittrc,1,/silent)

  ;; Create Final image
  aimg = dblarr(sz[0],sz[1])

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; LOOPING

  ;; Create the ordermask image
  ordermask = x_ordermask(sz[0], sz[1], ordr_str, trim=0.1, noextra=NOEXTRA)

  print, 'x_mkaimg: Looping on Orders'
  bad_anly = 0 
  for q=0L,nordr-1 do begin
      print, 'x_mkaimg: Order', ordr_str[q].order

      if ordr_str[q].flg_anly EQ 0 then begin
         print, 'x_mkaimg: Skipping...'
         continue
      endif
      
      ;; Set order
      mkaimg_ordr = ordr_str[q].order
      
      ;; xcen
      xcen = (ordr_str[q].lhedg + ordr_str[q].rhedg)/2.
      
      inorder = where(ordermask EQ ordr_str[q].order, nin)
      if nin GT 300000. then stop
      
      ystart = 1.0d*yimage[inorder]
      xstart = (ximage - xcen ## replicate(1,sz[0]))[inorder]
      ycol = dindgen(sz[1])
      slope_spline = spl_init(ycol, double(ordr_str[q].arc_m))
      
;
;     Four iterations is plenty, converges quickly after the second call 
;
      ywave = ystart
      if not keyword_set(ZERO_SLOPE) then begin
         for islope=1,4 do begin
            slope = spl_interp(ycol, ordr_str[q].arc_m, slope_spline, ywave)
            ywave = ystart - slope*xstart 
            oldslope = slope
         endfor
      endif else begin
         if q EQ 0 then begin
            print, 'x_mkaimg: Not tilting the arc lines'
            print, 'x_mkaimg: I hope you know what you are doing!'
         endif
      endelse
              
      p = 2.0d * (ywave - arc_2dfit.nrm[0])/arc_2dfit.nrm[1]
      t = (2.0d *(ordr_str[q].order - arc_2dfit.nrmt[0]) $
           /arc_2dfit.nrmt[1])
      worky = flegendre(p, arc_2dfit.ny)
      workt = flegendre(t[0], arc_2dfit.no)
;      aimg[inorder] =  $
;        worky # transpose(reform(arc_2dfit.res, arc_2dfit.no, $
;                                 arc_2dfit.ny) ## workt)
      aimg[inorder] =  alog10($
        worky # transpose(reform(arc_2dfit.res, arc_2dfit.no, $
                                 arc_2dfit.ny) ## workt) / ordr_str[q].order)

      ys = sort(ywave)
      nss = n_elements(ys)
      multi_forw = total( aimg[inorder[ys[1:*]]] -  aimg[inorder[ys[0:nss-2]]] GT 0)
      multi_back = total( aimg[inorder[ys[1:*]]] -  aimg[inorder[ys[0:nss-2]]] LT 0)
      ;if keyword_set(chk) then print, 'x_mkaimg: mono?  ', multi_forw, multi_back
      if (multi_forw GT 0) AND (multi_back GT 0) then  begin
         print, 'x_mkaimg: This order is not single-valued in wavelength.'
         print, 'The 2D wavelength fit is likely bad!'
         print, 'Setting flg_anly to 0'
         stop
         bad_anly = 1
      endif 
 
  endfor

  if keyword_set( CHK ) then begin
      xatv, aimg, /block, min=3.5, max=4.0
  endif

  if bad_anly GT 0 then begin
    print, 'x_mkaimg: Bad wavelength solution, not writing to disk'
    return, out_fil
  endif

  ;; Output
  print, 'x_mkaimg_work: Writing: ', out_fil
  mwrfits, aimg, out_fil, /create
  print, 'x_mkaimg_work: Compressing..'
  spawn, 'gzip -f '+out_fil
  
  print, 'x_mkaimg: All done..'
          
  return, out_fil
end

